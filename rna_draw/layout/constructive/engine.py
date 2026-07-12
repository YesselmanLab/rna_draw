"""`ConstructiveEngine`: bottom-up, checker-clean-by-construction layout.

M2 scope (see the constructive-engine runbook, milestone M2): ANY
pseudoknot-free structure up to `_MAX_NUCLEOTIDES` nt (the pure-Python
engine is too slow past that -- routes to the fallback engine until a
future C++ port, runbook Risk R2). A general post-order recursion over
`rna_draw.layout.structure_tree.StructureTree`: stems (collapsed through
runs of stacked pairs, `structure_tree.collapse_stem`) plus the loop each
one closes -- a hairpin, a bulge/interior loop, or a multiloop of ANY
degree and ANY nesting depth, all packed by the same exact circular
angular-interval packer (`geometry_helpers.pack_loop_angles`). The one
loop with an OPEN boundary -- the exterior (dangling 5'/3' tails + any
number of top-level branches, no closing pair) -- is packed along a line
instead (`geometry_helpers.pack_line_positions`); see `_place_exterior`.

The change from M1 that makes this general composition sound: a child
branch's slot on its parent's loop/line is no longer a fixed constant --
it is `envelope.branch_reach`, a topology-only bound (S5) on the child's
WHOLE subtree, computed bottom-up BEFORE any placement happens (as a side
effect of the top-level `envelope.branch_reach` calls in `_place_exterior`
-- see `envelope.py`'s module docstring for the disjointness proof this
composition relies on).

Every produced layout is verified against the frozen checker
(`rna_draw.overlap.check_overlaps`) before it is returned: a dirty result
is never handed back silently, it raises `EngineError` instead (see
`ConstructiveEngine.layout`). `EngineError` is also the ONLY exception this
module raises for a give-up path (never a bare `RuntimeError`), including
converting a pathologically deep structure's `RecursionError` (see
`_ensure_recursion_headroom`).
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass

from rna_draw.layout.base import EngineError, is_pseudoknot_free
from rna_draw.layout.structure_tree import (
    Branch,
    Loop,
    StructureTree,
    build_structure_tree,
    collapse_stem,
)
from rna_draw.overlap import OverlapParams, check_overlaps, rescale_coords
from rna_draw.parameters import DrawParameters
from rna_draw.render_rna import get_pairmap_from_secstruct

from . import compaction, envelope
from .geometry_helpers import (
    EXTERIOR_AXIS,
    Point,
    StemLadder,
    loop_member_point,
    midpoint,
    pack_bulge_linear,
    pack_line_positions,
    place_bulge_geometry,
    place_stem,
    stem_base_for_attachment,
)

# The pure-Python engine's recursive tree-DP is too slow (and grows Python
# stack depth with structural nesting) to be worth attempting past this
# length -- larger structures should route to a faster fallback engine
# instead (runbook Risk R2). Checked FIRST in `layout`, before any other
# work, so the engine never even starts on a giant structure.
#
# RAISED 1200 -> 4000 after the bulge-chain straight-continuation fix
# (`_place_bulge`/`envelope._bulge_branch_reach`) removed the dominant cause
# of slow/huge constructions on real rRNA: measured on the FULL
# `benchmarks/hard_set.json` 1200-4000nt bucket (100 structures, guard
# bypassed for the measurement) -- 89/100 checker-clean by construction, the
# other 11 a fast (<0.01s) `_MAX_REACH` bail on a DIFFERENT, still-circular
# cause (long chains of degree-2 multiloops) -- FIXED by
# `envelope.lateral_reach`/`_degree2_packing` (see `_MAX_REACH`'s docstring),
# 0 timeouts, 0 dirty, worst-case wall clock 7.3s (well under the benchmark
# harness's 30s per-structure kill). No structure past 4000nt has been
# measured -- do not raise further without new evidence.
_MAX_NUCLEOTIDES = 4000

# `_place_branch`/`_place_loop_members` recurse one Python stack frame per
# structural nesting level (not per nucleotide); this is a generous
# per-nucleotide upper bound on how deep that could plausibly go, used to
# raise the recursion limit defensively (never lowers an already-higher
# limit some other caller set).
_RECURSION_FRAMES_PER_NT = 4
_RECURSION_HEADROOM = 1000

# A hard cap on `envelope.branch_reach`, checked BEFORE any coordinate is
# written or the checker runs (see `_place_exterior`). ORIGINALLY this
# guarded a long CHAIN of single-child bulge/interior loops, which made
# `envelope.branch_reach` compound ~3x PER NESTING LEVEL via the circular
# envelope (radius >= dominant child's own reach, by `radius_floor`, then
# +2x that radius on top) -- FIXED by `_place_bulge`/
# `envelope._bulge_branch_reach`: a single-child loop is now placed as a
# straight continuation of the parent stem's axis instead of a fresh
# circular envelope, so a bulge CHAIN's reach grows LINEARLY (one constant
# term per level), not exponentially. Measured impact on
# `benchmarks/hard_set.json`'s <=1200nt buckets: clean-by-construction
# 37.7% -> 99.1% (347/350), the 3 residual failures NOT bulge chains (see
# below).
#
# The 3 residual <=1200nt failures (and 11 more in the 1200-4000nt bucket)
# were a DIFFERENT, still-circular cause: a long CHAIN of degree-2
# multiloops -- one big "continuing" branch plus one small side branch
# (e.g. a tiny hairpin), repeated many nesting levels deep in real rRNA --
# compounding `2 * packing.radius` per level the same way a bulge chain
# used to. FIXED by `envelope.lateral_reach`/`_degree2_packing`: a
# degree-2 loop's DOMINANT child is now (when it pays off; see
# `_degree2_packing`'s docstring) pinned as a collinear straight
# continuation too, sized by its directional `lateral_reach` instead of
# its isotropic `branch_reach`, so a degree-2 CHAIN's reach also grows
# LINEARLY. This guard remains as the never-silent backstop for anything
# still pathological (e.g. a giant 3+-way junction with two or more large
# children -- genuinely isotropic, out of scope for a straight-
# continuation fix) -- never removed, only its trigger rate driven down.
#
# THRESHOLD, calibrated against `benchmarks/hard_set.json` (bypassing this
# guard and timing construction+check directly): time scales with
# log(reach), staying well under a second up to ~2e6, ~2s at ~5e6, and
# crossing 10-19s by ~1.5e7-3e7 (now dominated by `_build_verified`'s
# checker re-runs at huge coordinate magnitudes, not construction itself,
# which stays cheap) -- `5e6` keeps every admitted structure comfortably
# fast (a couple of seconds, worst case) while still rejecting the
# multi-second tail before it ever starts. Raising it to ~3e7 (measured)
# would rescue the 3 residual <=1200nt failures above at ~11-13s each --
# not done here (bounded value for a real per-structure cost increase; see
# the handoff report for the full tradeoff).
_MAX_REACH = 5_000_000.0

# S8's bounded re-inflation schedule: the construction is clean BY PROOF for
# sibling-vs-sibling disjointness (`envelope.py`'s soundness argument) and
# for the two structural risks identified and fixed while building M2 --
# lopsided loops leaving their radius-fitting slack entirely on one side
# (`geometry_helpers._angles_at` now splits it), and an exterior branch's
# own "return" strand landing on the same side as an incoming backbone edge
# (`EXTERIOR_AXIS`'s sign). What's left is a DIFFERENT, more diffuse risk:
# a branch's own near-seam interior content sitting close enough to its
# OUTERMOST rung's far strand that the backbone edge continuing PAST that
# strand (to whatever comes next, outside this branch's own subtree) can
# graze it -- a genuine but small-margin near-miss (observed within ~2% of
# the checker's clearance), not a directional/structural flaw. Retrying
# with larger clearance margins is a legitimate, monotone way to close that
# gap (never accepted unless the checker actually confirms clean).
_MARGIN_SCALES: tuple[float, ...] = (1.0, 1.5, 2.0, 3.0, 5.0)


@dataclass
class _LayoutState:
    """Mutable accumulators threaded through the placement recursion.

    Args:
        tree: The structure tree being laid out.
        x: Nucleotide x-coordinates, filled in as each stem/loop is placed.
        y: Nucleotide y-coordinates, filled in as each stem/loop is placed.
        params: Target geometry (`NODE_R`/`PRIMARY_SPACE`/`PAIR_SPACE`).
        cache: S5 envelope memoization (`envelope.ReachCache`).
        margin_scale: Multiplier on every additive clearance margin for
            this attempt (see `_MARGIN_SCALES`).
    """

    tree: StructureTree
    x: list[float]
    y: list[float]
    params: DrawParameters
    cache: envelope.ReachCache
    margin_scale: float = 1.0


class ConstructiveEngine:
    """Bottom-up constructive layout, clean-by-construction (M2 scope)."""

    name = "constructive"

    def __init__(self, params: DrawParameters | None = None) -> None:
        """Store the drawing parameters this engine lays out with.

        Args:
            params: `NODE_R`/`PRIMARY_SPACE`/`PAIR_SPACE` to use; defaults
                to `DrawParameters()`.
        """
        self._params = params or DrawParameters()

    def layout(self, secstruct: str) -> tuple[list[float], list[float]]:
        """Lay out `secstruct`, guaranteed checker-clean or an `EngineError`.

        Args:
            secstruct: Dot-bracket secondary structure.

        Returns:
            `(x, y)`, one entry per nucleotide in sequence order, rescaled
            so the median backbone step equals `params.PRIMARY_SPACE` --
            unless rescaling would make an already checker-clean layout
            dirty (see `_rescale_or_keep_clean`), in which case the
            un-rescaled (still checker-clean) coordinates are returned.

        Raises:
            EngineError: If `secstruct` is longer than `_MAX_NUCLEOTIDES`,
                is a pseudoknot, or every attempt in `_MARGIN_SCALES`
                produced a checker-dirty layout (never returned silently).
        """
        n = len(secstruct)
        if n > _MAX_NUCLEOTIDES:
            raise EngineError(
                f"ConstructiveEngine declines structures over {_MAX_NUCLEOTIDES} nt "
                f"(got {n} nt); route to a faster fallback engine instead"
            )
        if n == 0:
            return [], []
        if n < 2:
            return [0.0] * n, [0.0] * n
        if not is_pseudoknot_free(secstruct):
            raise EngineError(f"ConstructiveEngine cannot lay out a pseudoknot: {secstruct!r}")

        pair_map = get_pairmap_from_secstruct(secstruct)
        tree = build_structure_tree(pair_map)
        _ensure_recursion_headroom(n)
        x, y = _build_verified(tree, pair_map, n, self._params, secstruct)
        return _rescale_or_keep_clean(x, y, pair_map, self._params.PRIMARY_SPACE)


def _build_verified(
    tree: StructureTree,
    pair_map: list[int],
    n: int,
    params: DrawParameters,
    secstruct: str,
) -> tuple[list[float], list[float]]:
    """Build `secstruct`'s layout, retrying at larger margins if needed (S8).

    Args:
        tree: The structure tree to lay out.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.
        n: The structure's length.
        params: Target geometry.
        secstruct: The original structure, for error messages.

    Returns:
        The first attempt's `(x, y)` that comes back checker-clean.

    Raises:
        EngineError: If a `RecursionError` fires, or every attempt in
            `_MARGIN_SCALES` is dirty (never returned silently).
    """
    last_overlaps = 0
    for margin_scale in _MARGIN_SCALES:
        state = _LayoutState(
            tree=tree,
            x=[0.0] * n,
            y=[0.0] * n,
            params=params,
            cache=envelope.ReachCache(),
            margin_scale=margin_scale,
        )
        try:
            _place_exterior(state)
        except RecursionError as exc:
            raise EngineError(
                f"ConstructiveEngine hit Python's recursion limit on {secstruct!r} "
                "-- refusing to crash; route to a fallback engine instead"
            ) from exc
        report = check_overlaps(state.x, state.y, pair_map, OverlapParams())
        if report.passed:
            return _compact_or_keep(tree, pair_map, state)
        last_overlaps = report.num_overlaps
    raise EngineError(
        f"ConstructiveEngine produced a dirty layout for {secstruct!r} after "
        f"{len(_MARGIN_SCALES)} attempt(s) ({last_overlaps} overlaps on the last) "
        "-- refusing to return it silently"
    )


def _compact_or_keep(
    tree: StructureTree, pair_map: list[int], state: _LayoutState
) -> tuple[list[float], list[float]]:
    """Try the M3 checker-gated compaction pass; keep it only if still clean.

    `state.x`/`state.y` are already checker-clean (the "sound" layout).
    `compaction.compact_layout` only ever applies a move it has itself
    locally verified against the frozen checker, so this whole-structure
    re-check is belt-and-suspenders, not the primary safety mechanism --
    but it is what makes the never-silent-overlap contract airtight even if
    a local check's scope assumption ever turns out to be wrong.

    Args:
        tree: The structure tree the sound layout was built from.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.
        state: The just-verified sound `_LayoutState` (its `cache` holds
            the sound `loop_packing` the compaction pass compares against).

    Returns:
        The compacted `(x, y)` if it stays checker-clean, else the
        pre-compaction sound `(x, y)`.
    """
    compact_x, compact_y = compaction.compact_layout(
        tree,
        state.x,
        state.y,
        pair_map,
        state.params,
        OverlapParams(),
        state.cache,
        state.margin_scale,
    )
    if check_overlaps(compact_x, compact_y, pair_map, OverlapParams()).passed:
        return compact_x, compact_y
    return state.x, state.y


def _ensure_recursion_headroom(n: int) -> None:
    """Raise Python's recursion limit if `n` residues could plausibly nest deep.

    Never lowers an already-higher limit some other caller set.

    Args:
        n: The structure's length.
    """
    needed = n * _RECURSION_FRAMES_PER_NT + _RECURSION_HEADROOM
    if sys.getrecursionlimit() < needed:
        sys.setrecursionlimit(needed)


def _rescale_or_keep_clean(
    x: list[float], y: list[float], pair_map: list[int], primary_space: float
) -> tuple[list[float], list[float]]:
    """Rescale to the target backbone step, but never at the cost of clean.

    `x`/`y` are already checker-clean (verified by `_build_verified`),
    built in absolute geometry units. `rescale_coords` uniformly scales by
    the *median* consecutive-nucleotide step to hit `PRIMARY_SPACE` exactly
    -- usually a near-identity touch-up, but for a structure whose backbone
    is mostly loop-to-loop hops rather than in-stem steps, that median can
    be skewed well away from `PRIMARY_SPACE`, and uniformly rescaling by it
    would shrink/grow the already-correct absolute clearances, breaking
    them (`OverlapParams` are absolute, not relative). So: rescale, but
    only keep it if it stays clean; otherwise fall back to the un-rescaled
    coordinates, which are already proven clean.

    Args:
        x: Nucleotide x-coordinates, already checker-clean.
        y: Nucleotide y-coordinates, already checker-clean.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.
        primary_space: Target median consecutive-nucleotide step.

    Returns:
        The rescaled `(x, y)` if that stays clean, else the original
        `(x, y)`.
    """
    rescaled_x, rescaled_y = rescale_coords(x, y, primary_space)
    if check_overlaps(rescaled_x, rescaled_y, pair_map, OverlapParams()).passed:
        return rescaled_x, rescaled_y
    return x, y


def _place_exterior(state: _LayoutState) -> None:
    """Place the exterior loop: dangling tails + any number of top-level
    branches, along an OPEN line (`pack_line_positions`) -- the one loop
    with no closing pair, so no reserved seam is needed (S7).

    Args:
        state: Layout accumulators; `state.x`/`state.y` are filled in place.
    """
    exterior = state.tree.exterior
    branch_by_start = {branch.start: branch for branch in exterior.children}
    extents = [_member_extent(state, member, branch_by_start) for member in exterior.members]
    _check_reach_bounded(extents)
    positions = pack_line_positions(extents, state.params.PRIMARY_SPACE)

    for member, x_pos in zip(exterior.members, positions):
        anchor = (x_pos, 0.0)
        branch = branch_by_start.get(member)
        if branch is None:
            state.x[member], state.y[member] = anchor
            continue
        _place_branch(state, branch, anchor, EXTERIOR_AXIS)


def _check_reach_bounded(extents: list[float]) -> None:
    """Fast, cheap guard against `_MAX_REACH`-scale envelope blow-up.

    Computing `envelope.branch_reach` for every top-level branch (already
    done by the time `extents` is built) is itself fast (pure arithmetic,
    no coordinates written, no checker run) -- this is the earliest point a
    pathological structure (see `_MAX_REACH`'s docstring) can be caught,
    well before the expensive part (placement, then the checker on
    resulting huge coordinates) would run.

    Args:
        extents: Every exterior member's required half-width, as built by
            `_place_exterior` (a top-level branch's own `envelope.branch_reach`,
            propagated bottom-up from the whole tree it contains -- so
            bounding these few values bounds every branch anywhere in the
            structure).

    Raises:
        EngineError: If any extent exceeds `_MAX_REACH`.
    """
    worst = max(extents, default=0.0)
    if worst > _MAX_REACH:
        raise EngineError(
            f"ConstructiveEngine envelope reach {worst:.3g} exceeds the "
            f"{_MAX_REACH:.3g} sanity cap (likely a giant isotropic (3+-way) "
            "junction compounding the envelope bound -- see _MAX_REACH's "
            "docstring); route to a fallback engine instead"
        )


def _member_extent(state: _LayoutState, member: int, branch_by_start: dict[int, Branch]) -> float:
    """The chord/linear half-width the exterior slot for `member` needs.

    Args:
        state: Layout accumulators.
        member: An exterior member index (a tail nt, or a top-level
            branch's own start index).
        branch_by_start: Maps a top-level branch's start index to itself.

    Returns:
        The branch's `envelope.branch_reach` if `member` starts one, else a
        bare unpaired nt's fixed `envelope.unpaired_half_width`.
    """
    branch = branch_by_start.get(member)
    if branch is None:
        return envelope.unpaired_half_width(state.params, state.margin_scale)
    return envelope.branch_reach(
        state.tree, branch.closing_pair, state.params, state.cache, state.margin_scale
    )


def _place_branch(state: _LayoutState, branch: Branch, attachment: Point, axis_dir: Point) -> None:
    """Place one branch: its straight (collapsed) stem, then the loop it
    closes -- a hairpin, a bulge/interior loop, or a multiloop of any
    degree, recursed into via `_place_loop_members`.

    Args:
        state: Layout accumulators.
        branch: The child branch to place.
        attachment: Where `branch.closing_pair[0]` must land exactly (a
            point on the parent loop's circle, or the exterior line) -- see
            `_place_collapsed_stem`.
        axis_dir: Unit vector pointing away from the parent (radially
            outward from a loop's center, or straight up from the exterior
            line).
    """
    loop, tip_a, tip_b = _place_collapsed_stem(state, branch.closing_pair, attachment, axis_dir)
    _place_loop_members(state, loop, tip_a, tip_b, axis_dir)


def _place_collapsed_stem(
    state: _LayoutState, closing_pair: tuple[int, int], attachment: Point, axis_dir: Point
) -> tuple[Loop, Point, Point]:
    """Collapse a run of stacked pairs from `closing_pair` and place it as
    one straight ladder, with `closing_pair[0]` landing EXACTLY at
    `attachment` (see `geometry_helpers.stem_base_for_attachment` for why
    that exact placement, not just "near" it, is load-bearing for S5's
    envelope soundness argument).

    Args:
        state: Layout accumulators; `state.x`/`state.y` are filled in place
            for every nucleotide in the collapsed stem.
        closing_pair: The outermost pair of the stem to place.
        attachment: Where `closing_pair[0]` (the branch's own attachment
            nucleotide) must land -- the point the parent packer assigned
            this branch's slot.
        axis_dir: Unit vector the stem rises along, base toward tip.

    Returns:
        `(loop, tip_a, tip_b)`: the loop this stem closes (not yet placed),
        and the innermost rung's two strand coordinates (that loop's own
        closing-pair coordinates).
    """
    depth, loop = collapse_stem(state.tree, closing_pair)
    base_point = stem_base_for_attachment(attachment, axis_dir, state.params.PAIR_SPACE)
    ladder = place_stem(
        depth, base_point, axis_dir, state.params.PRIMARY_SPACE, state.params.PAIR_SPACE
    )
    _write_ladder(state, ladder, closing_pair, depth)
    return loop, ladder.strand_a[-1], ladder.strand_b[-1]


def _write_ladder(
    state: _LayoutState, ladder: StemLadder, closing_pair: tuple[int, int], depth: int
) -> None:
    """Write a placed stem's coordinates into `state.x`/`state.y`.

    Args:
        state: Layout accumulators; mutated in place.
        ladder: The placed stem.
        closing_pair: `(i, j)`, the stem's outermost pair.
        depth: Number of rungs in `ladder`.
    """
    i, j = closing_pair
    for k in range(depth):
        state.x[i + k], state.y[i + k] = ladder.strand_a[k]
        state.x[j - k], state.y[j - k] = ladder.strand_b[k]


def _place_loop_members(
    state: _LayoutState, loop: Loop, tip_a: Point, tip_b: Point, axis_dir: Point
) -> None:
    """Place a loop's interior members (unpaired nts + child stem bases),
    dispatching on child count: a genuine multiloop (2+ children) or a
    terminal hairpin loop (0 children) packs them on a circle
    (`_place_circular_loop_members`); a bulge/interior loop (exactly 1
    child) places them as a straight continuation of the stem's own axis
    instead (`_place_bulge`) -- see `envelope.branch_reach`'s module
    docstring for why the circular case compounds `~3x` per nesting level
    on a long bulge chain while the straight case only adds a constant.

    Args:
        state: Layout accumulators; mutated in place.
        loop: The loop to place; `tip_a`/`tip_b` (its own closing pair) are
            already placed by the caller's stem.
        tip_a: The closing pair's near-strand coordinate.
        tip_b: The closing pair's far-strand coordinate.
        axis_dir: Unit vector pointing from the stem base toward this loop
            (the direction the loop opens away from its parent).
    """
    if len(loop.children) == 1:
        _place_bulge(state, loop, tip_a, tip_b, axis_dir)
    else:
        _place_circular_loop_members(state, loop, tip_a, tip_b, axis_dir)


def _place_circular_loop_members(
    state: _LayoutState, loop: Loop, tip_a: Point, tip_b: Point, axis_dir: Point
) -> None:
    """Place a multiloop's (2+ children) or hairpin's (0 children) interior
    members on a circle beyond its own closing-pair rung, reusing the
    packing already computed by S5's envelope pass (`envelope.loop_packing`)
    -- so the radius/angles here are guaranteed identical to the ones every
    child's `envelope.branch_reach` was measured against.

    Recurses into every child branch (`_place_branch`).

    Args:
        state: Layout accumulators; mutated in place.
        loop: The loop to place; `tip_a`/`tip_b` (its own closing pair) are
            already placed by the caller's stem.
        tip_a: The closing pair's near-strand coordinate.
        tip_b: The closing pair's far-strand coordinate.
        axis_dir: Unit vector pointing from the stem base toward this loop.
    """
    params = state.params
    interior = loop.members[1:-1]
    branch_by_start = {branch.start: branch for branch in loop.children}
    packing = envelope.loop_packing(state.tree, loop, params, state.cache, state.margin_scale)

    center = midpoint(tip_a, tip_b)
    loop_center = (
        center[0] + packing.radius * axis_dir[0],
        center[1] + packing.radius * axis_dir[1],
    )
    zero_dir = (-axis_dir[0], -axis_dir[1])

    for member, angle in zip(interior, packing.angles):
        anchor = loop_member_point(loop_center, zero_dir, packing.radius, angle)
        branch = branch_by_start.get(member)
        if branch is None:
            state.x[member], state.y[member] = anchor
            continue
        _place_branch(state, branch, anchor, _radial_dir(loop_center, anchor))


def _place_bulge(
    state: _LayoutState, loop: Loop, tip_a: Point, tip_b: Point, axis_dir: Point
) -> None:
    """Place a bulge/interior loop (exactly 1 child) as a straight
    continuation of the parent stem's own axis (`envelope.py`'s module
    docstring): the unpaired members form a short kink/offset to either
    side, and the child continues along the SAME `axis_dir`, so a whole
    CHAIN of such loops lays out as one straight run instead of a fresh
    circular envelope per level.

    Args:
        state: Layout accumulators; mutated in place.
        loop: The loop to place (`len(loop.children) == 1`); `tip_a`/`tip_b`
            (its own closing pair) are already placed by the caller's stem.
        tip_a: The closing pair's near-strand coordinate.
        tip_b: The closing pair's far-strand coordinate.
        axis_dir: Unit vector pointing from the stem base toward this loop.
    """
    child = loop.children[0]
    n_before, n_after = envelope.bulge_split(loop, child)
    packing = pack_bulge_linear(n_before, n_after)
    interior = loop.members[1:-1]
    unpaired_before, unpaired_after = interior[:n_before], interior[n_before + 1 :]

    near_points, far_points, attachment = place_bulge_geometry(
        tip_a, tip_b, axis_dir, packing, state.params.PRIMARY_SPACE, state.params.PAIR_SPACE
    )
    for member, point in zip(unpaired_before, near_points):
        state.x[member], state.y[member] = point
    for member, point in zip(unpaired_after, far_points):
        state.x[member], state.y[member] = point
    _place_branch(state, child, attachment, axis_dir)


def _radial_dir(center: Point, point: Point) -> Point:
    """Unit vector from `center` toward `point`.

    Args:
        center: The reference center (a loop's circle center).
        point: A point on that circle.

    Returns:
        The unit vector from `center` to `point`.

    Raises:
        EngineError: If `point` coincides with `center` (a degenerate,
        zero-radius loop -- should be unreachable given `pack_loop_angles`'
        positive `radius_floor`).
    """
    dx, dy = point[0] - center[0], point[1] - center[1]
    length = math.hypot(dx, dy)
    if length == 0.0:
        raise EngineError("degenerate loop: a member's anchor coincides with the loop center")
    return (dx / length, dy / length)


__all__ = ["ConstructiveEngine"]
