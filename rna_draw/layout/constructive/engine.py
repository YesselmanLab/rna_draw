"""`ConstructiveEngine`: bottom-up, checker-clean-by-construction layout.

M1 scope (see `.claude/plans/constructive-engine-runbook.md`): a single
multiloop of hairpins. Concretely, the engine only handles an exterior loop
with exactly one top-level branch (no dangling 5'/3' tails, no sibling
branches), whose closing stem -- collapsed through any run of stacked pairs
-- leads to one loop (a plain hairpin, or a multiloop whose every child
branch itself collapses straight to a plain hairpin). Every other shape
(pseudoknots, multi-branch/tailed exterior loops, bulges/interior loops,
nested multiloops) raises `EngineError`; those are later milestones (M2's
S5-S7).

Every produced layout is verified against the frozen checker
(`rna_draw.overlap.check_overlaps`) before it is returned: a dirty result
is never handed back silently, it raises `EngineError` instead (see
`ConstructiveEngine.layout`).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from rna_draw.layout.base import EngineError, is_pseudoknot_free
from rna_draw.layout.structure_tree import Branch, Loop, StructureTree, build_structure_tree
from rna_draw.overlap import OverlapParams, check_overlaps, rescale_coords
from rna_draw.parameters import DrawParameters
from rna_draw.render_rna import get_pairmap_from_secstruct

from .geometry_helpers import (
    Point,
    StemLadder,
    loop_member_point,
    midpoint,
    pack_loop_angles,
    place_stem,
)

_RADIUS_STEP = 1.0
_RADIUS_SEARCH_SLACK = 1.0

# Clearance margins added on top of each slot's provably-sufficient minimum
# chord half-width (see `_place_loop_members`'s docstring for the
# `chord >= hw_i + hw_j` proof) before angular packing: just enough to clear
# float noise and the small deviation a stem rung's actual disks have from
# their anchor point (offset `PAIR_SPACE / 2` tangentially). M1 favors
# correctness over maximal compactness -- see the runbook's M3 compactness
# pass for tightening this further.
_UNPAIRED_MARGIN = 0.5
_STEM_RUNG_MARGIN = 9.0


@dataclass
class _LayoutState:
    """Mutable accumulators threaded through the M1 placement recursion.

    Args:
        tree: The structure tree being laid out.
        x: Nucleotide x-coordinates, filled in as each stem/loop is placed.
        y: Nucleotide y-coordinates, filled in as each stem/loop is placed.
        params: Target geometry (`NODE_R`/`PRIMARY_SPACE`/`PAIR_SPACE`).
    """

    tree: StructureTree
    x: list[float]
    y: list[float]
    params: DrawParameters


class ConstructiveEngine:
    """Bottom-up constructive layout, clean-by-construction (M1 scope)."""

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
            EngineError: If `secstruct` is a pseudoknot, is not the M1
                single-multiloop-of-hairpins shape this engine supports, or
                the constructed layout fails the frozen checker (never
                returned silently).
        """
        n = len(secstruct)
        if n == 0:
            return [], []
        if n < 2:
            return [0.0] * n, [0.0] * n
        if not is_pseudoknot_free(secstruct):
            raise EngineError(f"ConstructiveEngine cannot lay out a pseudoknot: {secstruct!r}")

        pair_map = get_pairmap_from_secstruct(secstruct)
        state = _LayoutState(
            tree=build_structure_tree(pair_map),
            x=[0.0] * n,
            y=[0.0] * n,
            params=self._params,
        )
        _place_exterior(state)
        return _rescale_or_keep_clean(state, pair_map, secstruct)


def _rescale_or_keep_clean(
    state: _LayoutState, pair_map: list[int], secstruct: str
) -> tuple[list[float], list[float]]:
    """Rescale to the target backbone step, but never at the cost of clean.

    `place_stem`/`pack_loop_angles` already build every coordinate in
    `state.params`' absolute units, so `state.x`/`state.y` are checker-clean
    by construction (belt-and-suspenders-verified below) BEFORE any
    rescaling. `rescale_coords` uniformly scales by the *median*
    consecutive-nucleotide step to hit `PRIMARY_SPACE` exactly -- usually a
    near-identity touch-up, but for a structure whose backbone is mostly
    loop-to-loop hops rather than in-stem steps (e.g. many depth-1
    children), that median can be skewed well away from `PRIMARY_SPACE`,
    and uniformly rescaling by it would shrink/grow the already-correct
    absolute clearances, breaking them (`OverlapParams` are absolute, not
    relative). So: rescale, but only keep it if it stays clean; otherwise
    fall back to the un-rescaled coordinates, which are already proven
    clean.

    Args:
        state: Layout accumulators; `state.x`/`state.y` already fully placed.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.
        secstruct: The original structure, for the error message.

    Returns:
        The rescaled `(x, y)` if that stays clean, else the un-rescaled
        `(state.x, state.y)`.

    Raises:
        EngineError: If even the un-rescaled construction fails the frozen
            checker (never returned silently).
    """
    _verify_clean(state.x, state.y, pair_map, secstruct)
    x, y = rescale_coords(state.x, state.y, state.params.PRIMARY_SPACE)
    if check_overlaps(x, y, pair_map, OverlapParams()).passed:
        return x, y
    return state.x, state.y


def _verify_clean(x: list[float], y: list[float], pair_map: list[int], secstruct: str) -> None:
    """Belt-and-suspenders: never hand back a checker-dirty layout silently.

    Args:
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.
        secstruct: The original structure, for the error message.

    Raises:
        EngineError: If `check_overlaps` reports any witness.
    """
    report = check_overlaps(x, y, pair_map, OverlapParams())
    if not report.passed:
        raise EngineError(
            f"ConstructiveEngine produced a dirty layout for {secstruct!r} "
            f"({report.num_overlaps} overlaps) -- refusing to return it silently"
        )


def _place_exterior(state: _LayoutState) -> None:
    """Place the exterior loop's one branch (M1's only supported exterior shape).

    Args:
        state: Layout accumulators; `state.x`/`state.y` are filled in place.

    Raises:
        EngineError: If the exterior loop is not exactly one branch with no
            dangling tails (multi-branch/tailed exteriors are a later
            milestone, see S7).
    """
    exterior = state.tree.exterior
    if len(exterior.children) != 1 or len(exterior.members) != 1:
        raise EngineError(
            "ConstructiveEngine (M1) only supports an exterior loop with exactly "
            "one top-level branch and no dangling 5'/3' tails"
        )
    branch = exterior.children[0]
    loop, tip_a, tip_b = _place_collapsed_stem(state, branch.closing_pair, (0.0, 0.0), (1.0, 0.0))
    _place_loop_members(state, loop, tip_a, tip_b, (1.0, 0.0))


def _place_hairpin_branch(
    state: _LayoutState, branch: Branch, base_point: Point, axis_dir: Point
) -> None:
    """Place one multiloop child: its straight stem, then its hairpin loop.

    Args:
        state: Layout accumulators.
        branch: The child branch to place.
        base_point: Where this branch's outermost rung attaches (a point on
            the parent loop's circle).
        axis_dir: Unit vector pointing radially outward from the parent
            loop's center through `base_point`.

    Raises:
        EngineError: If `branch`, once its stacked pairs are collapsed,
            closes a loop with further branches -- M1 only supports plain
            hairpin children (nested multiloops are a later milestone).
    """
    loop, tip_a, tip_b = _place_collapsed_stem(state, branch.closing_pair, base_point, axis_dir)
    if loop.children:
        raise EngineError(
            "ConstructiveEngine (M1) only supports plain-hairpin multiloop "
            f"children; branch {branch.closing_pair} closes a loop with further "
            "branches (nested multiloops are a later milestone)"
        )
    _place_loop_members(state, loop, tip_a, tip_b, axis_dir)


def _place_collapsed_stem(
    state: _LayoutState, closing_pair: tuple[int, int], base_point: Point, axis_dir: Point
) -> tuple[Loop, Point, Point]:
    """Collapse a run of stacked pairs from `closing_pair` and place it as
    one straight ladder.

    Args:
        state: Layout accumulators; `state.x`/`state.y` are filled in place
            for every nucleotide in the collapsed stem.
        closing_pair: The outermost pair of the stem to place.
        base_point: Center of the stem's outermost rung.
        axis_dir: Unit vector the stem rises along, base toward tip.

    Returns:
        `(loop, tip_a, tip_b)`: the loop this stem closes (not yet placed),
        and the innermost rung's two strand coordinates (that loop's own
        closing-pair coordinates).
    """
    depth, loop = _collapse_stem(state.tree, closing_pair)
    ladder = place_stem(
        depth, base_point, axis_dir, state.params.PRIMARY_SPACE, state.params.PAIR_SPACE
    )
    _write_ladder(state, ladder, closing_pair, depth)
    return loop, ladder.strand_a[-1], ladder.strand_b[-1]


def _collapse_stem(tree: StructureTree, start_pair: tuple[int, int]) -> tuple[int, Loop]:
    """Walk a run of consecutive stacked pairs to the loop that branches.

    A loop is a pure "stack continuation" -- not a real branch point -- when
    it has exactly one child and no unpaired members (`Loop.members` is then
    just `[closing_pair[0], child.start, closing_pair[1]]`). A straight
    helix in a real RNA drawing spans the whole stack, not one segment per
    stacked pair, so this collapses the whole run into one `place_stem` call.

    Args:
        tree: The structure tree containing `start_pair`.
        start_pair: The outermost pair of the stem to walk.

    Returns:
        `(depth, loop)`: the total stacked-pair count, and the `Loop` at the
        bottom of the stack (the first real branch point: a hairpin, an
        interior loop/bulge, or a multiloop).
    """
    depth = 0
    pair = start_pair
    while True:
        depth += 1
        loop = tree.loop_by_closing_pair[pair]
        is_stack_continuation = len(loop.children) == 1 and len(loop.members) == 3
        if not is_stack_continuation:
            return depth, loop
        pair = loop.children[0].closing_pair


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
    """Place a loop's interior members (unpaired nts + child stem bases) on
    a circle beyond its own closing-pair rung.

    Each child branch is recursed into via `_place_hairpin_branch`, which
    enforces the M1 "plain hairpin children only" restriction -- so this
    function is exercised both for the root multiloop (children allowed)
    and for a hairpin's own terminal loop (never has children).

    Args:
        state: Layout accumulators; mutated in place.
        loop: The loop to place; `tip_a`/`tip_b` (its own closing pair) are
            already placed by the caller's stem.
        tip_a: The closing pair's near-strand coordinate.
        tip_b: The closing pair's far-strand coordinate.
        axis_dir: Unit vector pointing from the stem base toward this loop
            (the direction the loop opens away from its parent).
    """
    params = state.params
    interior = loop.members[1:-1]
    branch_by_start = {branch.start: branch for branch in loop.children}
    reserved = _stem_rung_half_width(params)
    slots = [
        reserved if member in branch_by_start else params.NODE_R + _UNPAIRED_MARGIN
        for member in interior
    ]
    radius_floor = max(reserved, params.NODE_R, *slots) + _RADIUS_SEARCH_SLACK
    packing = pack_loop_angles(slots, reserved, radius_floor, _RADIUS_STEP)

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
        _place_hairpin_branch(state, branch, anchor, _radial_dir(loop_center, anchor))


def _stem_rung_half_width(params: DrawParameters) -> float:
    """Chord-clearance half-width a stem's outermost rung needs on a loop circle.

    Args:
        params: Target geometry (`NODE_R`, `PAIR_SPACE`).

    Returns:
        Half the rung's physical `PAIR_SPACE` separation, plus a disk
        radius and a clearance margin (see the module-level
        `_STEM_RUNG_MARGIN` comment).
    """
    return params.PAIR_SPACE / 2.0 + params.NODE_R + _STEM_RUNG_MARGIN


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
