"""Provable area-minimization pass for the constructive layout engine.

The constructive engine (`engine.py`) is overlap-free BY CONSTRUCTION (the
isotropic bounding-disk certificate, `envelope.py`'s module docstring) but
sprawls: a dense structure puzzler can't place falls to it and renders as a
stretched strip. This module adds two, and only two, provably-clean area
levers on top of the sound (optionally already `compaction`-tightened)
layout:

1. Per-branch ROTATION `theta_b` about its own attachment pivot
   (`closing_pair[0]`'s placed coordinate). The disk `envelope.branch_disk`
   proves sibling-disjoint is centered EXACTLY at that pivot, so a rigid
   rotation about it is an isometry fixing the disk -- every
   sibling/cousin disjointness relation is preserved AUTOMATICALLY, no
   re-verification needed. Rotation still moves area because the disk's
   CONTENTS are anisotropic (an elongated stem in a round disk): reorienting
   it changes where its far tip lands. The ONE primitive the disk argument
   does not cover is the branch's own DEPARTING backbone capsule
   `(j_b, j_b + 1)` (`j_b = closing_pair[1]`, `j_b + 1` the next sibling in
   structure order) -- its far endpoint is fixed while its near endpoint
   rotates, so it is checked explicitly (`_capsule_clear`).
2. Exterior 2-D FOLD: the exterior loop places dangling tails and top-level
   branches along a single open line (`pack_line_positions`,
   `engine._place_exterior`) -- the dominant cause of a big rRNA's wide
   strip. This wraps that 1-D line into rows: ROWS are separated along
   `EXTERIOR_AXIS` by each row's own DIRECTIONAL up/down reach
   (`_row_y_offsets`) -- every exterior member's content only ever advances
   along that same axis, so this is a tighter, sound analogue of
   `pack_line_positions`' isotropic-disk argument, not the isotropic bound
   itself (an isotropic bound measured on a real hard-set structure
   defeated the whole fold: one tall branch's radius alone forced every
   other row arbitrarily far away). WITHIN a row, members are spaced along
   X by their own directional lateral half-width (`_lateral_extents`),
   spliced boustrophedon-style (`_row_layout_positions`) so consecutive
   rows connect with a short, near-vertical transition instead of a long
   diagonal. Because this tighter, re-derived packing is no longer a byte-
   for-byte instance of the ORIGINAL isotropic recipe (and, measured on a
   real hard-set structure, a per-connector certificate alone still misses
   genuine NON-adjacent-member collisions), a candidate row layout is
   verified with a real, WHOLE-structure `overlap.check_overlaps` call
   before being kept -- the same O(n) frozen checker the caller
   re-verifies with regardless, so this is airtight rather than "provably
   clean modulo one named primitive" the way rotation's targeted
   `_capsule_clear` check is (see `_try_apply_exterior_rows`'s own
   docstring).

**Scope, per the plan-critic's BLOCKING correction:** these are the only two
DOFs built here. A per-loop RADIUS re-pack is deliberately NOT built: it
would feed ANISOTROPIC (tangential) half-widths into `pack_loop_angles`,
which only bounds PERPENDICULAR extent, so chord-disjoint anchors would no
longer imply disjoint subtrees (see `compaction.py`'s "local GUESS not
PROOF" scope-limitation docstring) -- not certificate-clean. Radius
tightening is already available, checker-gated, via
`compaction.compact_layout`; this module does not reimplement it.

**Pipeline placement -- why the two DOFs are NOT one combined pass:**
`rotate_branches` needs the sibling-disk disjointness relation to be
LITERALLY true at the coordinates it rotates -- true of the sound build's
own output (that is exactly what `envelope.py` proves), but NOT true once
`compaction.compact_layout` has already re-packed siblings closer together
using tighter TANGENTIAL half-widths instead of full isotropic
`branch_reach` (measured directly: running rotation after compaction
produced real overlaps on a real hard-set structure). So `rotate_branches`
must run on the SOUND, pre-compaction layout. `fold_exterior`, in contrast,
does not lean on any PRIOR disjointness relation -- it re-derives a brand
new arrangement from each member's own CURRENT isotropic footprint
(measured fresh from placed coordinates, like `compaction._tangential_
half_width` does for its own re-pack), so it is sound regardless of when it
runs, and is most effective placed AFTER compaction: only once compaction
has tightened every subtree does the exterior's own 1-D line become the
dominant remaining sprawl (measured: on the SOUND layout the isotropic
bulk of a big rRNA dwarfs the exterior line entirely, hiding the fold's
own win; after compaction it dominates). So the intended pipeline is
`sound -> rotate_branches -> compaction.compact_layout -> fold_exterior`,
each stage independently checker-gated (`engine._compact_or_keep`).

**Safety, unconditional:** every accepted move here is followed by the
caller's whole-structure `check_overlaps` (`engine._compact_or_keep`), with
a monotone revert to the best already-known-clean layout on any failure --
so a bug in either targeted certificate can only ever cost a lost (reverted)
gain, never a silent overlap.
"""

from __future__ import annotations

import math
import time
from collections.abc import Sequence
from dataclasses import dataclass

from rna_draw import overlap
from rna_draw.layout.structure_tree import Branch, Loop, StructureTree, collapse_stem
from rna_draw.overlap import OverlapParams
from rna_draw.parameters import DrawParameters

from . import compaction, envelope
from .geometry_helpers import EXTERIOR_AXIS, Point, pack_line_positions, rotate

# Wall-clock safety net, mirroring `compaction.COMPACTION_TIME_BUDGET_S`:
# stop attempting further moves once this much time has elapsed, keeping
# every move already applied (each is independently certificate-verified,
# so a partial pass stays fully sound). Kept smaller than compaction's own
# 8.0s budget -- `_compact_or_keep` (engine.py) runs this pass AFTER
# compaction, so the two budgets are additive against the "few seconds on
# 4000nt" construction-time contract.
AREA_MIN_TIME_BUDGET_S = 4.0

# A coarse angular sweep (delta from the branch's current orientation) plus
# a handful of structure-aware seeds (`_candidate_axis_dirs`) -- the plan's
# "K ~= 8-16" range, kept at the low end since each candidate that survives
# the cheap bbox-area proxy costs one O(n) certificate check.
_ROTATION_SWEEP = 6

# Row-count candidates the exterior fold tries (`_row_count_candidates`
# always adds a size-aware `sqrt` guess on top of these).
_EXTERIOR_ROW_GUESSES = (2, 3, 4)

# `_lateral_extents`'/`_row_vertical_reach`'s own checker-clearance margin,
# matching `compaction._TANGENTIAL_MARGIN`'s rationale exactly (a capsule
# anchored near the measured boundary can still project its own half-width
# almost entirely sideways).
_FOLD_MARGIN = 20.0

_AREA_EPS = 1e-9

Extent = tuple[float, float, float, float]
_EMPTY_EXTENT: Extent = (math.inf, -math.inf, math.inf, -math.inf)


@dataclass
class _AreaMinState:
    """Accumulators threaded through the area-minimization traversal.

    Args:
        tree: The structure tree being laid out.
        x: Nucleotide x-coordinates, tightened in place.
        y: Nucleotide y-coordinates, tightened in place.
        pair_map: Entry `i` holds the partner index of `i`, or `-1`.
        params: Target geometry.
        overlap_params: Checker geometry for every targeted certificate check.
        margin_scale: The sound build's own margin multiplier.
        deadline: `time.monotonic()` value past which no further move is
            attempted, or `None` for no limit.
    """

    tree: StructureTree
    x: list[float]
    y: list[float]
    pair_map: Sequence[int]
    params: DrawParameters
    overlap_params: OverlapParams
    margin_scale: float
    deadline: float | None


def rotate_branches(
    tree: StructureTree,
    x: list[float],
    y: list[float],
    pair_map: Sequence[int],
    params: DrawParameters,
    overlap_params: OverlapParams,
    cache: envelope.ReachCache,
    margin_scale: float = 1.0,
    time_budget_s: float = AREA_MIN_TIME_BUDGET_S,
) -> tuple[list[float], list[float]]:
    """DOF 1: minimize bbox area by re-orienting each eligible branch.

    MUST be called on the SOUND (pre-`compaction.compact_layout`) layout --
    see the module docstring's "Pipeline placement" section for why the
    rotation certificate depends on that.

    Args:
        tree: The structure tree `x`/`y` was built from.
        x: Checker-clean nucleotide x-coordinates (not mutated; a working
            copy is tightened and returned).
        y: Checker-clean nucleotide y-coordinates (not mutated).
        pair_map: Entry `i` holds the partner index of `i`, or `-1`.
        params: Target geometry.
        overlap_params: Checker geometry for every targeted certificate check.
        cache: The sound build's `envelope.ReachCache` (its `degree2_pinned`
            gates which loops' children stay eligible).
        margin_scale: The sound build's own margin multiplier.
        time_budget_s: Wall-clock budget; `<= 0` disables the limit.

    Returns:
        `(x, y)`: a NEW pair of coordinate lists, with bbox area never
        larger than the input's (the caller still does one final
        whole-structure verification, `engine._compact_or_keep`).
    """
    deadline = time.monotonic() + time_budget_s if time_budget_s > 0 else None
    state = _AreaMinState(
        tree=tree,
        x=list(x),
        y=list(y),
        pair_map=pair_map,
        params=params,
        overlap_params=overlap_params,
        margin_scale=margin_scale,
        deadline=deadline,
    )
    _process_loop(state, cache, tree.exterior)
    return state.x, state.y


def fold_exterior(
    tree: StructureTree,
    x: list[float],
    y: list[float],
    pair_map: Sequence[int],
    params: DrawParameters,
    overlap_params: OverlapParams,
    margin_scale: float = 1.0,
    time_budget_s: float = AREA_MIN_TIME_BUDGET_S,
) -> tuple[list[float], list[float]]:
    """DOF 2: minimize bbox area by folding the exterior's 1-D line into rows.

    Sound at any pipeline position (re-derives each member's own extent
    fresh from its CURRENT placement, `_exterior_extents`); most effective
    AFTER `compaction.compact_layout` -- see the module docstring.

    Args:
        tree: The structure tree `x`/`y` was built from.
        x: Checker-clean nucleotide x-coordinates (not mutated).
        y: Checker-clean nucleotide y-coordinates (not mutated).
        pair_map: Entry `i` holds the partner index of `i`, or `-1`.
        params: Target geometry.
        overlap_params: Checker geometry for every targeted certificate check.
        margin_scale: The sound build's own margin multiplier.
        time_budget_s: Wall-clock budget; `<= 0` disables the limit.

    Returns:
        `(x, y)`: a NEW pair of coordinate lists, with bbox area never
        larger than the input's (the caller still does one final
        whole-structure verification, `engine._compact_or_keep`).
    """
    deadline = time.monotonic() + time_budget_s if time_budget_s > 0 else None
    state = _AreaMinState(
        tree=tree,
        x=list(x),
        y=list(y),
        pair_map=pair_map,
        params=params,
        overlap_params=overlap_params,
        margin_scale=margin_scale,
        deadline=deadline,
    )
    _apply_exterior_fold(state)
    return state.x, state.y


def _deadline_hit(state: _AreaMinState) -> bool:
    """Whether `state.deadline` has passed."""
    return state.deadline is not None and time.monotonic() > state.deadline


# --------------------------------------------------------------------------
# DOF 1: per-branch rotation about its own attachment pivot.
# --------------------------------------------------------------------------


def _loop_allows_rotation(cache: envelope.ReachCache, loop: Loop) -> bool:
    """Whether `loop`'s own children may be rotated (M2b(i) preserved).

    A bulge/interior loop's sole child (`nch == 1`) must keep continuing
    the parent stem's own axis, and a pinned degree-2 loop's two children
    (`cache.degree2_pinned`) must stay collinear -- both stay untouched,
    exactly like `compaction._compact_branch`'s own skip (`compaction.py:263`).

    Args:
        cache: The sound build's `envelope.ReachCache`.
        loop: The loop whose children are being considered.

    Returns:
        Whether `loop`'s own children are eligible for rotation.
    """
    if loop.closing_pair is None:
        return True  # the exterior has no straight-chain convention to break
    nch = len(loop.children)
    if nch == 1:
        return False
    if nch == 2:
        return not cache.degree2_pinned.get(loop.closing_pair, False)
    return True


def _process_loop(state: _AreaMinState, cache: envelope.ReachCache, loop: Loop) -> None:
    """Bottom-up: rotate every descendant branch before `loop`'s own children.

    Args:
        state: Area-minimization accumulators; mutated in place.
        cache: The sound build's `envelope.ReachCache`.
        loop: The loop whose children (and their descendants) to visit.
    """
    for child in loop.children:
        if _deadline_hit(state):
            return
        _, child_loop = collapse_stem(state.tree, child.closing_pair)
        _process_loop(state, cache, child_loop)
    if not _loop_allows_rotation(cache, loop):
        return
    for child in loop.children:
        _try_rotate_branch(state, child)


def _try_rotate_branch(state: _AreaMinState, branch: Branch) -> None:
    """Try re-orienting `branch`'s whole subtree about its own attachment
    pivot; keep the smallest-bbox-area candidate whose departing capsule
    stays clear, else leave `branch` exactly as placed.

    Args:
        state: Area-minimization accumulators; mutated in place.
        branch: The branch to try rotating.
    """
    if _deadline_hit(state):
        return
    lo, hi = branch.start, branch.end
    i_b, j_b = branch.closing_pair
    pivot = (state.x[i_b], state.y[i_b])
    old_dir = compaction._axis_dir_from_rung(
        (state.x[i_b], state.y[i_b]), (state.x[j_b], state.y[j_b]), state.params.PAIR_SPACE
    )
    outside = _outside_extent(state, lo, hi)
    current_area = _extent_area(_combine_extent(outside, _subtree_extent(state, lo, hi, pivot)))
    candidates = _candidate_axis_dirs(state, lo, hi, pivot, old_dir)
    ranked = sorted(
        candidates, key=lambda d: _candidate_area(state, lo, hi, pivot, old_dir, d, outside)
    )
    for new_dir in ranked:
        if _deadline_hit(state):
            return
        area = _candidate_area(state, lo, hi, pivot, old_dir, new_dir, outside)
        if area >= current_area - _AREA_EPS:
            return  # sorted ascending: no remaining candidate can improve
        if _apply_and_verify_rotation(state, lo, hi, j_b, pivot, old_dir, new_dir):
            return


def _candidate_axis_dirs(
    state: _AreaMinState, lo: int, hi: int, pivot: Point, current_dir: Point
) -> list[Point]:
    """A coarse angular sweep plus a structure-aware seed.

    Args:
        state: Area-minimization accumulators.
        lo: First index of the branch's own range.
        hi: Last index of the branch's own range.
        pivot: The branch's own attachment (rotation center).
        current_dir: The branch's current outward axis.

    Returns:
        Candidate outward axes to try (excludes the identity/no-op).
    """
    dirs = [
        rotate(current_dir, 2.0 * math.pi * k / _ROTATION_SWEEP) for k in range(1, _ROTATION_SWEEP)
    ]
    far_dir = _farthest_point_dir(state, lo, hi, pivot)
    if far_dir is not None:
        dirs.append(far_dir)
        dirs.append((-far_dir[0], -far_dir[1]))
    return dirs


def _farthest_point_dir(state: _AreaMinState, lo: int, hi: int, pivot: Point) -> Point | None:
    """Unit direction from `pivot` to the subtree's own farthest point.

    A structure-aware seed: aligning an elongated subtree's far tip
    radially (or its opposite, tangentially) is a natural low-area
    candidate the pure angular sweep might miss between grid points.

    Args:
        state: Area-minimization accumulators.
        lo: First index of the branch's own range.
        hi: Last index of the branch's own range.
        pivot: The branch's own attachment point.

    Returns:
        The unit direction, or `None` if the whole subtree sits on `pivot`.
    """
    best_d2 = 0.0
    best = None
    for k in range(lo, hi + 1):
        dx, dy = state.x[k] - pivot[0], state.y[k] - pivot[1]
        d2 = dx * dx + dy * dy
        if d2 > best_d2:
            best_d2, best = d2, (dx, dy)
    if best is None:
        return None
    length = math.sqrt(best_d2)
    return (best[0] / length, best[1] / length)


def _apply_and_verify_rotation(
    state: _AreaMinState, lo: int, hi: int, j_b: int, pivot: Point, old_dir: Point, new_dir: Point
) -> bool:
    """Apply one rotation candidate in place; keep it iff sound, else revert.

    Two boundary capsules cross out of `[lo, hi]`, and BOTH need a targeted
    check, not just the departing one: the incoming capsule `(lo - 1, lo)`
    has both its OWN endpoints fixed (`lo` is the rotation pivot itself),
    but its shaft still terminates INSIDE the branch's own disk (at the
    pivot) -- unlike a sibling's disk, which the chord proof places entirely
    OUTSIDE this branch's disk, the incoming capsule sits where the disk
    argument gives no separation, so a rotation can swing the branch's own
    content back into its path even though the capsule itself never moves.

    Args:
        state: Area-minimization accumulators; mutated in place.
        lo: First index of the branch's own range (also the rotation
            pivot's own index; the incoming capsule's near endpoint is
            `lo - 1`).
        hi: Last index of the branch's own range.
        j_b: The branch's closing pair's far index (the departing capsule's
            near endpoint is `j_b`; the far endpoint is `j_b + 1`).
        pivot: The rotation center (unchanged by this move).
        old_dir: The branch's outward axis before this candidate.
        new_dir: The branch's outward axis this candidate proposes.

    Returns:
        Whether the candidate was kept.
    """
    saved_x = state.x[lo : hi + 1]
    saved_y = state.y[lo : hi + 1]
    compaction._rigid_transform_range(state.x, state.y, lo, hi, pivot, old_dir, pivot, new_dir)
    if _boundary_capsules_clear(state, lo, j_b):
        return True
    state.x[lo : hi + 1] = saved_x
    state.y[lo : hi + 1] = saved_y
    return False


def _boundary_capsules_clear(state: _AreaMinState, lo: int, j_b: int) -> bool:
    """Whether both of a branch's boundary-crossing capsules stay clear.

    Args:
        state: Area-minimization accumulators (reads current coordinates).
        lo: The branch's own start index (incoming capsule near endpoint is
            `lo - 1`).
        j_b: The branch's closing pair's far index (departing capsule near
            endpoint).

    Returns:
        Whether the incoming `(lo - 1, lo)` and departing `(j_b, j_b + 1)`
        capsules are each either absent (at a structure boundary) or clear.
    """
    n = len(state.x)
    incoming_clear = lo == 0 or _capsule_clear(
        state.x, state.y, state.pair_map, state.overlap_params, lo - 1
    )
    departing_clear = j_b + 1 >= n or _capsule_clear(
        state.x, state.y, state.pair_map, state.overlap_params, j_b
    )
    return incoming_clear and departing_clear


# --------------------------------------------------------------------------
# Shared: bounding-box extent bookkeeping + the single-capsule certificate.
# --------------------------------------------------------------------------


def _outside_extent(state: _AreaMinState, lo: int, hi: int) -> Extent:
    """The bbox extent of every index NOT in `[lo, hi]`.

    Args:
        state: Area-minimization accumulators.
        lo: First index of the excluded range.
        hi: Last index of the excluded range.

    Returns:
        `(minx, maxx, miny, maxy)`, or `_EMPTY_EXTENT` if `[lo, hi]` is the
        whole structure.
    """
    if lo == 0 and hi == len(state.x) - 1:
        return _EMPTY_EXTENT
    xs = state.x[:lo] + state.x[hi + 1 :]
    ys = state.y[:lo] + state.y[hi + 1 :]
    return (min(xs), max(xs), min(ys), max(ys))


def _subtree_extent(state: _AreaMinState, lo: int, hi: int, pivot: Point) -> Extent:
    """The bbox extent of `[lo, hi]` at its CURRENT (untransformed) coords."""
    return _subtree_extent_at(state, lo, hi, pivot, (1.0, 0.0), (1.0, 0.0))


def _subtree_extent_at(
    state: _AreaMinState, lo: int, hi: int, pivot: Point, old_dir: Point, new_dir: Point
) -> Extent:
    """The bbox extent `[lo, hi]` would have after rotating `old_dir` -> `new_dir`.

    Args:
        state: Area-minimization accumulators (reads `state.x`/`state.y`).
        lo: First index of the range.
        hi: Last index of the range.
        pivot: Rotation center.
        old_dir: The range's current outward axis.
        new_dir: The candidate outward axis.

    Returns:
        `(minx, maxx, miny, maxy)` of the rotated range, without mutating
        `state`.
    """
    cos_t = old_dir[0] * new_dir[0] + old_dir[1] * new_dir[1]
    sin_t = old_dir[0] * new_dir[1] - old_dir[1] * new_dir[0]
    minx = miny = math.inf
    maxx = maxy = -math.inf
    for k in range(lo, hi + 1):
        rel_x, rel_y = state.x[k] - pivot[0], state.y[k] - pivot[1]
        nx = pivot[0] + rel_x * cos_t - rel_y * sin_t
        ny = pivot[1] + rel_x * sin_t + rel_y * cos_t
        minx, maxx = min(minx, nx), max(maxx, nx)
        miny, maxy = min(miny, ny), max(maxy, ny)
    return (minx, maxx, miny, maxy)


def _combine_extent(a: Extent, b: Extent) -> Extent:
    """The union bbox extent of two extents (an empty extent is the identity)."""
    return (min(a[0], b[0]), max(a[1], b[1]), min(a[2], b[2]), max(a[3], b[3]))


def _extent_area(extent: Extent) -> float:
    """Bounding-box area of an extent, floored to avoid a degenerate zero."""
    minx, maxx, miny, maxy = extent
    return max(maxx - minx, 1e-9) * max(maxy - miny, 1e-9)


def _candidate_area(
    state: _AreaMinState,
    lo: int,
    hi: int,
    pivot: Point,
    old_dir: Point,
    new_dir: Point,
    outside: Extent,
) -> float:
    """Cheap (O(subtree)) global bbox-area proxy for one rotation candidate."""
    subtree = _subtree_extent_at(state, lo, hi, pivot, old_dir, new_dir)
    return _extent_area(_combine_extent(outside, subtree))


def _bbox_area(x: Sequence[float], y: Sequence[float]) -> float:
    """True bounding-box area of the whole structure's current coordinates."""
    return max(max(x) - min(x), 1e-9) * max(max(y) - min(y), 1e-9)


def _capsule_clear(
    x: Sequence[float], y: Sequence[float], pair_map: Sequence[int], params: OverlapParams, k: int
) -> bool:
    """Whether the ordinary backbone capsule `(k, k + 1)` is clear of every
    OTHER primitive at trial coordinates `x`/`y`.

    A targeted single-capsule-vs-all superset test, built entirely from the
    frozen checker's own read-only building blocks (`overlap.build_primitives`,
    `overlap.is_excluded`, `overlap.test_pair`) -- their semantics are never
    touched. This is the one primitive a disk-preserving rotation, or an
    exterior-fold row transition, can leave uncovered (see the module
    docstring). A single fixed target against every other primitive is
    already linear in structure size, the same order as building the
    primitive list itself, so a spatial hash buys nothing here and is
    skipped.

    Args:
        x: Trial nucleotide x-coordinates for the whole structure.
        y: Trial nucleotide y-coordinates for the whole structure.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or -1.
        params: Checker geometry.
        k: The capsule's near endpoint; the far endpoint is `k + 1`.

    Returns:
        Whether this one capsule overlaps nothing else.
    """
    primitives = overlap.build_primitives(x, y, pair_map, params)
    target = primitives[len(x) + k]
    for other in primitives:
        if other is target:
            continue
        if overlap.is_excluded(target, other, pair_map):
            continue
        if overlap.test_pair(target, other, params.tol) is not None:
            return False
    return True


# --------------------------------------------------------------------------
# DOF 2: exterior 2-D fold.
# --------------------------------------------------------------------------


def _apply_exterior_fold(state: _AreaMinState) -> None:
    """Try folding the exterior's members from a 1-D line into a compact
    2-D row grid, keeping the fold only if it strictly shrinks the true
    bbox area and every exterior connecting capsule stays clear
    (`_try_apply_exterior_rows`).

    Args:
        state: Area-minimization accumulators; mutated in place.
    """
    exterior = state.tree.exterior
    members = exterior.members
    if len(members) <= 1 or _deadline_hit(state):
        return
    branch_by_start = {b.start: b for b in exterior.children}
    local_bboxes = _member_local_bboxes(state, members, branch_by_start)
    pad = state.params.NODE_R + _FOLD_MARGIN * state.margin_scale
    lateral_extents = _lateral_extents(local_bboxes, pad)
    chosen = _best_row_layout(state, lateral_extents, local_bboxes, pad)
    if chosen is None:
        return
    rows, anchors = chosen
    if len(rows) <= 1:
        return  # identical to the current 1-line layout; nothing to gain
    _try_apply_exterior_rows(state, members, branch_by_start, anchors)


def _lateral_extents(local_bboxes: list[Extent], pad: float) -> list[float]:
    """Each member's own X-only (lateral) half-width, for WITHIN-row spacing.

    Every exterior member's content is placed along `EXTERIOR_AXIS`
    (straight down), so how far it reaches SIDEWAYS is a strictly smaller,
    directional quantity than an isotropic disk radius -- using it (instead
    of `envelope.branch_reach`, or a farthest-point-in-any-direction bound)
    is what lets two same-row members sit as close as their actual sideways
    footprints allow. Sound by the same argument `_row_y_offsets` uses:
    `pack_line_positions`'s own proof only needs the X-distance between two
    slots to exceed the sum of their X-half-widths, which holds regardless
    of what SHAPE (disk, or this axis-aligned rectangle bound) fills that
    half-width.

    Args:
        local_bboxes: Each member's own anchor-relative bbox
            (`_member_local_bboxes`).
        pad: A nucleotide disk radius plus checker-clearance margin, added
            uniformly (mirroring `compaction._TANGENTIAL_MARGIN`'s own
            rationale: a capsule's half-width can project almost entirely
            sideways at the boundary).

    Returns:
        One half-width per member.
    """
    return [max(-dxmin, dxmax, 0.0) + pad for dxmin, dxmax, _dymin, _dymax in local_bboxes]


def _row_count_candidates(m: int) -> list[int]:
    """A small set of candidate row counts to try (clipped to `[1, m]`)."""
    guesses = {1, max(1, round(math.sqrt(m))), *_EXTERIOR_ROW_GUESSES}
    return sorted(r for r in guesses if 1 <= r <= m)


def _partition_rows(lateral_extents: list[float], row_count: int) -> list[range]:
    """Split members into `row_count` order-preserving, WIDTH-balanced chunks.

    A real exterior mixes a few huge top-level branches with many tiny
    bare tail nts (measured on a real hard-set structure); splitting by
    equal MEMBER COUNT then leaves one row dominated by a single giant
    branch and others nearly empty, which does not actually shrink the
    bbox. Splitting by cumulative lateral width instead balances what
    matters for area.

    Args:
        lateral_extents: Each member's own X-only half-width, in order.
        row_count: Target number of rows (`>= 1`).

    Returns:
        `row_count` contiguous ranges of member positions covering
        `range(len(lateral_extents))` (the last row may be short if
        `row_count` does not evenly divide the total width).
    """
    m = len(lateral_extents)
    target = sum(2.0 * e for e in lateral_extents) / row_count
    rows: list[range] = []
    start = 0
    acc = 0.0
    for i, extent in enumerate(lateral_extents):
        acc += 2.0 * extent
        at_row_budget = acc >= target and i > start
        rows_remaining = len(rows) < row_count - 1
        if at_row_budget and rows_remaining and i + 1 < m:
            rows.append(range(start, i + 1))
            start, acc = i + 1, 0.0
    rows.append(range(start, m))
    return rows


def _row_vertical_reach(
    rows: list[range], local_bboxes: list[Extent], pad: float
) -> list[tuple[float, float]]:
    """Each row's own `(down, up)` reach along `EXTERIOR_AXIS`.

    Every exterior member's content is placed by continuing to advance
    along `EXTERIOR_AXIS` from its own anchor (`_place_branch`'s recursion
    never steps backward), so `down` (how far a row's worst member reaches
    AWAY from its own line) and `up` (how far it reaches back TOWARD the
    previous row, normally ~0) are the DIRECTIONAL quantities that actually
    bound row-to-row clearance -- tighter than treating a row like an
    isotropic disk (which the measured CRW_55320 case showed washes out
    the fold's whole benefit: one tall branch's isotropic radius alone
    would force every other row arbitrarily far away, even though that
    branch never extends sideways or upward at all).

    Args:
        rows: Member-position ranges, one per row.
        local_bboxes: Each member's own anchor-relative bbox
            (`_member_local_bboxes`).
        pad: A nucleotide disk radius plus checker-clearance margin, added
            uniformly (see `_lateral_extents`).

    Returns:
        One `(down, up)` pair per row, both `>= pad`.
    """
    reach = []
    for row in rows:
        down = max(max(-local_bboxes[p][2] for p in row), 0.0) + pad
        up = max(max(local_bboxes[p][3] for p in row), 0.0) + pad
        reach.append((down, up))
    return reach


def _row_y_offsets(row_reach: list[tuple[float, float]], primary_space: float) -> list[float]:
    """Cumulative along-axis offset for each row, gapped by adjacent reach.

    SOUND by the same "sum of consecutive gaps" argument `pack_line_positions`
    uses (its own docstring): row `k`'s offset exceeds row `k - 1`'s by at
    least `down[k - 1] + up[k]`, so for ANY `i < j` (not just adjacent), the
    total gap is at least `down[i] + up[j]` (the first and last terms of a
    sum of non-negative gaps) -- enough that row `i`'s downward reach and
    row `j`'s upward reach (each already inflated by `pad`,
    `_row_vertical_reach`) never meet.

    Args:
        row_reach: Each row's own `(down, up)` (`_row_vertical_reach`).
        primary_space: Minimum consecutive-row spacing floor.

    Returns:
        One along-axis offset per row, starting at `0.0`, strictly increasing.
    """
    offsets = [0.0]
    for i in range(1, len(row_reach)):
        gap = max(primary_space, row_reach[i - 1][0] + row_reach[i][1])
        offsets.append(offsets[-1] + gap)
    return offsets


def _row_layout_positions(
    rows: list[range],
    lateral_extents: list[float],
    local_bboxes: list[Extent],
    pad: float,
    primary_space: float,
) -> list[Point]:
    """Each member's `(x, row_y)` anchor for one row partition.

    ROWS are spaced along `EXTERIOR_AXIS` by their own directional up/down
    reach (`_row_y_offsets`); WITHIN a row, members are spaced along X by
    their own directional lateral half-width (`_lateral_extents`) -- both
    are tighter, directional analogues of `pack_line_positions`' isotropic
    disk model, sound for the same reason (see each helper's own docstring).

    Rows alternate direction (boustrophedon/"snake"), each one CONTINUING
    from the X position the previous row ended at, rather than every row
    independently restarting at `x = 0`: this is what keeps the new
    row-TRANSITION capsule (the one primitive `_capsule_clear` must verify,
    see the module docstring) short and close to vertical instead of a long
    diagonal that cuts back across the whole structure -- measured on a
    real hard-set structure, restarting every row at 0 left EVERY
    transition capsule clipping unrelated content; alternating direction
    fixed it. Reversing a row's own X-order is still sound: it is a pure
    reflection + translation of that row's own `pack_line_positions` output,
    which preserves every pairwise X-distance within the row exactly, so
    the within-row disjointness proof is untouched.

    Args:
        rows: Member-position ranges, one per row, in structure order.
        lateral_extents: Each member's own X-only half-width, aligned to position.
        local_bboxes: Each member's own anchor-relative bbox.
        pad: A nucleotide disk radius plus checker-clearance margin.
        primary_space: Minimum consecutive-slot spacing floor.

    Returns:
        One `(x, y)` anchor per member position, in the same flat order as
        `lateral_extents`.
    """
    row_reach = _row_vertical_reach(rows, local_bboxes, pad)
    row_offsets = _row_y_offsets(row_reach, primary_space)
    anchors: list[Point] = []
    carry_x = 0.0
    for row_index, (row, offset) in enumerate(zip(rows, row_offsets)):
        row_y = offset * EXTERIOR_AXIS[1]
        local_xs = pack_line_positions([lateral_extents[p] for p in row], primary_space)
        direction = -1.0 if row_index % 2 else 1.0
        xs = [carry_x + direction * local_x for local_x in local_xs]
        anchors.extend((x_pos, row_y) for x_pos in xs)
        carry_x = xs[-1]
    return anchors


def _member_local_bboxes(
    state: _AreaMinState, members: list[int], branch_by_start: dict[int, Branch]
) -> list[Extent]:
    """Each member's own bbox extent RELATIVE to its current anchor.

    Placing a member is a pure TRANSLATION (`_place_exterior_member`'s own
    `_rigid_transform_range` call keeps `axis_dir` fixed at `EXTERIOR_AXIS`
    in and out), so this local shape is EXACT and translation-invariant --
    unlike the isotropic `extents` (a single radius, used only for the
    disjointness packing math), this gives an accurate bbox-area proxy for
    ranking candidate row layouts, since a real branch's true footprint is
    usually far more anisotropic (tall and thin) than an isotropic disk.

    Args:
        state: Area-minimization accumulators (reads current coordinates).
        members: `exterior.members`, in structure order.
        branch_by_start: Maps a top-level branch's start index to itself.

    Returns:
        One `(dx_min, dx_max, dy_min, dy_max)` per member, relative to that
        member's own current anchor point.
    """
    local: list[Extent] = []
    for member in members:
        branch = branch_by_start.get(member)
        if branch is None:
            local.append((0.0, 0.0, 0.0, 0.0))  # a single point; `pad` covers its own disk
            continue
        anchor = (state.x[branch.start], state.y[branch.start])
        local.append(_anchor_relative_extent(state, branch.start, branch.end, anchor))
    return local


def _anchor_relative_extent(state: _AreaMinState, lo: int, hi: int, anchor: Point) -> Extent:
    """The bbox extent of `[lo, hi]`'s CURRENT coordinates, relative to `anchor`."""
    minx = min(state.x[k] - anchor[0] for k in range(lo, hi + 1))
    maxx = max(state.x[k] - anchor[0] for k in range(lo, hi + 1))
    miny = min(state.y[k] - anchor[1] for k in range(lo, hi + 1))
    maxy = max(state.y[k] - anchor[1] for k in range(lo, hi + 1))
    return (minx, maxx, miny, maxy)


def _relative_extent(extent: Extent, anchor: Point) -> Extent:
    """Shift an `anchor`-relative extent back to absolute coordinates."""
    dxmin, dxmax, dymin, dymax = extent
    return (anchor[0] + dxmin, anchor[0] + dxmax, anchor[1] + dymin, anchor[1] + dymax)


def _row_layout_true_area(anchors: list[Point], local_bboxes: list[Extent]) -> float:
    """The TRUE bbox area a candidate row layout would produce, exactly
    (translation preserves every member's own shape -- see
    `_member_local_bboxes`)."""
    combined = _EMPTY_EXTENT
    for anchor, local in zip(anchors, local_bboxes):
        combined = _combine_extent(combined, _relative_extent(local, anchor))
    return _extent_area(combined)


def _best_row_layout(
    state: _AreaMinState, lateral_extents: list[float], local_bboxes: list[Extent], pad: float
) -> tuple[list[range], list[Point]] | None:
    """The row partition with the smallest TRUE bbox area, if any beats 1 row.

    Args:
        state: Area-minimization accumulators.
        lateral_extents: Each member's own X-only half-width (for
            within-row disjointness packing).
        local_bboxes: Each member's own anchor-relative bbox (for ranking
            and row-to-row disjointness packing).
        pad: A nucleotide disk radius plus checker-clearance margin.

    Returns:
        `(rows, anchors)` for the winning partition, or `None` if no
        candidate could be evaluated (deadline hit mid-search).
    """
    m = len(lateral_extents)
    best: tuple[list[range], list[Point]] | None = None
    best_area = math.inf
    for row_count in _row_count_candidates(m):
        if _deadline_hit(state):
            break
        rows = _partition_rows(lateral_extents, row_count)
        anchors = _row_layout_positions(
            rows, lateral_extents, local_bboxes, pad, state.params.PRIMARY_SPACE
        )
        area = _row_layout_true_area(anchors, local_bboxes)
        if area < best_area:
            best_area, best = area, (rows, anchors)
    return best


def _try_apply_exterior_rows(
    state: _AreaMinState,
    members: list[int],
    branch_by_start: dict[int, Branch],
    anchors: list[Point],
) -> None:
    """Apply the candidate row layout in place; keep it iff sound, else revert.

    Gated by a real, WHOLE-structure `overlap.check_overlaps` call, not a
    targeted per-capsule certificate: unlike rotation's disk-preserving
    argument (exactly one uncovered primitive per move, cheap to name and
    check individually), the WITHIN-row X-packing's tighter (lateral-only,
    not isotropic) extents only bound each member's own CONTENT relative to
    its immediate structure-order neighbor -- measured on a real hard-set
    structure, this leaves a genuine gap for NON-adjacent members (a small
    branch's departure point can land on the "wrong" side of its own
    anchor, or two members several rows apart can drift close), which a
    connector-only check does not catch. `check_overlaps` already runs in
    O(n) via the frozen checker's own spatial hash (the same order as
    scanning every connector capsule individually, so this closes the gap
    at no extra asymptotic cost) and is the SAME frozen predicate the
    caller (`engine._compact_or_keep`) re-verifies with regardless -- using
    it here directly, rather than a bespoke partial certificate, is both
    simpler and airtight for this DOF.

    Args:
        state: Area-minimization accumulators; mutated in place.
        members: `exterior.members`, in structure order.
        branch_by_start: Maps a top-level branch's start index to itself.
        anchors: One `(x, y)` anchor per member position.
    """
    before_area = _bbox_area(state.x, state.y)
    saved_x, saved_y = list(state.x), list(state.y)
    for member, anchor in zip(members, anchors):
        _place_exterior_member(state, member, branch_by_start.get(member), anchor)
    improved = _bbox_area(state.x, state.y) < before_area - _AREA_EPS
    report = overlap.check_overlaps(state.x, state.y, state.pair_map, state.overlap_params)
    if improved and report.passed:
        return
    state.x[:] = saved_x
    state.y[:] = saved_y


def _place_exterior_member(
    state: _AreaMinState, member: int, branch: Branch | None, anchor: Point
) -> None:
    """Move one exterior member (a bare nt, or a whole branch subtree) to `anchor`."""
    if branch is None:
        state.x[member], state.y[member] = anchor
        return
    old_anchor = (state.x[branch.start], state.y[branch.start])
    compaction._rigid_transform_range(
        state.x, state.y, branch.start, branch.end, old_anchor, EXTERIOR_AXIS, anchor, EXTERIOR_AXIS
    )


__all__ = ["rotate_branches", "fold_exterior", "AREA_MIN_TIME_BUDGET_S"]
