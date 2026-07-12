"""Provable area-minimization pass for the constructive layout engine.

The constructive engine (`engine.py`) is overlap-free BY CONSTRUCTION (the
isotropic bounding-disk certificate, `envelope.py`'s module docstring) but
sprawls: a dense structure puzzler can't place falls to it and renders as a
stretched strip. This module adds one provably-clean area lever on top of
the sound layout:

Per-branch ROTATION `theta_b` about its own attachment pivot
(`closing_pair[0]`'s placed coordinate). The disk `envelope.branch_disk`
proves sibling-disjoint is centered EXACTLY at that pivot, so a rigid
rotation about it is an isometry fixing the disk -- every sibling/cousin
disjointness relation is preserved AUTOMATICALLY, no re-verification
needed. Rotation still moves area because the disk's CONTENTS are
anisotropic (an elongated stem in a round disk): reorienting it changes
where its far tip lands. The ONE primitive the disk argument does not
cover is the branch's own DEPARTING backbone capsule `(j_b, j_b + 1)`
(`j_b = closing_pair[1]`, `j_b + 1` the next sibling in structure order) --
its far endpoint is fixed while its near endpoint rotates, so it is
checked explicitly (`_capsule_clear`).

**Scope:** a per-loop RADIUS re-pack is deliberately NOT built: it would
feed ANISOTROPIC (tangential) half-widths into `pack_loop_angles`, which
only bounds PERPENDICULAR extent, so chord-disjoint anchors would no
longer imply disjoint subtrees (see `compaction.py`'s "local GUESS not
PROOF" scope-limitation docstring) -- not certificate-clean. Radius
tightening is already available, checker-gated, via
`compaction.compact_layout`; this module does not reimplement it.

**Pipeline placement:** `rotate_branches` needs the sibling-disk
disjointness relation to be LITERALLY true at the coordinates it rotates --
true of the sound build's own output (that is exactly what `envelope.py`
proves), but NOT true once `compaction.compact_layout` has already
re-packed siblings closer together using tighter TANGENTIAL half-widths
instead of full isotropic `branch_reach` (measured directly: running
rotation after compaction produced real overlaps on a real hard-set
structure). So `rotate_branches` must run on the SOUND, pre-compaction
layout: `sound -> rotate_branches -> compaction.compact_layout`, each
stage independently checker-gated (`engine._compact_or_keep`).

**Safety, unconditional:** every accepted move here is followed by the
caller's whole-structure `check_overlaps` (`engine._compact_or_keep`), with
a monotone revert to the best already-known-clean layout on any failure --
so a bug in the targeted certificate can only ever cost a lost (reverted)
gain, never a silent overlap.

**History:** an earlier version also carried an exterior 2-D FOLD DOF
(wrapping the exterior's open line into rows). It was verified sound but
was rejected by the safety check (never beat the 1-line layout, or failed
the whole-structure overlap check) on every real hard-set structure and
every synthetic case engineered to favor it -- a confirmed no-op -- and
was removed for leanness.
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
from .geometry_helpers import Point, rotate

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


def _capsule_clear(
    x: Sequence[float], y: Sequence[float], pair_map: Sequence[int], params: OverlapParams, k: int
) -> bool:
    """Whether the ordinary backbone capsule `(k, k + 1)` is clear of every
    OTHER primitive at trial coordinates `x`/`y`.

    A targeted single-capsule-vs-all superset test, built entirely from the
    frozen checker's own read-only building blocks (`overlap.build_primitives`,
    `overlap.is_excluded`, `overlap.test_pair`) -- their semantics are never
    touched. This is the one primitive a disk-preserving rotation can leave
    uncovered (see the module docstring). A single fixed target against
    every other primitive is already linear in structure size, the same
    order as building the primitive list itself, so a spatial hash buys
    nothing here and is skipped.

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


__all__ = ["rotate_branches", "AREA_MIN_TIME_BUDGET_S"]
