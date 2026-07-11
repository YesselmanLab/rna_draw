"""Checker-gated monotone compaction pass (M3, see the constructive-engine
runbook's milestone M3, S10-S11).

`ConstructiveEngine` (M2) is clean BY CONSTRUCTION but SPRAWLING: every
multiloop packs its children using `envelope.branch_reach`, an ISOTROPIC
disk radius bounding a child's WHOLE subtree -- a huge over-estimate of the
angular footprint an elongated subtree (e.g. a long straight stem) actually
needs (see `envelope.py`'s module docstring for why that isotropic bound
was the right choice for SOUNDNESS; it is not the right choice for
compactness).

This module tightens the sound layout AFTER the fact, replacing each
child's isotropic reach with a TANGENTIAL half-width measured from its own
ALREADY-PLACED (and already checker-clean) coordinates -- the sideways
extent perpendicular to the child's own outward axis, which is usually far
smaller than the isotropic bound for anything that isn't itself round.
Re-packing a loop's children with these tighter half-widths (the same
`geometry_helpers.pack_loop_angles` exact-interval packer the sound build
uses) yields a smaller candidate loop radius; each child's WHOLE subtree is
then moved onto its new slot with a single RIGID transform (translate +
rotate, no scaling), which preserves every internal distance exactly -- so
a child that is already checker-clean internally STAYS checker-clean after
being moved, by construction (Euclidean distances are transform-invariant).

Each candidate is verified with a LOCAL `check_overlaps` call scoped to
just this loop's own closing-pair index range (which fully contains the
loop's members and every descendant subtree -- nothing OUTSIDE that range
is written), not a whole-structure re-check. A move is applied only if the
local check stays clean and strictly smaller than the sound radius;
otherwise the loop is left exactly as the sound build left it (monotone: a
rejected move never regresses anything). Traversal is bottom-up (children
compacted before their parent), so an outer loop's tangential measurement
already reflects its children's own compaction gains. The exterior loop
(an open line, not a circle) gets the same treatment via
`geometry_helpers.pack_line_positions`, checked once against the whole
structure (there is no smaller enclosing range for an open boundary).

**Known scope limitation (why the local check is still only a local
GUESS, not a local PROOF):** re-packing a child onto a new slot ROTATES
it (its axis now points at a different angle), not just moves it closer.
The local check only re-verifies content WITHIN this loop's own index
range -- it cannot see a rotation that redirects a deeply-nested
grandchild toward a completely unrelated sibling TOP-LEVEL branch several
levels up the tree, outside this loop's own range. This is rare (measured:
a handful of the largest, many-top-level-branch hard-set structures) and
NEVER unsound -- the caller (`engine._compact_or_keep`) still re-verifies
the WHOLE structure once at the end as the real belt-and-suspenders check,
and falls back to the untouched pre-compaction (sound) layout if that
somehow fails -- but it means those specific structures get NONE of
compaction's gain rather than most of it (all-or-nothing at the whole-
structure level). Recovering partial credit there (e.g. a global blend/
bisect between the sound and compacted layouts) was prototyped but not
kept: it cost 10-15x more wall-clock time on the largest structures for a
sub-2% aggregate size reduction, a poor trade against the "few seconds on
4000nt" budget -- see the runbook's Run log / handoff report.
"""

from __future__ import annotations

import time
from collections.abc import Sequence
from dataclasses import dataclass

from rna_draw.layout.structure_tree import Branch, Loop, StructureTree, collapse_stem
from rna_draw.overlap import OverlapParams, check_overlaps
from rna_draw.parameters import DrawParameters

from . import envelope
from .geometry_helpers import (
    EXTERIOR_AXIS,
    LoopPacking,
    Point,
    loop_member_point,
    midpoint,
    pack_line_positions,
    pack_loop_angles,
    rotate,
    rotate90_ccw,
    translate,
)

# A slot's tangential half-width also needs a checker-clearance margin: the
# X-projection alone only bounds a NUCLEOTIDE's disk (`NODE_R`, added
# separately below); a CAPSULE (backbone/pair) anchored near the boundary
# but oriented steeply relative to the packing axis can still project its
# own half-width (`backbone_half_width`/`pair_half_width`, both
# `0.75 * NODE_R` by default) almost entirely sideways -- measured
# empirically (a real hard-set structure) that `NODE_R` alone left as
# little as ~1e-12 units of true clearance, rejecting an otherwise-good
# move. This margin is still much smaller than `envelope.RHO_MARGIN`
# (17.5, deliberately generous for the sound build's PROOF-based bound)
# because this pass gates every move with a real local `check_overlaps`
# call: an imperfect margin only means a smaller gain (or a rejected
# move), never a silent overlap -- so it is fine to be a little generous.
_TANGENTIAL_MARGIN = 20.0

_RADIUS_STEP = 1.0
_RADIUS_SEARCH_SLACK = 1.0

# Wall-clock safety net: stop attempting further loop compaction once this
# much time has elapsed since the pass started, keeping every move already
# applied (each was independently checker-verified, so a partial pass is
# still fully sound -- this only bounds compaction's OWN cost, never
# construction's; see the runbook's "a few seconds on 4000nt" requirement).
COMPACTION_TIME_BUDGET_S = 8.0


@dataclass
class _CompactionState:
    """Accumulators threaded through the compaction traversal.

    Args:
        tree: The structure tree being compacted.
        x: Nucleotide x-coordinates, tightened in place.
        y: Nucleotide y-coordinates, tightened in place.
        pair_map: Entry `i` holds the partner index of `i`, or `-1`.
        params: Target geometry (`NODE_R`/`PRIMARY_SPACE`/`PAIR_SPACE`).
        overlap_params: Checker geometry for every local gate check.
        margin_scale: The sound build's own margin multiplier (S8's
            bounded re-inflation schedule) -- reused so `reserved`/
            `unpaired_half_width` match what was actually proven clean.
        deadline: `time.monotonic()` value past which no further loop is
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


@dataclass
class _RepackPlan:
    """A candidate tighter packing for one multiloop, ready to apply.

    A child's CURRENT axis is deliberately not cached here: it is
    recomputed fresh from `state.x`/`state.y` inside `_apply_multiloop_packing`
    instead, since that is always correct regardless of what has already
    been applied to `state`, whereas a value computed earlier (e.g. in
    `_multiloop_half_widths`) would go stale the moment anything upstream
    of this plan rotates the child.

    Args:
        lo: First index of the loop's own closing-pair range.
        hi: Last index of the loop's own closing-pair range.
        interior: The loop's interior member indices, in structure order.
        branch_by_start: Maps a child branch's start index to itself.
        center: Midpoint of the loop's own closing-pair rung.
        axis_dir: The loop's own outward axis (unit vector).
        packing: The candidate tighter `LoopPacking`.
    """

    lo: int
    hi: int
    interior: list[int]
    branch_by_start: dict[int, Branch]
    center: Point
    axis_dir: Point
    packing: LoopPacking


def compact_layout(
    tree: StructureTree,
    x: list[float],
    y: list[float],
    pair_map: Sequence[int],
    params: DrawParameters,
    overlap_params: OverlapParams,
    cache: envelope.ReachCache,
    margin_scale: float = 1.0,
    time_budget_s: float = COMPACTION_TIME_BUDGET_S,
) -> tuple[list[float], list[float]]:
    """Tighten a sound, checker-clean layout without ever accepting a dirty move.

    Args:
        tree: The structure tree the sound layout was built from.
        x: Sound nucleotide x-coordinates (not mutated; a working copy is
            tightened and returned).
        y: Sound nucleotide y-coordinates (not mutated).
        pair_map: Entry `i` holds the partner index of `i`, or `-1`.
        params: Target geometry.
        overlap_params: Checker geometry for every local gate check.
        cache: The sound build's `envelope.ReachCache` (its `loop_packing`
            gives each loop's sound radius, compared against every
            compaction candidate so a move is only ever applied if it is
            strictly smaller).
        margin_scale: The sound build's own margin multiplier, reused for
            `reserved`/`unpaired_half_width` so they stay consistent.
        time_budget_s: Wall-clock budget; `<= 0` disables the limit.

    Returns:
        `(x, y)`: a NEW pair of coordinate lists, at least as compact as
        the input and still checker-clean on every locally-verified move
        (the caller still does one final whole-structure verification).
    """
    deadline = time.monotonic() + time_budget_s if time_budget_s > 0 else None
    state = _CompactionState(
        tree=tree,
        x=list(x),
        y=list(y),
        pair_map=pair_map,
        params=params,
        overlap_params=overlap_params,
        margin_scale=margin_scale,
        deadline=deadline,
    )
    _compact_exterior(state, cache)
    return state.x, state.y


def _deadline_hit(state: _CompactionState) -> bool:
    """Whether `state.deadline` has passed."""
    return state.deadline is not None and time.monotonic() > state.deadline


def _compact_exterior(state: _CompactionState, cache: envelope.ReachCache) -> None:
    """Compact every top-level branch, then the exterior's own line spacing."""
    exterior = state.tree.exterior
    branch_by_start = {b.start: b for b in exterior.children}
    for member in exterior.members:
        if _deadline_hit(state):
            return
        branch = branch_by_start.get(member)
        if branch is not None:
            _compact_branch(state, cache, branch)
    if not _deadline_hit(state):
        _compact_exterior_spacing(state)


def _compact_branch(state: _CompactionState, cache: envelope.ReachCache, branch: Branch) -> None:
    """Compact a branch's children first (bottom-up), then its own loop."""
    if _deadline_hit(state):
        return
    _, loop = collapse_stem(state.tree, branch.closing_pair)
    if len(loop.children) >= 2:
        for child in loop.children:
            _compact_branch(state, cache, child)
        if not _deadline_hit(state):
            _compact_multiloop(state, cache, loop)
    elif len(loop.children) == 1:
        _compact_branch(state, cache, loop.children[0])


def _axis_dir_from_rung(tip_a: Point, tip_b: Point, pair_space: float) -> Point:
    """Recover a placed stem rung's `axis_dir` from its two coordinates.

    `place_stem` offsets both strands `pair_space / 2` to either side of
    `axis_dir`'s perpendicular, so `tip_b - tip_a == pair_space * perp`
    exactly for ANY rung of the stem -- this inverts that relationship.

    Args:
        tip_a: The rung's near-strand coordinate.
        tip_b: The rung's far-strand coordinate, `pair_space` away.
        pair_space: The rung's known separation.

    Returns:
        The stem's `axis_dir` (base toward tip).
    """
    perp = ((tip_b[0] - tip_a[0]) / pair_space, (tip_b[1] - tip_a[1]) / pair_space)
    return (perp[1], -perp[0])


def _tangential_half_width(
    state: _CompactionState, anchor: Point, axis_dir: Point, lo: int, hi: int
) -> float:
    """The chord-clearance half-width `[lo, hi]`'s subtree actually needs.

    Measured as the largest perpendicular (tangential) distance from the
    line through `anchor` along `axis_dir` to any of the subtree's ALREADY
    -PLACED points -- usually far smaller than the isotropic
    `envelope.branch_reach` bound for anything that isn't round (see the
    module docstring).

    Args:
        state: Compaction accumulators (reads `state.x`/`state.y`).
        anchor: The subtree's attachment point.
        axis_dir: The subtree's own outward axis (unit vector).
        lo: First nucleotide index of the subtree's index range.
        hi: Last nucleotide index of the subtree's index range.

    Returns:
        The tangential extent plus a checker-clearance margin.
    """
    perp = rotate90_ccw(axis_dir)
    x, y = state.x, state.y
    max_perp = 0.0
    for k in range(lo, hi + 1):
        rel_x, rel_y = x[k] - anchor[0], y[k] - anchor[1]
        proj = abs(rel_x * perp[0] + rel_y * perp[1])
        if proj > max_perp:
            max_perp = proj
    return max_perp + state.params.NODE_R + _TANGENTIAL_MARGIN * state.margin_scale


def _rigid_transform_range(
    x: list[float],
    y: list[float],
    lo: int,
    hi: int,
    old_anchor: Point,
    old_dir: Point,
    new_anchor: Point,
    new_dir: Point,
) -> None:
    """Move `[lo, hi]` rigidly (translate + rotate) so `old_anchor`/`old_dir`
    maps exactly onto `new_anchor`/`new_dir`, preserving every internal
    distance (a checker-clean subtree stays checker-clean).

    Args:
        x: X-coordinates, mutated in place for indices `[lo, hi]`.
        y: Y-coordinates, mutated in place for indices `[lo, hi]`.
        lo: First index of the range to move.
        hi: Last index of the range to move.
        old_anchor: The range's current attachment point.
        old_dir: The range's current outward axis (unit vector).
        new_anchor: Where the attachment point must land.
        new_dir: The new outward axis (unit vector).
    """
    cos_t = old_dir[0] * new_dir[0] + old_dir[1] * new_dir[1]
    sin_t = old_dir[0] * new_dir[1] - old_dir[1] * new_dir[0]
    for k in range(lo, hi + 1):
        rel_x, rel_y = x[k] - old_anchor[0], y[k] - old_anchor[1]
        x[k] = new_anchor[0] + rel_x * cos_t - rel_y * sin_t
        y[k] = new_anchor[1] + rel_x * sin_t + rel_y * cos_t


def _local_clean(
    state: _CompactionState, trial_x: list[float], trial_y: list[float], lo: int, hi: int
) -> bool:
    """Whether `[lo, hi]`, checked in isolation, is checker-clean.

    Args:
        state: Compaction accumulators (reads `state.pair_map`/`overlap_params`).
        trial_x: Candidate x-coordinates for the WHOLE structure.
        trial_y: Candidate y-coordinates for the WHOLE structure.
        lo: First index of the range to check (fully contains everything
            that could have changed -- see the module docstring).
        hi: Last index of the range to check.

    Returns:
        Whether `check_overlaps` on the extracted, re-indexed local
        sub-structure passes.
    """
    local_x = trial_x[lo : hi + 1]
    local_y = trial_y[lo : hi + 1]
    local_pair_map = [
        (state.pair_map[k] - lo) if state.pair_map[k] != -1 else -1 for k in range(lo, hi + 1)
    ]
    return check_overlaps(local_x, local_y, local_pair_map, state.overlap_params).passed


def _multiloop_half_widths(
    state: _CompactionState, interior: list[int], branch_by_start: dict[int, Branch]
) -> list[float]:
    """Each interior member's tangential half-width.

    Args:
        state: Compaction accumulators.
        interior: The loop's interior member indices, in structure order.
        branch_by_start: Maps a child branch's start index to itself.

    Returns:
        One half-width per `interior` entry: a child's
        `_tangential_half_width`, or a bare member's
        `envelope.unpaired_half_width`.
    """
    half_widths: list[float] = []
    x, y = state.x, state.y
    for member in interior:
        child = branch_by_start.get(member)
        if child is None:
            half_widths.append(envelope.unpaired_half_width(state.params, state.margin_scale))
            continue
        c_i, c_j = child.closing_pair
        axis = _axis_dir_from_rung((x[c_i], y[c_i]), (x[c_j], y[c_j]), state.params.PAIR_SPACE)
        anchor = (x[child.start], y[child.start])
        half_widths.append(_tangential_half_width(state, anchor, axis, child.start, child.end))
    return half_widths


def _compact_multiloop(state: _CompactionState, cache: envelope.ReachCache, loop: Loop) -> None:
    """Try to re-pack one multiloop's children onto a tighter circle."""
    if loop.closing_pair is None:
        return
    old_packing = cache.loop_packing.get(loop.closing_pair)
    if old_packing is None:
        return
    i, j = loop.closing_pair
    x, y = state.x, state.y
    tip_a, tip_b = (x[i], y[i]), (x[j], y[j])
    axis_dir = _axis_dir_from_rung(tip_a, tip_b, state.params.PAIR_SPACE)
    center = midpoint(tip_a, tip_b)

    interior = loop.members[1:-1]
    branch_by_start = {b.start: b for b in loop.children}
    reserved = envelope.stem_rung_half_width(state.params, state.margin_scale)
    half_widths = _multiloop_half_widths(state, interior, branch_by_start)

    radius_floor = max([reserved, state.params.NODE_R, *half_widths]) + _RADIUS_SEARCH_SLACK
    new_packing = pack_loop_angles(half_widths, reserved, radius_floor, _RADIUS_STEP)
    if new_packing.radius >= old_packing.radius:
        return
    plan = _RepackPlan(i, j, interior, branch_by_start, center, axis_dir, new_packing)
    _apply_multiloop_packing(state, plan)


def _apply_multiloop_packing(state: _CompactionState, plan: _RepackPlan) -> None:
    """Build the candidate tighter multiloop layout, apply it if locally clean."""
    x, y = state.x, state.y
    trial_x, trial_y = list(x), list(y)
    new_center = translate(plan.center, plan.axis_dir, plan.packing.radius)
    zero_dir = (-plan.axis_dir[0], -plan.axis_dir[1])

    for member, angle in zip(plan.interior, plan.packing.angles):
        new_anchor = loop_member_point(new_center, zero_dir, plan.packing.radius, angle)
        child = plan.branch_by_start.get(member)
        if child is None:
            trial_x[member], trial_y[member] = new_anchor
            continue
        c_i, c_j = child.closing_pair
        old_dir = _axis_dir_from_rung((x[c_i], y[c_i]), (x[c_j], y[c_j]), state.params.PAIR_SPACE)
        old_anchor = (x[child.start], y[child.start])
        new_dir = rotate(zero_dir, angle)
        _rigid_transform_range(
            trial_x, trial_y, child.start, child.end, old_anchor, old_dir, new_anchor, new_dir
        )

    if _local_clean(state, trial_x, trial_y, plan.lo, plan.hi):
        for k in range(plan.lo, plan.hi + 1):
            x[k], y[k] = trial_x[k], trial_y[k]


def _exterior_extents(
    state: _CompactionState, members: list[int], branch_by_start: dict[int, Branch]
) -> list[float]:
    """Each exterior member's tangential half-width (see `_multiloop_half_widths`)."""
    x, y = state.x, state.y
    extents = []
    for member in members:
        branch = branch_by_start.get(member)
        if branch is None:
            extents.append(envelope.unpaired_half_width(state.params, state.margin_scale))
            continue
        anchor = (x[branch.start], y[branch.start])
        extents.append(
            _tangential_half_width(state, anchor, EXTERIOR_AXIS, branch.start, branch.end)
        )
    return extents


def _member_span(state: _CompactionState, members: list[int]) -> float:
    """Current x-extent spanned by `members`' present positions."""
    if not members:
        return 0.0
    xs = [state.x[m] for m in members]
    return max(xs) - min(xs)


def _compact_exterior_spacing(state: _CompactionState) -> None:
    """Try to re-pack the exterior's top-level branches onto a tighter line."""
    exterior = state.tree.exterior
    branch_by_start = {b.start: b for b in exterior.children}
    if not branch_by_start:
        return
    members = exterior.members
    positions = pack_line_positions(
        _exterior_extents(state, members, branch_by_start), state.params.PRIMARY_SPACE
    )
    new_span = positions[-1] - positions[0] if positions else 0.0
    if new_span >= _member_span(state, members):
        return
    _apply_exterior_positions(state, members, branch_by_start, positions)


def _apply_exterior_positions(
    state: _CompactionState,
    members: list[int],
    branch_by_start: dict[int, Branch],
    positions: list[float],
) -> None:
    """Build the tighter exterior layout, apply it if checker-clean overall."""
    x, y = state.x, state.y
    trial_x, trial_y = list(x), list(y)
    for member, new_x_pos in zip(members, positions):
        branch = branch_by_start.get(member)
        if branch is None:
            trial_x[member], trial_y[member] = new_x_pos, 0.0
            continue
        old_anchor = (x[branch.start], y[branch.start])
        _rigid_transform_range(
            trial_x,
            trial_y,
            branch.start,
            branch.end,
            old_anchor,
            EXTERIOR_AXIS,
            (new_x_pos, 0.0),
            EXTERIOR_AXIS,
        )
    if _local_clean(state, trial_x, trial_y, 0, len(x) - 1):
        x[:] = trial_x
        y[:] = trial_y


__all__ = ["compact_layout", "COMPACTION_TIME_BUDGET_S"]
