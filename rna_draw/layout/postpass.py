"""Rigid local post-pass for residual overlap removal.

Given a computed layout that the frozen checker (`rna_draw.overlap`) still
finds a handful of overlaps in, this module tries to clear them by rigid
moves: rotating/translating whole branches (a stem plus everything nested
inside it), pivoting on the loop each overlapping pair shares, or inflating
a crowded loop outward (see `inflate_loop`). It never bends a helix or a
loop out of shape -- `rotate_range`/`translate_range` are exactly rigid, and
`inflate_loop` only ever scales a loop's own unpaired members radially and
translates each child branch as a whole, never distorting a stem -- and it
never makes a layout worse: a candidate move is kept only if it STRICTLY
reduces the checker's `num_overlaps` (see `_apply_first_improving`);
otherwise it is discarded and the next candidate is tried.

Instrumentation on the hard gate's near-clean structures (see
`benchmarks/hard_gate.py`) found the residual overlaps are almost entirely
local small-loop crowding -- unpaired-loop-nt vs unpaired-loop-nt, or a loop
nt vs its own closing-pair capsule, at sequence distance 1-2 -- not
cross-branch collisions. Rotating or translating a whole branch cannot fix
crowding *within* one loop's own boundary; `inflate_loop` (spread the loop's
unpaired members outward, dragging each child branch along rigidly) is what
targets that case.

The checker itself is radius-dependent (`OverlapParams.node_r`), and the
hard gate scores a structure by the MINIMUM overlap count over an adaptive
radius range -- so a structure that is near-clean at the floor radius but
not at the default radius still counts as near-clean to the gate. The
post-pass therefore runs its own gate/skip decision and its acceptance
checks at that floor radius, `POSTPASS_PARAMS` (`node_r=8.0`), not the
renderer's default `node_r=10.0`; see that constant's docstring.

Watch-out worth restating in code: a rotation, translation, or inflation
that clears the witness that suggested it can introduce a NEW overlap
elsewhere (e.g. with an uninvolved third sibling branch). The strict
acceptance check catches and reverts that case, so the search loop must
keep trying other candidates and witnesses rather than stopping at the
first rejected one -- see `_candidate_moves` and `_try_witnesses`.

The exterior loop (5'/3' tails and top-level branches) has no rotational
center in the same sense an interior loop does, so exterior-loop witnesses
are handled by translation only (`_exterior_translate_candidates`), pushing
a branch away from whatever it clashed with.
"""

from __future__ import annotations

import math
import statistics
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field

from rna_draw.geometry import PrimitiveId
from rna_draw.layout.base import empty_report
from rna_draw.overlap import OverlapParams, OverlapReport, Witness, check_overlaps

from .structure_tree import (
    Branch,
    Loop,
    StructureTree,
    build_structure_tree,
    lca_loop_and_branches,
    loop_center,
)


def _default_postpass_params() -> OverlapParams:
    """Fresh copy of `POSTPASS_PARAMS`, the post-pass's floor-radius geometry.

    Returns:
        A new `OverlapParams(node_r=8.0, backbone_half_width=6.0,
        pair_half_width=6.0)`, matching `POSTPASS_PARAMS`.
    """
    return OverlapParams(node_r=8.0, backbone_half_width=6.0, pair_half_width=6.0)


# The floor of the hard gate's adaptive radius range (`benchmarks/hard_gate.py`
# scores a structure by the MINIMUM overlap count over node_r in [8, 10]), and
# thus the radius at which a near-clean structure is at its cleanest. Both the
# gate/skip decision in `PostPassEngine` (`benchmarks/engines.py`) and this
# module's own acceptance checks (`PostPassConfig.params`, below) MUST use
# this, not the renderer's default `node_r=10.0` -- checking at r=10 sees
# 14-33 overlaps on a structure the gate already scores as near-clean at r=8,
# so every near-clean structure would be (wrongly) skipped as too dirty. See
# module docstring.
POSTPASS_PARAMS = _default_postpass_params()


@dataclass
class PostPassConfig:
    """Knobs for `remove_overlaps` (grouped so the entry point stays small).

    Args:
        params: Overlap-checker geometry the pass runs its acceptance
            checks against; defaults to `POSTPASS_PARAMS`, the floor radius
            a near-clean structure is cleanest at (see module docstring).
            Must match `POSTPASS_PARAMS` for the gate's clean/non-clean
            call on the result to agree with this pass's own view of it.
        max_moves: Maximum number of accepted rigid moves before giving up.
        angle_steps_deg: Candidate rotation angles (degrees) tried, in
            order, when a branch's LCA loop has a rotational center.
        translate_fractions: Fractions of the LCA loop's mean radius used
            as candidate outward-translation distances for interior loops.
        exterior_translate_steps: Fractions of the median backbone step
            used as candidate translation distances for exterior-loop
            witnesses (see module docstring).
        inflate_factors: Candidate radial scale factors tried for
            `inflate_loop` (see that function and the module docstring).
    """

    params: OverlapParams = field(default_factory=_default_postpass_params)
    max_moves: int = 12
    angle_steps_deg: tuple[float, ...] = (10, -10, 20, -20, 5, -5, 30, -30)
    translate_fractions: tuple[float, ...] = (0.5, 1.0)
    exterior_translate_steps: tuple[float, ...] = (0.5, 1.0)
    inflate_factors: tuple[float, ...] = (1.1, 1.2, 1.35, 1.5, 1.75, 2.0)


@dataclass
class PostPassResult:
    """Outcome of `remove_overlaps`.

    Args:
        x: Nucleotide x-coordinates after the pass (unchanged input if the
            layout was already clean or no move could improve it).
        y: Nucleotide y-coordinates after the pass.
        report_before: Overlap report of the input layout.
        report_after: Overlap report of the returned layout. Always
            satisfies `report_after.num_overlaps <= report_before.num_overlaps`.
        moves_applied: Number of rigid moves accepted.
    """

    x: list[float]
    y: list[float]
    report_before: OverlapReport
    report_after: OverlapReport
    moves_applied: int


def rotate_range(
    x: Sequence[float],
    y: Sequence[float],
    start: int,
    end: int,
    cx: float,
    cy: float,
    angle_rad: float,
) -> tuple[list[float], list[float]]:
    """Rigidly rotate `x[start:end+1]`/`y[start:end+1]` about `(cx, cy)`.

    Args:
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.
        start: First index of the slice to rotate (inclusive).
        end: Last index of the slice to rotate (inclusive).
        cx: Pivot x-coordinate.
        cy: Pivot y-coordinate.
        angle_rad: Rotation angle in radians (counter-clockwise positive).

    Returns:
        `(nx, ny)`: new coordinate lists, equal to the inputs outside
        `[start, end]` and exactly rotated (all pairwise distances within
        the slice preserved) inside it.
    """
    nx, ny = list(x), list(y)
    cos_a, sin_a = math.cos(angle_rad), math.sin(angle_rad)
    for i in range(start, end + 1):
        dx, dy = x[i] - cx, y[i] - cy
        nx[i] = cx + dx * cos_a - dy * sin_a
        ny[i] = cy + dx * sin_a + dy * cos_a
    return nx, ny


def translate_range(
    x: Sequence[float], y: Sequence[float], start: int, end: int, dx: float, dy: float
) -> tuple[list[float], list[float]]:
    """Rigidly translate `x[start:end+1]`/`y[start:end+1]` by `(dx, dy)`.

    Args:
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.
        start: First index of the slice to translate (inclusive).
        end: Last index of the slice to translate (inclusive).
        dx: Translation along x.
        dy: Translation along y.

    Returns:
        `(nx, ny)`: new coordinate lists, equal to the inputs outside
        `[start, end]` and shifted by `(dx, dy)` inside it.
    """
    nx, ny = list(x), list(y)
    for i in range(start, end + 1):
        nx[i] = x[i] + dx
        ny[i] = y[i] + dy
    return nx, ny


def _spread_unpaired_members(
    loop: Loop,
    x: Sequence[float],
    y: Sequence[float],
    pair_map: Sequence[int],
    center: tuple[float, float],
    factor: float,
) -> tuple[list[float], list[float]]:
    """Radially scale `loop`'s unpaired members about `center` by `factor`.

    Args:
        loop: The loop being inflated.
        x: Current nucleotide x-coordinates.
        y: Current nucleotide y-coordinates.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.
        center: `(cx, cy)`, the loop's center (its own closing pair, if
            any, is deliberately excluded -- it is paired, so `pair_map`
            already skips it: that stem belongs to the parent loop, not
            this one).
        factor: Radial scale factor (`>1` spreads outward).

    Returns:
        `(nx, ny)`: new coordinate lists, equal to the inputs except at
        each unpaired member of `loop`, moved to `center + factor *
        (pos - center)`.
    """
    cx, cy = center
    nx, ny = list(x), list(y)
    for m in loop.members:
        if pair_map[m] != -1:
            continue
        nx[m] = cx + factor * (x[m] - cx)
        ny[m] = cy + factor * (y[m] - cy)
    return nx, ny


def _translate_branch_for_inflation(
    branch: Branch,
    x: Sequence[float],
    y: Sequence[float],
    center: tuple[float, float],
    factor: float,
) -> tuple[list[float], list[float]]:
    """Rigidly translate `branch` by its attachment midpoint's inflation.

    The branch attaches to its parent loop at its closing pair `(i, j)`,
    which sits on the loop boundary. Translating the whole branch by the
    displacement its attachment midpoint would get under the same radial
    scale keeps the stem-base chord `(i, j)` rigid -- both ends move by the
    identical vector -- rather than scaling `i` and `j` independently,
    which would distort the stem width (see module docstring).

    Args:
        branch: The child branch to translate.
        x: Current nucleotide x-coordinates.
        y: Current nucleotide y-coordinates.
        center: `(cx, cy)`, the parent loop's center.
        factor: Radial scale factor applied to the attachment midpoint
            (`>1` pushes the branch further from `center`).

    Returns:
        `(nx, ny)`: new coordinate lists, equal to the inputs outside
        `[branch.start, branch.end]` and rigidly translated inside it.
    """
    cx, cy = center
    i, j = branch.closing_pair
    mid_x, mid_y = (x[i] + x[j]) / 2.0, (y[i] + y[j]) / 2.0
    scaled_x, scaled_y = cx + factor * (mid_x - cx), cy + factor * (mid_y - cy)
    return translate_range(x, y, branch.start, branch.end, scaled_x - mid_x, scaled_y - mid_y)


def inflate_loop(
    loop: Loop,
    x: Sequence[float],
    y: Sequence[float],
    pair_map: Sequence[int],
    factor: float,
) -> tuple[list[float], list[float]]:
    """Radially inflate one loop, keeping every stem straight and rigid.

    Spreads the loop's own unpaired members outward from its center
    (`_spread_unpaired_members`) and rigidly translates each child branch
    by its attachment midpoint's same radial scaling
    (`_translate_branch_for_inflation`), so a crowded loop's nucleotides
    move apart without bending any helix or stretching any stem. This is
    the move local small-loop crowding needs -- see module docstring;
    rotating or translating a whole branch cannot fix overlaps *within*
    one loop's own boundary.

    Args:
        loop: The loop to inflate.
        x: Current nucleotide x-coordinates.
        y: Current nucleotide y-coordinates.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.
        factor: Radial scale factor (`>1` spreads the loop outward).

    Returns:
        `(nx, ny)`: new coordinate lists. `loop`'s unpaired members and
        every child branch move; `loop`'s own closing pair (if any) and
        every nucleotide outside `loop`'s members/children are untouched.
    """
    center = loop_center(loop.members, x, y)
    nx, ny = _spread_unpaired_members(loop, x, y, pair_map, center, factor)
    for branch in loop.children:
        nx, ny = _translate_branch_for_inflation(branch, nx, ny, center, factor)
    return nx, ny


def _end_indices(pid: PrimitiveId, pair_map: Sequence[int]) -> tuple[int, int]:
    """Both nucleotide indices a primitive id identifies, lower index first.

    Args:
        pid: A witness-side primitive id (`kind` in `{"nt", "bb", "pair"}`).
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.

    Returns:
        `(lower, upper)`; equal for `kind == "nt"` (a single nucleotide).

    Raises:
        ValueError: If `pid.kind` is not one of the three known kinds.
    """
    if pid.kind == "nt":
        return pid.index, pid.index
    if pid.kind == "bb":
        return pid.index, pid.index + 1
    if pid.kind == "pair":
        partner = pair_map[pid.index]
        return (pid.index, partner) if pid.index < partner else (partner, pid.index)
    raise ValueError(f"unknown witness primitive kind: {pid.kind!r}")


def witness_nucleotides(witness: Witness, pair_map: Sequence[int]) -> tuple[int, int]:
    """Representative (lower-index) nucleotide for each side of a witness.

    Args:
        witness: A checker witness.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.

    Returns:
        `(ka, kb)`, the lower-index nucleotide of `id_a` and of `id_b`.
    """
    ka, _ = _end_indices(witness.id_a, pair_map)
    kb, _ = _end_indices(witness.id_b, pair_map)
    return ka, kb


def _alt_witness_nucleotides(witness: Witness, pair_map: Sequence[int]) -> tuple[int, int]:
    """The other-endpoint nucleotide for each side of a witness (R6 retry).

    Args:
        witness: A checker witness.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.

    Returns:
        `(ka, kb)`, the upper-index nucleotide of `id_a` and of `id_b`
        (identical to `witness_nucleotides` for a `"nt"`-kind side, which
        has only one endpoint).
    """
    _, ka = _end_indices(witness.id_a, pair_map)
    _, kb = _end_indices(witness.id_b, pair_map)
    return ka, kb


def _resolve_witness(
    tree: StructureTree, witness: Witness, pair_map: Sequence[int]
) -> tuple[Loop, Branch | None, int, Branch | None, int]:
    """Resolve a witness to its LCA loop and each side's branch (if any).

    Tries the lower-index representative nucleotide of each side first. If
    that resolves to no branch on EITHER side (both nucleotides are bare
    loop members), retries with each side's other capsule endpoint, since a
    `"bb"`/`"pair"` capsule straddling a loop boundary can have its lower
    endpoint on the loop while its upper endpoint is inside a branch.

    Args:
        tree: The structure tree built from `pair_map`.
        witness: The overlap witness to resolve.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.

    Returns:
        `(loop, branch_a, ka, branch_b, kb)`: the LCA loop, side A's
        branch (or `None`) and representative nucleotide, and side B's.
    """
    ka, kb = witness_nucleotides(witness, pair_map)
    loop, branch_a, branch_b = lca_loop_and_branches(tree, ka, kb)
    if branch_a is not None or branch_b is not None:
        return loop, branch_a, ka, branch_b, kb

    alt_ka, alt_kb = _alt_witness_nucleotides(witness, pair_map)
    if (alt_ka, alt_kb) == (ka, kb):
        return loop, branch_a, ka, branch_b, kb
    alt_loop, alt_a, alt_b = lca_loop_and_branches(tree, alt_ka, alt_kb)
    if alt_a is not None or alt_b is not None:
        return alt_loop, alt_a, alt_ka, alt_b, alt_kb
    return loop, branch_a, ka, branch_b, kb


def _median_backbone_step(x: Sequence[float], y: Sequence[float]) -> float:
    """Median consecutive-nucleotide distance, the exterior move's unit.

    Args:
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.

    Returns:
        The median backbone step, or `0.0` for fewer than 2 nucleotides.
    """
    if len(x) < 2:
        return 0.0
    steps = [math.hypot(x[i + 1] - x[i], y[i + 1] - y[i]) for i in range(len(x) - 1)]
    return statistics.median(steps)


def _loop_radius(
    members: Sequence[int], x: Sequence[float], y: Sequence[float], center: tuple[float, float]
) -> float:
    """Mean distance of a loop's members to its center.

    Args:
        members: Nucleotide indices on the loop's boundary.
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.
        center: `(cx, cy)`, the loop's center.

    Returns:
        The mean member-to-center distance (`0.0` for an empty loop).
    """
    if not members:
        return 0.0
    cx, cy = center
    return sum(math.hypot(x[m] - cx, y[m] - cy) for m in members) / len(members)


def _rotate_candidates(
    loop: Loop, branch: Branch, x: Sequence[float], y: Sequence[float], config: PostPassConfig
) -> Iterator[tuple[list[float], list[float]]]:
    """Yield candidate rotations of `branch` about `loop`'s center.

    Args:
        loop: The branch's parent loop (supplies the rotation pivot).
        branch: The branch to rotate.
        x: Current nucleotide x-coordinates.
        y: Current nucleotide y-coordinates.
        config: Move-sweep parameters.

    Yields:
        `(x, y)` candidate coordinate lists, one per `config.angle_steps_deg`.
    """
    cx, cy = loop_center(loop.members, x, y)
    for degrees in config.angle_steps_deg:
        yield rotate_range(x, y, branch.start, branch.end, cx, cy, math.radians(degrees))


def _translate_candidates(
    loop: Loop, branch: Branch, x: Sequence[float], y: Sequence[float], config: PostPassConfig
) -> Iterator[tuple[list[float], list[float]]]:
    """Yield candidate radial translations of `branch` away from `loop`'s center.

    Args:
        loop: The branch's parent loop (supplies the center and radius).
        branch: The branch to translate.
        x: Current nucleotide x-coordinates.
        y: Current nucleotide y-coordinates.
        config: Move-sweep parameters.

    Yields:
        `(x, y)` candidate coordinate lists, one per
        `config.translate_fractions`; nothing if the branch sits exactly
        on the loop center (no outward direction is defined).
    """
    center = loop_center(loop.members, x, y)
    cx, cy = center
    bx, by = loop_center(list(range(branch.start, branch.end + 1)), x, y)
    dx, dy = bx - cx, by - cy
    distance = math.hypot(dx, dy)
    if distance == 0.0:
        return
    ux, uy = dx / distance, dy / distance
    radius = _loop_radius(loop.members, x, y, center)
    for fraction in config.translate_fractions:
        step = fraction * radius
        yield translate_range(x, y, branch.start, branch.end, ux * step, uy * step)


def _exterior_translate_candidates(
    branch: Branch, other_k: int, x: Sequence[float], y: Sequence[float], config: PostPassConfig
) -> Iterator[tuple[list[float], list[float]]]:
    """Yield candidate translations of `branch` away from `other_k`.

    The exterior loop has no rotational center, so exterior-loop witnesses
    MUST translate rather than rotate (see module docstring).

    Args:
        branch: The top-level branch to translate.
        other_k: Nucleotide index of the witness's other side, used as the
            point to push `branch` away from.
        x: Current nucleotide x-coordinates.
        y: Current nucleotide y-coordinates.
        config: Move-sweep parameters.

    Yields:
        `(x, y)` candidate coordinate lists, one per
        `config.exterior_translate_steps`; nothing if `branch`'s centroid
        coincides exactly with `other_k`.
    """
    bx, by = loop_center(list(range(branch.start, branch.end + 1)), x, y)
    dx, dy = bx - x[other_k], by - y[other_k]
    distance = math.hypot(dx, dy)
    if distance == 0.0:
        return
    ux, uy = dx / distance, dy / distance
    unit_step = _median_backbone_step(x, y)
    for fraction in config.exterior_translate_steps:
        step = fraction * unit_step
        yield translate_range(x, y, branch.start, branch.end, ux * step, uy * step)


def _branch_candidates(
    loop: Loop,
    branch: Branch | None,
    other_k: int,
    x: Sequence[float],
    y: Sequence[float],
    config: PostPassConfig,
) -> Iterator[tuple[list[float], list[float]]]:
    """Yield every candidate move for one witness side's branch, if any.

    Args:
        loop: The witness's LCA loop.
        branch: This side's branch, or `None` if this side is a bare loop
            member (no move is generated).
        other_k: The witness's other-side representative nucleotide.
        x: Current nucleotide x-coordinates.
        y: Current nucleotide y-coordinates.
        config: Move-sweep parameters.

    Yields:
        `(x, y)` candidate coordinate lists: exterior translations if
        `loop` is the exterior loop, else rotations then radial
        translations about `loop`'s center.
    """
    if branch is None:
        return
    if loop.closing_pair is None:
        yield from _exterior_translate_candidates(branch, other_k, x, y, config)
        return
    yield from _rotate_candidates(loop, branch, x, y, config)
    yield from _translate_candidates(loop, branch, x, y, config)


def _loops_containing(tree: StructureTree, ka: int, kb: int) -> Iterator[Loop]:
    """Yield every loop in `tree` with `ka` or `kb` directly on its boundary.

    Args:
        tree: The structure tree built from `pair_map`.
        ka: First nucleotide index.
        kb: Second nucleotide index.

    Yields:
        Each `Loop` in `tree.loops` (exterior included) whose `members`
        contains `ka` or `kb`, in tree-build order.
    """
    for loop in tree.loops:
        if ka in loop.members or kb in loop.members:
            yield loop


def _inflate_candidates(
    tree: StructureTree,
    witness: Witness,
    x: Sequence[float],
    y: Sequence[float],
    pair_map: Sequence[int],
    config: PostPassConfig,
) -> Iterator[tuple[list[float], list[float]]]:
    """Yield candidate inflations of every loop either witness side sits on.

    Local small-loop crowding (see module docstring) is fixed by spreading
    the crowded loop itself, not by rotating or translating a branch, so
    this tries every loop the witness's representative nucleotides are a
    direct member of, at each of `config.inflate_factors`.

    Args:
        tree: The structure tree built from `pair_map`.
        witness: The overlap witness driving this candidate search.
        x: Current nucleotide x-coordinates.
        y: Current nucleotide y-coordinates.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.
        config: Move-sweep parameters.

    Yields:
        `(x, y)` candidate coordinate lists, one per `(loop, factor)` pair.
    """
    ka, kb = witness_nucleotides(witness, pair_map)
    for loop in _loops_containing(tree, ka, kb):
        for factor in config.inflate_factors:
            yield inflate_loop(loop, x, y, pair_map, factor)


def _candidate_moves(
    tree: StructureTree,
    witness: Witness,
    x: Sequence[float],
    y: Sequence[float],
    pair_map: Sequence[int],
    config: PostPassConfig,
) -> Iterator[tuple[list[float], list[float]]]:
    """Yield every candidate rigid move that might clear `witness`.

    Args:
        tree: The structure tree built from `pair_map`.
        witness: The overlap witness driving this candidate search.
        x: Current nucleotide x-coordinates.
        y: Current nucleotide y-coordinates.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.
        config: Move-sweep parameters.

    Yields:
        `(x, y)` candidate coordinate lists: side A's branch moves, side
        B's, then loop-inflation candidates for the loops either side's
        representative nucleotide sits on.
    """
    loop, branch_a, ka, branch_b, kb = _resolve_witness(tree, witness, pair_map)
    yield from _branch_candidates(loop, branch_a, kb, x, y, config)
    yield from _branch_candidates(loop, branch_b, ka, x, y, config)
    yield from _inflate_candidates(tree, witness, x, y, pair_map, config)


def _apply_first_improving(
    x: Sequence[float],
    y: Sequence[float],
    pair_map: Sequence[int],
    config: PostPassConfig,
    candidates: Iterator[tuple[list[float], list[float]]],
    current_count: int,
) -> tuple[list[float], list[float], int] | None:
    """Evaluate candidates lazily; accept the first that strictly improves.

    Args:
        x: Current nucleotide x-coordinates. Unused directly (each
            candidate already carries its own full coordinate lists); kept
            in the signature to match the module's other move-search
            helpers and read naturally at call sites.
        y: Current nucleotide y-coordinates. Unused directly, see `x`.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.
        config: Supplies the overlap-checker params to re-check each
            candidate against.
        candidates: Candidate `(x, y)` coordinate lists to try, in order.
        current_count: The overlap count a candidate must beat.

    Returns:
        `(x, y, count)` for the first candidate whose overlap count is
        strictly less than `current_count` (and every coordinate is
        finite), or `None` if no candidate improves on it.
    """
    for nx, ny in candidates:
        if not (all(math.isfinite(v) for v in nx) and all(math.isfinite(v) for v in ny)):
            continue
        count = check_overlaps(nx, ny, pair_map, config.params).num_overlaps
        if count < current_count:
            return nx, ny, count
    return None


def _try_witnesses(
    tree: StructureTree,
    witnesses: Sequence[Witness],
    x: Sequence[float],
    y: Sequence[float],
    pair_map: Sequence[int],
    config: PostPassConfig,
    current_count: int,
) -> tuple[list[float], list[float], int] | None:
    """Try each witness's candidate moves in turn for one improving move.

    A rotation that clears one witness can create a new overlap elsewhere
    (see module docstring); trying every witness here, not just the first,
    is what keeps the search from giving up too early.

    Args:
        tree: The structure tree built from `pair_map`.
        witnesses: The current overlap report's witnesses, in report order.
        x: Current nucleotide x-coordinates.
        y: Current nucleotide y-coordinates.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.
        config: Move-sweep parameters.
        current_count: The overlap count a candidate must beat.

    Returns:
        `(x, y, count)` for the first witness whose candidates include an
        improving move, or `None` if none of them do.
    """
    for witness in witnesses:
        candidates = _candidate_moves(tree, witness, x, y, pair_map, config)
        applied = _apply_first_improving(x, y, pair_map, config, candidates, current_count)
        if applied is not None:
            return applied
    return None


def _degenerate_report(
    x: Sequence[float], y: Sequence[float], pair_map: Sequence[int], params: OverlapParams
) -> OverlapReport:
    """Overlap report for inputs too small for the checker to accept.

    Args:
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.
        params: Overlap-checker geometry parameters.

    Returns:
        A zero-witness report for a length-0 input (which `check_overlaps`
        rejects), else the real `check_overlaps` report.
    """
    if len(x) == 0:
        return empty_report()
    return check_overlaps(x, y, pair_map, params)


def remove_overlaps(
    x: Sequence[float],
    y: Sequence[float],
    pair_map: Sequence[int],
    config: PostPassConfig | None = None,
) -> PostPassResult:
    """Remove residual overlaps from a layout by rigid local moves.

    Computes the baseline report; while the layout is still dirty and
    under the move budget, picks a witness, builds candidate rigid moves
    for the branches it implicates, and applies the first one that
    strictly reduces the overlap count. Stops at zero overlaps, at
    `config.max_moves` accepted moves, or when no witness's candidates
    improve on the current count. Never returns a layout worse than the
    input (see module docstring).

    Args:
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired. Must be well-nested (pseudoknot-free).
        config: Move-sweep and checker parameters; defaults to
            `PostPassConfig()`.

    Returns:
        A `PostPassResult` with `report_after.num_overlaps <=
        report_before.num_overlaps`, guaranteed. A clean or fewer-than-2
        nucleotide input is returned unchanged (no moves).
    """
    config = config or PostPassConfig()
    if len(x) < 2:
        report = _degenerate_report(x, y, pair_map, config.params)
        return PostPassResult(list(x), list(y), report, report, 0)

    report_before = check_overlaps(x, y, pair_map, config.params)
    if report_before.passed:
        return PostPassResult(list(x), list(y), report_before, report_before, 0)

    tree = build_structure_tree(pair_map)
    cur_x, cur_y = list(x), list(y)
    cur_count = report_before.num_overlaps
    witnesses = report_before.witnesses
    moves = 0
    while cur_count > 0 and moves < config.max_moves:
        applied = _try_witnesses(tree, witnesses, cur_x, cur_y, pair_map, config, cur_count)
        if applied is None:
            break
        cur_x, cur_y, cur_count = applied
        moves += 1
        witnesses = check_overlaps(cur_x, cur_y, pair_map, config.params).witnesses

    report_after = check_overlaps(cur_x, cur_y, pair_map, config.params)
    return PostPassResult(cur_x, cur_y, report_before, report_after, moves)


__all__ = [
    "POSTPASS_PARAMS",
    "PostPassConfig",
    "PostPassResult",
    "rotate_range",
    "translate_range",
    "inflate_loop",
    "witness_nucleotides",
    "remove_overlaps",
]
