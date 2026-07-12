"""Phase 2b: bias the nested layout to bring crossing partners close together.

`layout_pseudoknot`'s Phase 2 (`engine._layout_nested`) lays out the
pseudoknot-free nested subset WITHOUT regard to where each crossing
partner (a nucleotide `extraction.max_nested_subset` dropped from a
crossing stem) ends up -- so most crossing pairs end up far apart, forcing
Phase 3 (`placement.place_crossings`) to either bow a long PK-A polyline
across the figure or leave the pair unplaced entirely. This module is a
pure PRE-placement coordinate transform: given the nested layout's
coordinates, reflect whole nested subtrees about their attachment points
to shrink crossing-pair distances, so more crossings become short in-plane
PK-B connectors.

This module never touches `overlap.py`/`geometry.py`/`spatial_hash.py` --
it treats the frozen `check_overlaps` as the SOLE arbiter of correctness.
Every candidate move is re-verified by a whole-layout `check_overlaps`
call on a fresh coordinate copy; a move is adopted only if that check
still passes AND the objective strictly decreases. See
`bias_crossing_proximity`'s docstring for the never-worse contract.

MEASURED FINDING (real corpus, `benchmarks/pseudoknot_gate.py`, 200
structures): with only LCA-level candidates (the move set this module
ships with), the search is a safe no-op end to end -- 0/200 structures
get any accepted flip, though 184/200 DO have at least one candidate that
would lower the objective if the checker allowed it. Traced cause: a
branch's own rung axis (`_branch_frame`) is typically ALIGNED with the
direction its incoming/outgoing exterior (or multiloop) neighbor already
occupies (both a real engine's linear exterior layout and its radial loop
packing route a branch's neighbors along that same rung-perpendicular
line), so reflecting to the mirror side collides with them almost by
construction -- not merely "coarse" (the plan-critic's flagged risk), but
close to systematic for LCA-level (outermost) candidates specifically.
The reflection/search machinery itself is verified correct and effective
on hand-built layouts where a branch's neighbors are NOT on that line
(`tests/test_pseudoknot_proximity.py`'s `TestBiasCrossingProximityUnit`).
Recommended next step (deferred, not built here): also try DEEPER
candidates (a crossing endpoint's `enclosing_pairs` stack below the LCA
branch), whose smaller, more localized subtrees are less likely to sweep
across a neighbor's own layout line.
"""

from __future__ import annotations

import math

from rna_draw.layout.structure_tree import (
    Branch,
    StructureTree,
    build_structure_tree,
    lca_loop_and_branches,
)
from rna_draw.overlap import OverlapParams, check_overlaps

from .parsing import Stem

Point = tuple[float, float]


def crossing_endpoint_distance(x: list[float], y: list[float], crossing_stems: list[Stem]) -> float:
    """The proximity OBJECTIVE: summed crossing-pair endpoint distance.

    Args:
        x: Nucleotide x-coordinates (nested layout, or a trial copy).
        y: Nucleotide y-coordinates.
        crossing_stems: Every crossing stem (`extraction.max_nested_subset`'s
            second return value).

    Returns:
        `sum(hypot(x[a]-x[b], y[a]-y[b]))` over every pair of every
        crossing stem. Zero (and this module a no-op) for a pk-free input.
    """
    total = 0.0
    for stem in crossing_stems:
        for a, b in stem.pairs():
            total += math.hypot(x[a] - x[b], y[a] - y[b])
    return total


def _branch_frame(x: list[float], y: list[float], branch: Branch) -> tuple[Point, Point] | None:
    """The reflection frame for `branch`: its attachment pivot and axis.

    Args:
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.
        branch: The candidate branch (owns the contiguous index range
            `[branch.start, branch.end]`).

    Returns:
        `(pivot, axis_dir)`: `pivot` is `branch.start`'s own coordinate
        (fixed by the reflection); `axis_dir` is the unit vector from
        `branch.start` toward `branch.end` rotated -90 degrees (the
        stem's outward axis, base toward tip -- matches
        `constructive.compaction._axis_dir_from_rung`'s direction-only
        recovery, replicated here rather than importing that private
        symbol). `None` if `branch.start` and `branch.end` coincide (a
        degenerate zero-length rung; skip this candidate).
    """
    start, end = branch.start, branch.end
    dx, dy = x[end] - x[start], y[end] - y[start]
    length = math.hypot(dx, dy)
    if length == 0.0:
        return None
    perp_unit = (dx / length, dy / length)
    axis_dir = (perp_unit[1], -perp_unit[0])
    return (x[start], y[start]), axis_dir


def _reflect_range(
    x: list[float], y: list[float], lo: int, hi: int, pivot: Point, axis_dir: Point
) -> None:
    """Reflect the coordinate slice `[lo, hi]` across the line through
    `pivot` along `axis_dir`, in place.

    Args:
        x: Nucleotide x-coordinates, mutated in place on `[lo, hi]` only.
        y: Nucleotide y-coordinates, mutated in place on `[lo, hi]` only.
        lo: First index of the slice (inclusive).
        hi: Last index of the slice (inclusive).
        pivot: A point ON the reflection line (stays fixed).
        axis_dir: Unit vector along the reflection line.
    """
    px, py = pivot
    ax, ay = axis_dir
    for k in range(lo, hi + 1):
        relx, rely = x[k] - px, y[k] - py
        along = relx * ax + rely * ay
        perp = relx * (-ay) + rely * ax
        x[k] = px + along * ax - perp * (-ay)
        y[k] = py + along * ay - perp * ax


def _candidate_branches(tree: StructureTree, crossing_stems: list[Stem]) -> list[Branch]:
    """The MOVE SET: the LCA-level branches separating each crossing pair.

    Only a branch containing exactly ONE of a crossing pair's endpoints can
    change that pair's distance (a branch containing BOTH endpoints only
    undergoes an internal isometry -- the intra-branch distance is
    unchanged). `lca_loop_and_branches` gives exactly the two highest-
    leverage branches that separate a pair's endpoints.

    Args:
        tree: `build_structure_tree(base_pair_map)`.
        crossing_stems: Every crossing stem.

    Returns:
        Every distinct (`closing_pair`-deduplicated) LCA-level branch
        across all crossing pairs.
    """
    seen: dict[tuple[int, int], Branch] = {}
    for stem in crossing_stems:
        for a, b in stem.pairs():
            _, branch_a, branch_b = lca_loop_and_branches(tree, a, b)
            for branch in (branch_a, branch_b):
                if branch is not None:
                    seen[branch.closing_pair] = branch
    return list(seen.values())


def _flip_candidate(
    x: list[float], y: list[float], branch: Branch
) -> tuple[list[float], list[float]] | None:
    """Reflect `branch`'s owned slice on fresh coordinate copies.

    Args:
        x: Current nucleotide x-coordinates (never mutated).
        y: Current nucleotide y-coordinates (never mutated).
        branch: The candidate branch to flip.

    Returns:
        `(trial_x, trial_y)`, fresh copies with `branch`'s slice reflected,
        or `None` if `branch`'s frame is degenerate (see `_branch_frame`).
    """
    frame = _branch_frame(x, y, branch)
    if frame is None:
        return None
    trial_x, trial_y = list(x), list(y)
    _reflect_range(trial_x, trial_y, branch.start, branch.end, *frame)
    return trial_x, trial_y


def _best_improving_flip(
    x: list[float],
    y: list[float],
    candidates: list[Branch],
    base_pair_map: list[int],
    crossing_stems: list[Stem],
    params: OverlapParams,
    base_obj: float,
) -> tuple[list[float], list[float], float] | None:
    """The best single flip that strictly improves the objective and stays clean.

    Args:
        x: Current nucleotide x-coordinates.
        y: Current nucleotide y-coordinates.
        candidates: The move set (`_candidate_branches`), each re-evaluated
            fresh from `x`/`y` (an ancestor's earlier accepted flip moves a
            descendant's own attachment frame).
        base_pair_map: The nested layout's own pair_map (checker input;
            unaffected by reflection, an isometry never changes which
            primitives are drawn).
        crossing_stems: Every crossing stem (objective input).
        params: Geometry to gate against.
        base_obj: `crossing_endpoint_distance(x, y, crossing_stems)`.

    Returns:
        `(new_x, new_y, new_obj)` for the best strictly-improving,
        checker-clean flip, or `None` if no candidate both passes
        `check_overlaps` and strictly lowers `base_obj`.
    """
    best: tuple[list[float], list[float], float] | None = None
    for branch in candidates:
        trial = _flip_candidate(x, y, branch)
        if trial is None:
            continue
        trial_x, trial_y = trial
        trial_obj = crossing_endpoint_distance(trial_x, trial_y, crossing_stems)
        if trial_obj >= base_obj or (best is not None and trial_obj >= best[2]):
            continue
        if not check_overlaps(trial_x, trial_y, base_pair_map, params).passed:
            continue
        best = (trial_x, trial_y, trial_obj)
    return best


def bias_crossing_proximity(
    x: list[float],
    y: list[float],
    base_pair_map: list[int],
    crossing_stems: list[Stem],
    params: OverlapParams,
) -> tuple[list[float], list[float]]:
    """Greedily reflect nested subtrees to bring crossing partners closer.

    Never-worse contract: the entry guard requires the input already be
    checker-clean, and every accepted move is re-verified clean by the
    frozen `check_overlaps` on the WHOLE nested layout AND strictly lowers
    `crossing_endpoint_distance` -- so the result is always checker-clean
    (when the input was) with an objective `<=` the input's. A safe no-op
    (returns `x, y` unchanged) whenever the input is dirty, has no
    candidates (pk-free input, or every crossing endpoint is a bare LCA-
    loop member), or no candidate flip both passes the checker and
    improves.

    Args:
        x: Nested layout x-coordinates.
        y: Nested layout y-coordinates.
        base_pair_map: The nested layout's own pair_map.
        crossing_stems: Every crossing stem removed by
            `extraction.max_nested_subset`.
        params: Geometry to gate against (the SAME radius the nested
            layout resolved and `place_crossings`/the renderer will use).

    Returns:
        `(x, y)`, possibly improved (never worse, never unclean).
    """
    if not check_overlaps(x, y, base_pair_map, params).passed:
        return x, y
    tree = build_structure_tree(base_pair_map)
    candidates = _candidate_branches(tree, crossing_stems)
    if not candidates:
        return x, y
    cur_x, cur_y = x, y
    obj = crossing_endpoint_distance(cur_x, cur_y, crossing_stems)
    for _ in range(2 * len(candidates)):
        best = _best_improving_flip(
            cur_x, cur_y, candidates, base_pair_map, crossing_stems, params, obj
        )
        if best is None:
            break
        cur_x, cur_y, obj = best
    return cur_x, cur_y


__all__ = ["bias_crossing_proximity", "crossing_endpoint_distance"]
