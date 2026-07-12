"""Phase 3: place crossing stems as in-plane PK-B connectors or PK-A routed
lines, escalating per crossing stem, checker-gated the whole way.

KEY GEOMETRIC DECISION (locked, see the plan): crossing nucleotides are
NEVER re-placed -- they already sit where the nested layout (Phase 2) put
them. A crossing stem is purely an ADDITIONAL connector between two
already-placed nucleotides: PK-B (`_try_in_plane`) adds it to `pair_map`
and re-checks the frozen `check_overlaps` unchanged; PK-A (`_route_stem`)
routes a non-overlapping polyline per pair instead, via `routing`'s two
escalating strategies, validated by `validate.polyline_is_clean` against
the frozen geometry predicates.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from rna_draw.geometry import Capsule
from rna_draw.layout.base import RoutedLine
from rna_draw.overlap import OverlapParams, Primitive, build_primitives, check_overlaps

from . import routing
from .parsing import Stem
from .validate import polyline_capsules

Point = tuple[float, float]

# A crossing stem this long (bp count) is never kept as a straight in-plane
# PK-B connector even if geometrically clean (MUST-FIX #4): v1 ships with
# no proximity bias (Phase 2b is deferred), so a long straight chord can be
# checker-clean yet still slash conventionally-unwanted distances across
# the figure. 6 matches the measured corpus p90 crossing-stem length (plan
# `Corpus facts`, 316 real pseudoknots); longer crossing stems always route
# as PK-A lines. Tunable.
MAX_PK_B_LENGTH = 6


@dataclass
class PlacementResult:
    """The outcome of placing every crossing stem.

    Args:
        pk_b_pairs: Flattened pairs from every crossing stem placed as a
            clean in-plane connector.
        pk_a_lines: One `RoutedLine` per pair from every crossing stem
            that escalated to a routed line.
        unplaced: Crossing stems for which at least one pair's PK-A search
            never cleared (both `routing` tiers exhausted; recorded
            honestly if it ever happens -- never a silent overlap. Real
            corpus dense/large structures can genuinely have no
            straight-segment route between two deeply embedded
            nucleotides; the nested layout stays clean either way).
    """

    pk_b_pairs: list[tuple[int, int]] = field(default_factory=list)
    pk_a_lines: list[RoutedLine] = field(default_factory=list)
    unplaced: list[Stem] = field(default_factory=list)


@dataclass
class _PlacementState:
    """Mutable accumulators threaded through the crossing-placement escalation.

    Args:
        x: Nested layout's nucleotide x-coordinates (never mutated).
        y: Nested layout's nucleotide y-coordinates (never mutated).
        params: Geometry to gate against (the nested layout's own radius).
        pair_map: The current augmented pair_map; grows by one stem's
            pairs each time a crossing stem is accepted as PK-B.
        committed_lines: Every already-accepted PK-A line's own capsules,
            checked against so later lines never cross earlier ones.
        center: The layout's own bounding-box center (`routing`'s shared
            reference point for both tiers).
        extent: The layout's own bounding-box diagonal (tier 1's bow-
            magnitude unit).
        base_radius: The layout's enclosing-circle radius (tier 2's ring
            unit); see `routing.enclosing_circle`.
        result: The `PlacementResult` being built.
    """

    x: list[float]
    y: list[float]
    params: OverlapParams
    pair_map: list[int]
    committed_lines: list[Capsule] = field(default_factory=list)
    center: Point = (0.0, 0.0)
    extent: float = 1.0
    base_radius: float = 1.0
    result: PlacementResult = field(default_factory=PlacementResult)


def place_crossings(
    x: list[float],
    y: list[float],
    base_pair_map: list[int],
    crossing_stems: list[Stem],
    params: OverlapParams,
) -> PlacementResult:
    """Escalate every crossing stem: in-plane PK-B first, else routed PK-A lines.

    Args:
        x: Nucleotide x-coordinates from the nested layout (never moved).
        y: Nucleotide y-coordinates from the nested layout.
        base_pair_map: The nested layout's own pair_map (retained pairs).
        crossing_stems: Stems removed by `extraction.max_nested_subset`.
        params: Geometry to gate against -- the SAME radius the nested
            layout resolved and the renderer will draw at.

    Returns:
        The accumulated `PlacementResult`.
    """
    center, base_radius = routing.enclosing_circle(x, y, params)
    state = _PlacementState(
        x=x,
        y=y,
        params=params,
        pair_map=list(base_pair_map),
        center=center,
        extent=max(_extent(x, y), params.node_r),
        base_radius=base_radius,
    )
    for stem in crossing_stems:
        _place_one_stem(state, stem)
    return state.result


def _extent(x: list[float], y: list[float]) -> float:
    """The layout's own bounding-box diagonal length."""
    width = max(x) - min(x)
    height = max(y) - min(y)
    return (width**2 + height**2) ** 0.5


def _place_one_stem(state: _PlacementState, stem: Stem) -> None:
    """PK-B first; PK-A per-pair fallback if PK-B is rejected or too long."""
    if _try_in_plane(state, stem):
        state.pair_map = _augment(state.pair_map, stem.pairs())
        state.result.pk_b_pairs.extend(stem.pairs())
        return
    _route_stem(state, stem)


def _try_in_plane(state: _PlacementState, stem: Stem) -> bool:
    """PK-B: whether adding `stem`'s pairs to `pair_map` stays checker-clean."""
    if stem.length > MAX_PK_B_LENGTH:
        return False
    candidate = _augment(state.pair_map, stem.pairs())
    return check_overlaps(state.x, state.y, candidate, state.params).passed


def _augment(pair_map: list[int], pairs: list[tuple[int, int]]) -> list[int]:
    """Return `pair_map` with `pairs` added (both directions)."""
    augmented = list(pair_map)
    for i, j in pairs:
        augmented[i] = j
        augmented[j] = i
    return augmented


def _route_stem(state: _PlacementState, stem: Stem) -> None:
    """Route every pair of a PK-B-rejected crossing stem as its own PK-A line."""
    base_primitives = build_primitives(state.x, state.y, state.pair_map, state.params)
    any_unplaced = False
    for i, j in stem.pairs():
        line = _route_pair(state, i, j, base_primitives)
        if line is None:
            any_unplaced = True
            continue
        state.result.pk_a_lines.append(line)
        line_uid = len(state.result.pk_a_lines)
        state.committed_lines.extend(
            polyline_capsules(line.points, i, j, state.params.pair_half_width, line_uid)
        )
    if any_unplaced:
        state.result.unplaced.append(stem)


def _route_pair(
    state: _PlacementState, i: int, j: int, base_primitives: list[Primitive]
) -> RoutedLine | None:
    """Try `routing`'s tier 1 (cheap midpoint bow), else tier 2 (ring route)."""
    p_i, p_j = (state.x[i], state.y[i]), (state.x[j], state.y[j])
    line = routing.find_midpoint_bow_route(
        i,
        j,
        p_i,
        p_j,
        state.center,
        state.extent,
        state.params,
        base_primitives,
        state.committed_lines,
        state.pair_map,
    )
    if line is not None:
        return line
    return routing.find_ring_route(
        i,
        j,
        p_i,
        p_j,
        state.center,
        state.base_radius,
        state.params,
        base_primitives,
        state.committed_lines,
        state.pair_map,
    )


__all__ = ["PlacementResult", "place_crossings", "MAX_PK_B_LENGTH"]
