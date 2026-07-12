"""Phase 3 PK-A GUARANTEED-FLOOR routing (tier 4 of `placement._route_pair`'s
quality ladder, see `routing`'s module docstring for tiers 0-2).

The floor exists because a single straight leg checked over its FULL
length -- `routing`'s old `_find_exit_angle` -- can never clear a deeply
embedded endpoint: the blocked-bearing arcs of its surrounding primitives
sum past 2*pi at every scale (small ring: the ARC cuts the interior; large
ring: the EXIT LEG is blocked end-to-end). This module decouples the two
concerns instead: a SHORT local-escape leg (independent of the ring
radius) gets the endpoint clear of its immediate neighbors, then a
separate straight leg carries it the rest of the way out to the ring
(`escape_to_ring`). An endpoint too tightly boxed in for even that falls
back to `creep_out`: a bounded, greedy, monotonically-outward walk that
takes one step at a time in whichever clean bearing has the most room.

`find_floor_route` assembles `escape(i)` (or `creep_out(i)`) + the shorter
ring arc + reversed `escape(j)` (or `creep_out(j)`) into one candidate
polyline and validates the WHOLE thing via `routing._accept_if_clean` --
the same frozen-checker-backed predicate every other tier uses. Nothing
here ever loosens or reimplements that predicate; `_clearance`/
`_min_clearance` below only RANK already-clean candidates (by how much
room they leave), using the frozen `geometry` distance primitives, never
the overlap decision itself.
"""

from __future__ import annotations

from math import cos, hypot, sin

from rna_draw.geometry import (
    Capsule,
    Disk,
    PrimitiveId,
    point_segment_distance,
    segment_segment_distance,
)
from rna_draw.layout.base import RoutedLine
from rna_draw.overlap import OverlapParams, Primitive, is_excluded

from . import routing
from .validate import capsule_is_clean

Point = tuple[float, float]

# Local-escape length ladder, in units of one nucleotide "footprint" (disk
# diameter plus a pair-capsule width) -- short by design and independent
# of the ring radius; see the module docstring's "decouple" rationale.
_ESCAPE_LENGTH_FACTORS: tuple[float, ...] = (1.0, 2.0, 4.0, 8.0, 16.0, 32.0)
_ESCAPE_BEARING_STEPS = 32

# `creep_out`'s step-length ladder (same unit), bearing resolution, and
# bounded iteration cap -- a near-guarantee, not a closed proof for
# finite-width legs (see the plan's Risks/STOP section): kissing disks
# could in principle seal a gap narrower than the stroke. Escalating
# lengths (rather than one fixed step) matters: a small step almost
# always finds SOME clean local direction, but only a larger one can
# clear an entire nearby cluster in a single hop -- see `_creep_step`.
_CREEP_LENGTH_FACTORS: tuple[float, ...] = (1.0, 4.0, 16.0)
_CREEP_BEARING_STEPS = 16
_CREEP_MAX_STEPS = 10
# How many of the walk's own last waypoints a new step must clear (in
# units of `_footprint_unit`) -- breaks a two-point oscillation between
# similarly-good candidates that would otherwise never make net progress.
_CREEP_AVOID_RECENT = 4
_CREEP_AVOID_RADIUS_FACTOR = 0.5


def _footprint_unit(params: OverlapParams) -> float:
    """One nucleotide "footprint" length: the local-escape/creep step unit."""
    return 2.0 * (params.node_r + params.pair_half_width)


def escape_to_ring(
    real_index: int,
    point: Point,
    center: Point,
    ring_radius: float,
    params: OverlapParams,
    base_primitives: list[Primitive],
    committed_lines: list[Capsule],
    pair_map: list[int],
) -> list[Point] | None:
    """Two-stage exit: a short local escape, then a straight leg out to the ring.

    Escalates the escape leg's length over `_ESCAPE_LENGTH_FACTORS`; at
    each length, `_best_clear_bearing` picks the clean local direction
    with the most room, then a straight leg from THAT point radially out
    to the `radius`-ring is tried (`_snap_to_ring`). Succeeding at a short
    length keeps the overall route short; escalating covers an endpoint
    whose immediate neighbors block every short escape.

    Args:
        real_index: The endpoint's own nucleotide index (excluded from its
            own leg's cleanliness check, matching `routing`'s other tiers).
        point: The endpoint's layout position.
        center: The layout's bounding-box center (shared with `routing`).
        ring_radius: Radius of the ring this exit must reach.
        params: Geometry to gate against.
        base_primitives: Every disk/backbone/pair capsule the layout draws.
        committed_lines: Every already-accepted PK-A line's own capsules.
        pair_map: The current augmented pair_map.

    Returns:
        `[point, escape_point, ring_point]` if some escalation step
        cleared both stages, else `None` (the caller falls back to
        `creep_out`).
    """
    base_bearing = routing._radial_angle(center, point)
    unit = _footprint_unit(params)
    for factor in _ESCAPE_LENGTH_FACTORS:
        escape_point = _best_clear_bearing(
            point,
            factor * unit,
            base_bearing,
            _ESCAPE_BEARING_STEPS,
            frozenset({real_index}),
            params,
            base_primitives,
            committed_lines,
            pair_map,
        )
        if escape_point is None:
            continue
        exit_point = _snap_to_ring(
            escape_point, center, ring_radius, params, base_primitives, committed_lines, pair_map
        )
        if exit_point is not None:
            return [point, escape_point, exit_point]
    return None


def creep_out(
    real_index: int,
    point: Point,
    center: Point,
    ring_radius: float,
    params: OverlapParams,
    base_primitives: list[Primitive],
    committed_lines: list[Capsule],
    pair_map: list[int],
) -> list[Point] | None:
    """Fallback for a genuinely surrounded endpoint: a greedy outward walk.

    At each waypoint, first tries to snap straight out to the ring
    (`_snap_to_ring`; usually succeeds as soon as the walk clears its
    local obstacles); otherwise takes one step toward whichever CLEAN
    candidate -- across every length in `_CREEP_LENGTH_FACTORS` and every
    fanned bearing -- ends up farthest from `center` (`_creep_step`), also
    steering clear of the walk's own last few waypoints so two
    similarly-good candidates cannot trap it in a back-and-forth
    oscillation. The bounded `_CREEP_MAX_STEPS` cap is what makes this a
    near-guarantee rather than a closed proof (see the plan's Risks/STOP
    section): most real packings have enough positive-clearance local room
    to escape well within the cap, but a pathological one could still
    exhaust it.

    Args:
        real_index: The endpoint's own nucleotide index.
        point: The endpoint's layout position.
        center: The layout's bounding-box center.
        ring_radius: Radius of the ring this walk must reach.
        params: Geometry to gate against.
        base_primitives: Every disk/backbone/pair capsule the layout draws.
        committed_lines: Every already-accepted PK-A line's own capsules.
        pair_map: The current augmented pair_map.

    Returns:
        The walk's waypoints (`point` through a ring-radius exit point),
        or `None` if `_CREEP_MAX_STEPS` was exhausted without a clean
        step -- recorded honestly by the caller, never a silent overlap.
    """
    unit = _footprint_unit(params)
    min_first_length = params.node_r + params.pair_half_width + params.tol
    waypoints = [point]
    recent = [point]
    current = point
    ends = frozenset({real_index})
    for _ in range(_CREEP_MAX_STEPS):
        exit_point = _snap_to_ring(
            current, center, ring_radius, params, base_primitives, committed_lines, pair_map
        )
        if exit_point is not None:
            return [*waypoints, exit_point]
        next_point = _creep_step(
            current,
            center,
            ends,
            min_first_length,
            unit,
            recent,
            params,
            base_primitives,
            committed_lines,
            pair_map,
        )
        if next_point is None:
            return None
        waypoints.append(next_point)
        recent = [*recent, next_point][-_CREEP_AVOID_RECENT:]
        current = next_point
        ends = frozenset()  # only the walk's own first leg touches the real disk
    return None


def find_floor_route(
    i: int,
    j: int,
    p_i: Point,
    p_j: Point,
    center: Point,
    ring_radius: float,
    params: OverlapParams,
    base_primitives: list[Primitive],
    committed_lines: list[Capsule],
    pair_map: list[int],
) -> RoutedLine | None:
    """Tier 4 (GUARANTEED FLOOR): escape(i) + shorter arc + reversed escape(j).

    Args:
        i: 5' nucleotide index.
        j: 3' nucleotide index.
        p_i: `i`'s layout position.
        p_j: `j`'s layout position.
        center: The layout's bounding-box center.
        ring_radius: THIS floor line's own concentric ring radius (see
            `placement`'s per-line radius growth: strictly nested radii
            keep every floor line's arc disjoint from every other's by
            construction).
        params: Geometry to gate against.
        base_primitives: Every disk/backbone/pair capsule the layout draws.
        committed_lines: Every already-accepted PK-A line's own capsules.
        pair_map: The current augmented pair_map.

    Returns:
        A `RoutedLine` iff the WHOLE assembled polyline validates clean
        (`routing._accept_if_clean`), else `None`.
    """
    exit_i = _exit_points(
        i, p_i, center, ring_radius, params, base_primitives, committed_lines, pair_map
    )
    exit_j = _exit_points(
        j, p_j, center, ring_radius, params, base_primitives, committed_lines, pair_map
    )
    if exit_i is None or exit_j is None:
        return None
    points = routing._assemble_two_stage_route(center, ring_radius, exit_i, exit_j)
    return routing._accept_if_clean(
        i, j, points, params, base_primitives, committed_lines, pair_map
    )


def _exit_points(
    real_index: int,
    point: Point,
    center: Point,
    ring_radius: float,
    params: OverlapParams,
    base_primitives: list[Primitive],
    committed_lines: list[Capsule],
    pair_map: list[int],
) -> list[Point] | None:
    """`escape_to_ring`, falling back to `creep_out` for a surrounded endpoint."""
    exit_path = escape_to_ring(
        real_index, point, center, ring_radius, params, base_primitives, committed_lines, pair_map
    )
    if exit_path is not None:
        return exit_path
    return creep_out(
        real_index, point, center, ring_radius, params, base_primitives, committed_lines, pair_map
    )


def _snap_to_ring(
    point: Point,
    center: Point,
    ring_radius: float,
    params: OverlapParams,
    base_primitives: list[Primitive],
    committed_lines: list[Capsule],
    pair_map: list[int],
) -> Point | None:
    """A clean straight leg from `point` to the ring point at its own bearing."""
    exit_point = routing._ring_point(center, ring_radius, routing._radial_angle(center, point))
    leg = Capsule(
        pid=PrimitiveId("pkfloorexit", 0),
        x0=point[0],
        y0=point[1],
        x1=exit_point[0],
        y1=exit_point[1],
        half_width=params.pair_half_width,
        ends=frozenset(),
    )
    others: list[Primitive] = [*base_primitives, *committed_lines]
    if capsule_is_clean(leg, others, pair_map, params.tol):
        return exit_point
    return None


def _step_leg(
    point: Point, bearing: float, length: float, ends: frozenset[int], half_width: float
) -> tuple[Point, Capsule]:
    """A candidate step point at `length` along `bearing` from `point`, and its capsule."""
    candidate = (point[0] + length * cos(bearing), point[1] + length * sin(bearing))
    leg = Capsule(
        pid=PrimitiveId("pkfloorleg", 0),
        x0=point[0],
        y0=point[1],
        x1=candidate[0],
        y1=candidate[1],
        half_width=half_width,
        ends=ends,
    )
    return candidate, leg


def _best_clear_bearing(
    point: Point,
    length: float,
    base_bearing: float,
    n_bearings: int,
    ends: frozenset[int],
    params: OverlapParams,
    base_primitives: list[Primitive],
    committed_lines: list[Capsule],
    pair_map: list[int],
) -> Point | None:
    """The clean `length`-step (if any) with the greatest obstacle clearance."""
    others: list[Primitive] = [*base_primitives, *committed_lines]
    best_point: Point | None = None
    best_clearance = -1.0
    for bearing in routing._fan_out(base_bearing, n_bearings):
        candidate, leg = _step_leg(point, bearing, length, ends, params.pair_half_width)
        if not capsule_is_clean(leg, others, pair_map, params.tol):
            continue
        clearance = _min_clearance(leg, others, pair_map)
        if clearance > best_clearance:
            best_point, best_clearance = candidate, clearance
    return best_point


def _creep_step(
    current: Point,
    center: Point,
    ends: frozenset[int],
    min_first_length: float,
    unit: float,
    recent: list[Point],
    params: OverlapParams,
    base_primitives: list[Primitive],
    committed_lines: list[Capsule],
    pair_map: list[int],
) -> Point | None:
    """The globally best (farthest from `center`) clean step across every length.

    Tries EVERY length in `_CREEP_LENGTH_FACTORS`, not just the first that
    works: a small step almost always finds some clean local direction,
    but only a larger one can clear an entire nearby cluster in a single
    hop, so ranking candidates from every length together (instead of
    stopping at the first length with any success) is what lets the walk
    escape a tight pocket instead of shuffling in small, locally-optimal
    but globally-stuck steps.
    """
    others: list[Primitive] = [*base_primitives, *committed_lines]
    base_bearing = routing._radial_angle(center, current)
    best_point: Point | None = None
    best_radius = -1.0
    for factor in _CREEP_LENGTH_FACTORS:
        length = factor * unit
        if ends and length < min_first_length:
            continue
        for bearing in routing._fan_out(base_bearing, _CREEP_BEARING_STEPS):
            candidate, leg = _step_leg(current, bearing, length, ends, params.pair_half_width)
            if _too_close_to_recent(candidate, recent, unit):
                continue
            if not capsule_is_clean(leg, others, pair_map, params.tol):
                continue
            radius = hypot(candidate[0] - center[0], candidate[1] - center[1])
            if radius > best_radius:
                best_point, best_radius = candidate, radius
    return best_point


def _too_close_to_recent(point: Point, recent: list[Point], unit: float) -> bool:
    """Whether `point` sits within one recent-avoidance radius of a just-visited waypoint.

    Breaks a two-point oscillation: without this, two nearly-tied
    candidates on opposite sides of a pinch point can keep re-electing
    each other forever (bounded only by `_CREEP_MAX_STEPS`, never
    escaping).
    """
    threshold = _CREEP_AVOID_RADIUS_FACTOR * unit
    return any(hypot(point[0] - r[0], point[1] - r[1]) < threshold for r in recent)


def _clearance(leg: Capsule, other: Primitive) -> float:
    """How far `leg` sits past `other`'s required touching distance.

    Positive: clear by that many layout units. Negative: `other` would
    flag as an overlap by that much. Built ONLY from the frozen
    `geometry` distance primitives (`point_segment_distance`/
    `segment_segment_distance`); the overlap ACCEPT/REJECT decision
    itself stays `capsule_is_clean`/`is_excluded` everywhere -- this is
    purely a tie-breaker for ranking several already-clean candidates.
    """
    if isinstance(other, Disk):
        dist = point_segment_distance(other.cx, other.cy, leg.x0, leg.y0, leg.x1, leg.y1)
        return dist - (other.radius + leg.half_width)
    dist = segment_segment_distance(
        (leg.x0, leg.y0), (leg.x1, leg.y1), (other.x0, other.y0), (other.x1, other.y1)
    )
    return dist - (other.half_width + leg.half_width)


def _min_clearance(leg: Capsule, candidates: list[Primitive], pair_map: list[int]) -> float:
    """The tightest clearance from `leg` to any non-excluded `candidates` member."""
    margins = [
        _clearance(leg, other) for other in candidates if not is_excluded(leg, other, pair_map)
    ]
    return min(margins) if margins else float("inf")


__all__ = ["escape_to_ring", "creep_out", "find_floor_route"]
