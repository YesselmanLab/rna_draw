"""Phase 3 PK-A candidate-route geometry: two escalating strategies for
routing a non-overlapping line from nucleotide `i` to `j` when a straight
in-plane connector (PK-B) is rejected.

Tier 1 (`find_midpoint_bow_route`): a single-apex "elbow" out from the
`(i, j)` chord's own midpoint, escalating BOW MAGNITUDE (relative to the
layout's own extent) and DIRECTION -- cheap and usually sufficient for
short/local crossings (the common case; corpus crossing stems are short,
see the plan's Corpus facts).

Tier 2 (`find_ring_route`): route each endpoint independently out to a
circle that provably encloses the WHOLE layout (so it is clear of every
real primitive by construction), then around the shorter arc to the other
endpoint's own exit point. More robust for a nucleotide deeply embedded in
a dense structure, where tier 1's fixed-midpoint bow direction can keep
clipping a close neighbor regardless of magnitude.

Both tiers generate AND check candidates (via `validate`), returning the
first clean `RoutedLine` or `None`; `placement.py` only orchestrates which
tier to try, in what order, and what to do once both fail.
"""

from __future__ import annotations

from math import atan2, cos, hypot, pi, sin

from rna_draw.geometry import Capsule, PrimitiveId
from rna_draw.layout.base import RoutedLine
from rna_draw.overlap import OverlapParams, Primitive

from .validate import capsule_is_clean, polyline_capsules, polyline_is_clean

Point = tuple[float, float]

# Tier 1 (midpoint bow): magnitude ladder as fractions/multiples of the
# layout's own bounding-box diagonal, and how many bearings to fan out
# around the away-from-center one at each magnitude.
_BOW_FACTORS: tuple[float, ...] = (0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8)
_BOW_DIRECTION_COUNT = 16

# Tier 2 (ring route): local-exit bearing resolution, and radius
# escalation multipliers on the layout's own enclosing-circle radius (see
# `enclosing_circle`). Below 1.0 the ring is smaller than the true
# enclosing circle -- not provably clear by construction, but still
# checked, and often succeeds as a purely local sidestep far short of a
# full trip out to the guaranteed-safe floor.
_EXIT_ANGLE_STEPS = 24
_RING_SCALES: tuple[float, ...] = (0.03, 0.06, 0.1, 0.2, 0.35, 0.55, 1.0, 1.3, 1.7, 2.2, 3.0, 4.0)
_ARC_STEP = pi / 9  # 20 degrees
_OUTER_MARGIN_FACTOR = 1.5


def enclosing_circle(x: list[float], y: list[float], params: OverlapParams) -> tuple[Point, float]:
    """A circle, centered on the layout's bbox center, that encloses every
    primitive the layout can draw -- points outside it are provably clear
    of every disk/backbone/pair capsule.
    """
    center = ((max(x) + min(x)) / 2, (max(y) + min(y)) / 2)
    half_diag = hypot(max(x) - center[0], max(y) - center[1])
    margin = _OUTER_MARGIN_FACTOR * (
        params.node_r + max(params.backbone_half_width, params.pair_half_width)
    )
    return center, max(half_diag + margin, params.node_r)


def find_midpoint_bow_route(
    i: int,
    j: int,
    p_i: Point,
    p_j: Point,
    center: Point,
    extent: float,
    params: OverlapParams,
    base_primitives: list[Primitive],
    committed_lines: list[Capsule],
    pair_map: list[int],
) -> RoutedLine | None:
    """Tier 1: search (magnitude, bearing) for a clean single-apex elbow.

    Magnitude alone cannot guarantee clearing: a fixed-direction elbow's
    approach at each endpoint stays the same bearing as the apex recedes,
    so a single "always bow away from center" direction can keep clipping
    a neighbor sitting close to `i` or `j` no matter how large the bow
    gets -- hence the bearing fan-out at every magnitude.
    """
    base_angle = _away_angle(center, p_i, p_j)
    bearings = _fan_out(base_angle, _BOW_DIRECTION_COUNT)
    for factor in _BOW_FACTORS:
        magnitude = factor * extent
        for bearing in bearings:
            direction = (cos(bearing), sin(bearing))
            points = _elbow_points(p_i, p_j, direction, magnitude)
            line = _accept_if_clean(
                i, j, points, params, base_primitives, committed_lines, pair_map
            )
            if line is not None:
                return line
    return None


def find_ring_route(
    i: int,
    j: int,
    p_i: Point,
    p_j: Point,
    center: Point,
    base_radius: float,
    params: OverlapParams,
    base_primitives: list[Primitive],
    committed_lines: list[Capsule],
    pair_map: list[int],
) -> RoutedLine | None:
    """Tier 2: route `i -> ring -> j` at an escalating ring radius.

    Each radius independently searches a local exit bearing for `i` and
    for `j` (each checked only against `base_primitives`, since the arc is
    provably clear of them by construction at scale >= 1.0), then verifies
    the WHOLE assembled path -- including the arc -- against both
    `base_primitives` and `committed_lines`.
    """
    for scale in _RING_SCALES:
        radius = base_radius * scale
        angle_i = _find_exit_angle(i, p_i, center, radius, params, base_primitives, pair_map)
        angle_j = _find_exit_angle(j, p_j, center, radius, params, base_primitives, pair_map)
        if angle_i is None or angle_j is None:
            continue
        points = _assemble_ring_route(center, radius, p_i, angle_i, p_j, angle_j)
        line = _accept_if_clean(i, j, points, params, base_primitives, committed_lines, pair_map)
        if line is not None:
            return line
    return None


def _accept_if_clean(
    i: int,
    j: int,
    points: list[Point],
    params: OverlapParams,
    base_primitives: list[Primitive],
    committed_lines: list[Capsule],
    pair_map: list[int],
) -> RoutedLine | None:
    """Build `points`' capsules and return a `RoutedLine` iff they're clean."""
    segments = polyline_capsules(points, i, j, params.pair_half_width, line_uid=0)
    if polyline_is_clean(segments, base_primitives, committed_lines, pair_map, params.tol):
        return RoutedLine(i=i, j=j, points=points)
    return None


def _away_angle(center: Point, p_i: Point, p_j: Point) -> float:
    """Bearing from `center` to the `(p_i, p_j)` chord's midpoint."""
    mid = ((p_i[0] + p_j[0]) / 2, (p_i[1] + p_j[1]) / 2)
    return _radial_angle(center, mid)


def _radial_angle(center: Point, point: Point) -> float:
    """Bearing from `center` to `point`, in radians."""
    dx, dy = point[0] - center[0], point[1] - center[1]
    if dx == 0.0 and dy == 0.0:
        return 0.0
    return atan2(dy, dx)


def _fan_out(base_angle: float, count: int) -> list[float]:
    """`count` bearings evenly spaced around `base_angle`, tried first."""
    step = 2 * pi / count
    return [base_angle + k * step for k in range(count)]


def _elbow_points(p_i: Point, p_j: Point, direction: Point, magnitude: float) -> list[Point]:
    """A 3-point outward "elbow" from `p_i` to `p_j` via one apex.

    A straight-out, straight-back elbow (rather than a sampled curve)
    deliberately avoids a smooth curve's failure mode here: a Bezier
    parametrized uniformly in `t` slows to a near-stop at its midpoint as
    the control point recedes, clustering samples together there and
    producing spurious CAPSULE self-overlap between non-adjacent segments
    of the very same (non-self-intersecting) curve. Two straight segments
    sharing one apex point have no such interior curvature to misjudge.
    """
    mid = ((p_i[0] + p_j[0]) / 2, (p_i[1] + p_j[1]) / 2)
    apex = (mid[0] + direction[0] * magnitude, mid[1] + direction[1] * magnitude)
    return [p_i, apex, p_j]


def _exit_angle_candidates(base_angle: float) -> list[float]:
    """Bearings to try for a local exit, starting at `base_angle` (the
    shortest, purely radial exit) and fanning out to cover the full circle.
    """
    step = 2 * pi / _EXIT_ANGLE_STEPS
    angles = [base_angle]
    for k in range(1, _EXIT_ANGLE_STEPS // 2 + 1):
        angles.append(base_angle + k * step)
        angles.append(base_angle - k * step)
    return angles


def _ring_point(center: Point, radius: float, angle: float) -> Point:
    """A point on the ring of `radius` around `center` at `angle`."""
    return (center[0] + radius * cos(angle), center[1] + radius * sin(angle))


def _find_exit_angle(
    real_index: int,
    point: Point,
    center: Point,
    radius: float,
    params: OverlapParams,
    base_primitives: list[Primitive],
    pair_map: list[int],
) -> float | None:
    """The first clean bearing (radial-first) from `point` out to the ring."""
    for angle in _exit_angle_candidates(_radial_angle(center, point)):
        exit_point = _ring_point(center, radius, angle)
        leg = _leg_capsule(point, exit_point, real_index, params.pair_half_width)
        if capsule_is_clean(leg, base_primitives, pair_map, params.tol):
            return angle
    return None


def _leg_capsule(point: Point, exit_point: Point, real_index: int, half_width: float) -> Capsule:
    """A capsule for one local exit leg, tagged so its real end is excluded.

    `ends={real_index}` makes `is_excluded` treat the leg's own starting
    nucleotide's disk as an intentional touch (the leg literally starts at
    that disk's center); the ring-side end has no real nucleotide, so it
    carries no id -- nothing else in the layout can coincidentally share
    an id with it.
    """
    return Capsule(
        pid=PrimitiveId("pkleg", real_index),
        x0=point[0],
        y0=point[1],
        x1=exit_point[0],
        y1=exit_point[1],
        half_width=half_width,
        ends=frozenset({real_index}),
    )


def _assemble_ring_route(
    center: Point, radius: float, p_i: Point, angle_i: float, p_j: Point, angle_j: float
) -> list[Point]:
    """Build the full `i -> exit_i -> (arc) -> exit_j -> j` waypoint list."""
    exit_i = _ring_point(center, radius, angle_i)
    exit_j = _ring_point(center, radius, angle_j)
    arc = _arc_waypoints(center, radius, angle_i, angle_j)
    return [p_i, exit_i, *arc, exit_j, p_j]


def _arc_waypoints(center: Point, radius: float, angle_a: float, angle_b: float) -> list[Point]:
    """Ring points strictly between `angle_a` and `angle_b`, shorter way round.

    A tiny epsilon keeps the LAST generated waypoint strictly short of
    `angle_b` even when `delta` is an exact multiple of `_ARC_STEP` (e.g.
    antipodal `i`/`j` exit angles) -- without it, that edge case emits a
    waypoint numerically equal to `exit_j`, producing a zero-length
    segment that spuriously self-conflicts (`polyline_is_clean`'s sibling
    check on two capsules meeting at the same point without a shared id).
    """
    delta = (angle_b - angle_a + pi) % (2 * pi) - pi
    step = _ARC_STEP if delta >= 0 else -_ARC_STEP
    n_steps = int((abs(delta) - 1e-9) // _ARC_STEP)
    return [_ring_point(center, radius, angle_a + step * k) for k in range(1, n_steps + 1)]


__all__ = ["enclosing_circle", "find_midpoint_bow_route", "find_ring_route"]
