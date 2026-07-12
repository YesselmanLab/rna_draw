"""Phase 3 PK-A candidate-route geometry: the first three tiers of
`placement._route_pair`'s quality ladder (see that module) for routing a
non-overlapping line from nucleotide `i` to `j` when a straight in-plane
connector (PK-B) is rejected. Tier 3, the GUARANTEED FLOOR, lives in
`.floor` (used here too, by tier 2's exits).

Tier 0 (`find_direct_route`): the bare `[p_i, p_j]` chord. The nicest
possible route -- most pairs' direct chord is already checker-clean.

Tier 1 (`find_midpoint_bow_route`): a single-apex "elbow" out from the
`(i, j)` chord's own midpoint, escalating BOW MAGNITUDE (relative to the
layout's own extent) and DIRECTION -- cheap and usually sufficient for
short/local crossings (the common case; corpus crossing stems are short,
see the plan's Corpus facts). Kept at v1's own full magnitude range (see
`_BOW_FACTORS`); the floor tier is purely ADDITIVE coverage on top.

Tier 2 (`find_ring_route`): route each endpoint independently out to a
LOCAL ring (usually far smaller than the layout's own enclosing circle) --
v1's own cheap single-leg exit first, `floor.escape_to_ring`'s two-stage
exit as the fallback -- then around the shorter arc to the other
endpoint's own exit point. A purely local sidestep, tried at small ring
scales before paying for the floor's own always-safe enclosing-circle-
scale ring.

Every tier generates AND checks candidates (via `validate`), returning the
first clean `RoutedLine` or `None`; `placement.py` only orchestrates which
tier to try, in what order, and what to do once every tier fails.
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
# around the away-from-center one at each magnitude. MEASURED: trimming
# the largest factor (12.8) regressed real corpus structures whose only
# clean route was specifically that magnitude -- the floor tier does not
# yet reliably replace every such case (see `.floor`'s module docstring
# and this milestone's STOP-criterion note), so this ladder is kept at
# v1's own full range; tiers 0 (direct) and 3/4 (ring/floor) are the new,
# purely ADDITIVE coverage layered on top.
_BOW_FACTORS: tuple[float, ...] = (0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8)
_BOW_DIRECTION_COUNT = 16

# Tier 2 (ring route): local single-leg-exit bearing resolution (v1's own
# exit search, kept as the cheap first try -- see `_single_leg_exit`), and
# radius escalation multipliers on the layout's own enclosing-circle
# radius (see `enclosing_circle`). Below 1.0 the ring is smaller than the
# true enclosing circle -- not provably clear by construction, but still
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


def find_direct_route(
    i: int,
    j: int,
    p_i: Point,
    p_j: Point,
    params: OverlapParams,
    base_primitives: list[Primitive],
    committed_lines: list[Capsule],
    pair_map: list[int],
) -> RoutedLine | None:
    """Tier 0 (best case): the single straight `i -> j` segment, unmodified.

    The nicest possible route -- tried first so most pairs (whose direct
    chord happens to already be checker-clean) never pay for a bow or a
    ring detour at all.
    """
    return _accept_if_clean(i, j, [p_i, p_j], params, base_primitives, committed_lines, pair_map)


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
    """Tier 2: route `i -> ring -> j` at an escalating LOCAL ring radius.

    Each radius independently exits `i` and `j`: `_single_leg_exit` (v1's
    own cheap single-straight-leg search) is tried FIRST so tier 2 never
    regresses below v1's own reach; `floor.escape_to_ring` (a short local
    escape, then a straight leg to the ring -- decoupling local clearance
    from the ring's own size) is the fallback, catching endpoints the
    single leg cannot clear. Either way the WHOLE assembled path --
    including the arc -- is verified against both `base_primitives` and
    `committed_lines`.

    """
    for scale in _RING_SCALES:
        radius = base_radius * scale
        exit_i = _endpoint_exit(
            i, p_i, center, radius, params, base_primitives, committed_lines, pair_map
        )
        exit_j = _endpoint_exit(
            j, p_j, center, radius, params, base_primitives, committed_lines, pair_map
        )
        if exit_i is None or exit_j is None:
            continue
        points = _assemble_two_stage_route(center, radius, exit_i, exit_j)
        line = _accept_if_clean(i, j, points, params, base_primitives, committed_lines, pair_map)
        if line is not None:
            return line
    return None


def _endpoint_exit(
    real_index: int,
    point: Point,
    center: Point,
    radius: float,
    params: OverlapParams,
    base_primitives: list[Primitive],
    committed_lines: list[Capsule],
    pair_map: list[int],
) -> list[Point] | None:
    """`_single_leg_exit` (cheap, v1-equivalent), else `floor.escape_to_ring`.

    Local import of `.floor` breaks the routing<->floor module cycle
    (`floor` also needs `routing`'s own bearing/ring helpers); same
    pattern as `engine._layout_nested`'s local import.
    """
    exit_path = _single_leg_exit(
        real_index, point, center, radius, params, base_primitives, pair_map
    )
    if exit_path is not None:
        return exit_path

    from . import floor

    return floor.escape_to_ring(
        real_index, point, center, radius, params, base_primitives, committed_lines, pair_map
    )


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


def _ring_point(center: Point, radius: float, angle: float) -> Point:
    """A point on the ring of `radius` around `center` at `angle`."""
    return (center[0] + radius * cos(angle), center[1] + radius * sin(angle))


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


def _single_leg_exit(
    real_index: int,
    point: Point,
    center: Point,
    radius: float,
    params: OverlapParams,
    base_primitives: list[Primitive],
    pair_map: list[int],
) -> list[Point] | None:
    """v1's own exit search: the first clean SINGLE straight leg out to the ring.

    Checked over its FULL length, so it can never clear a deeply embedded
    endpoint (see `.floor`'s module docstring) -- but for the common
    non-embedded case it is the cheapest exit, so `find_ring_route` tries
    this first and only falls back to `floor.escape_to_ring` when it fails.
    """
    for angle in _exit_angle_candidates(_radial_angle(center, point)):
        exit_point = _ring_point(center, radius, angle)
        leg = _leg_capsule(point, exit_point, real_index, params.pair_half_width)
        if capsule_is_clean(leg, base_primitives, pair_map, params.tol):
            return [point, exit_point]
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


def _assemble_two_stage_route(
    center: Point, radius: float, exit_i: list[Point], exit_j: list[Point]
) -> list[Point]:
    """Join two endpoints' own two-stage exit paths via the shorter ring arc.

    `exit_i`/`exit_j` are `[point, ..., ring_point]` waypoint lists (from
    `floor.escape_to_ring`/`floor.creep_out`); both end exactly on the
    `radius` ring, so the arc between their two ring points is well-formed
    regardless of how many interior waypoints either exit path has.
    """
    angle_i = _radial_angle(center, exit_i[-1])
    angle_j = _radial_angle(center, exit_j[-1])
    arc = _arc_waypoints(center, radius, angle_i, angle_j)
    return [*exit_i, *arc, *reversed(exit_j)]


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


__all__ = [
    "enclosing_circle",
    "find_direct_route",
    "find_midpoint_bow_route",
    "find_ring_route",
]
