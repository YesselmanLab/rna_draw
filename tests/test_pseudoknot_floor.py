"""Tests for `rna_draw.layout.pseudoknot.floor`: the GUARANTEED-FLOOR tier
(`escape_to_ring`/`creep_out`/`find_floor_route`) at the unit level.

`_boxed_in_rings` builds a point genuinely surrounded by two concentric,
gapped rings of disks whose gaps are offset from each other (see the
module docstring's root-cause writeup): no SINGLE straight leg -- the
short local escape stage alone, or the old single-leg exit `routing`'s
tier 2 tries first -- can clear both rings at once, so `escape_to_ring`
must fail; only a BENT path (a short escape through the inner gap, then a
turn toward the outer gap) can, which is exactly what `creep_out` walks.
"""

from __future__ import annotations

from math import atan2, cos, degrees, radians, sin

from rna_draw.geometry import Disk, PrimitiveId
from rna_draw.layout.pseudoknot import floor
from rna_draw.layout.pseudoknot.validate import polyline_capsules, polyline_is_clean
from rna_draw.overlap import OverlapParams, Primitive

Point = tuple[float, float]

PARAMS = OverlapParams()
PAIR_MAP = [-1] * 10


def _ring_disks(
    center: Point, radius: float, gap_deg: float, disk_radius: float, id_base: int
) -> list[Primitive]:
    """A ring of 16 disks around `center`, with one 40-degree gap at `gap_deg`."""
    disks: list[Primitive] = []
    for k in range(16):
        angle_deg = k * 360.0 / 16
        delta = (angle_deg - gap_deg + 180) % 360 - 180
        if abs(delta) < 20.0:
            continue
        angle = radians(angle_deg)
        cx = center[0] + radius * cos(angle)
        cy = center[1] + radius * sin(angle)
        pid = PrimitiveId("obstacle", id_base + k)
        disks.append(
            Disk(pid=pid, cx=cx, cy=cy, radius=disk_radius, ends=frozenset({-(id_base + k)}))
        )
    return disks


def _sealed_ring(
    center: Point, radius: float, disk_radius: float, n: int, id_base: int
) -> list[Primitive]:
    """A GAPLESS ring: adjacent disks overlap angularly, sealing every bearing.

    No straight leg nor bent walk can ever escape this (there is
    genuinely no positive-clearance gap) -- the pathological case the
    plan's STOP criterion documents; `escape_to_ring`/`creep_out` must
    both honestly return `None` rather than fabricate a dirty route.
    """
    disks: list[Primitive] = []
    for k in range(n):
        angle = radians(k * 360.0 / n)
        cx = center[0] + radius * cos(angle)
        cy = center[1] + radius * sin(angle)
        pid = PrimitiveId("sealed", id_base + k)
        disks.append(
            Disk(pid=pid, cx=cx, cy=cy, radius=disk_radius, ends=frozenset({-(id_base + k)}))
        )
    return disks


def _away_bearing_deg(center: Point, point: Point) -> float:
    """Bearing (degrees) from `center` to `point` -- the ring gaps' reference."""
    return degrees(atan2(point[1] - center[1], point[0] - center[0]))


def _boxed_in_rings(center: Point, point: Point, id_base: int) -> list[Primitive]:
    """Two concentric rings around `point`, gaps 90 degrees apart.

    No single straight bearing clears both: whichever bearing threads the
    inner ring's gap does not line up with the (offset) outer ring's own
    gap, so a straight continuation in that same direction hits the outer
    ring -- exactly the "two-stage exit" failure mode `creep_out` (a real
    bent path) exists to solve.
    """
    gap = _away_bearing_deg(center, point)
    inner = _ring_disks(point, 60.0, gap, 14.0, id_base)
    outer = _ring_disks(point, 140.0, gap + 90.0, 14.0, id_base + 1000)
    return inner + outer


class TestEscapeToRingFailsWhenBoxedIn:
    def test_no_clean_single_leg_exit_exists(self) -> None:
        center = (0.0, 0.0)
        point = (0.0, 0.0)
        base_primitives = _boxed_in_rings(center, point, id_base=0)
        exit_path = floor.escape_to_ring(
            0, point, center, 400.0, PARAMS, base_primitives, [], PAIR_MAP
        )
        assert exit_path is None


class TestCreepOutEscapesTheBoxedInPoint:
    def test_creep_out_produces_a_clean_polyline(self) -> None:
        center = (0.0, 0.0)
        point = (0.0, 0.0)
        base_primitives = _boxed_in_rings(center, point, id_base=0)
        waypoints = floor.creep_out(0, point, center, 400.0, PARAMS, base_primitives, [], PAIR_MAP)
        assert waypoints is not None
        assert waypoints[0] == point
        segments = polyline_capsules(waypoints, 0, 0, PARAMS.pair_half_width, line_uid=1)
        assert polyline_is_clean(segments, base_primitives, [], PAIR_MAP, PARAMS.tol)

    def test_creep_out_reaches_at_least_the_target_radius(self) -> None:
        center = (0.0, 0.0)
        point = (0.0, 0.0)
        base_primitives = _boxed_in_rings(center, point, id_base=0)
        waypoints = floor.creep_out(0, point, center, 400.0, PARAMS, base_primitives, [], PAIR_MAP)
        assert waypoints is not None
        last = waypoints[-1]
        dist = (last[0] - center[0]) ** 2 + (last[1] - center[1]) ** 2
        assert dist >= (400.0 - 1e-6) ** 2

    def test_creep_out_returns_none_when_fully_sealed(self) -> None:
        center = (0.0, 0.0)
        point = (0.0, 0.0)
        sealed = _sealed_ring(point, radius=50.0, disk_radius=12.0, n=20, id_base=0)
        assert floor.creep_out(0, point, center, 200.0, PARAMS, sealed, [], PAIR_MAP) is None


class TestFindFloorRouteBetweenTwoBoxedInEndpoints:
    def test_whole_assembled_route_validates_clean(self) -> None:
        # Two independently boxed-in endpoints, far enough apart that
        # neither's local ring obstacles interfere with the other's exit
        # or with the connecting arc.
        p_i = (0.0, 0.0)
        p_j = (2000.0, 0.0)
        center = (1000.0, 0.0)
        base_primitives = [
            *_boxed_in_rings(center, p_i, id_base=0),
            *_boxed_in_rings(center, p_j, id_base=2000),
        ]
        line = floor.find_floor_route(
            0, 1, p_i, p_j, center, 3000.0, PARAMS, base_primitives, [], PAIR_MAP
        )
        assert line is not None
        assert line.points[0] == p_i
        assert line.points[-1] == p_j
        segments = polyline_capsules(line.points, 0, 1, PARAMS.pair_half_width, line_uid=1)
        assert polyline_is_clean(segments, base_primitives, [], PAIR_MAP, PARAMS.tol)

    def test_returns_none_when_an_endpoint_is_fully_sealed(self) -> None:
        # A gapless ring (adjacent disks overlap angularly, no seam at
        # all) has NO valid escape at any bearing or length -- a
        # legitimate `None`, honestly reported rather than fabricated.
        p_i, p_j = (0.0, 0.0), (500.0, 0.0)
        center = (0.0, 0.0)
        sealed = _sealed_ring(p_i, radius=50.0, disk_radius=12.0, n=20, id_base=0)
        line = floor.find_floor_route(0, 1, p_i, p_j, center, 200.0, PARAMS, sealed, [], PAIR_MAP)
        assert line is None


class TestEscapeToRingSucceedsWithNoObstacles:
    def test_short_direct_escape(self) -> None:
        center = (0.0, 0.0)
        point = (0.0, 0.0)
        exit_path = floor.escape_to_ring(0, point, center, 100.0, PARAMS, [], [], PAIR_MAP)
        assert exit_path is not None
        assert exit_path[0] == point
        dist = (exit_path[-1][0] - center[0]) ** 2 + (exit_path[-1][1] - center[1]) ** 2
        assert abs(dist - 100.0**2) < 1e-6
