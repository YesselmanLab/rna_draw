"""Tests for `rna_draw.layout.pseudoknot.placement`/`validate`/`routing`:
the PK-B/PK-A escalation, the geometry-predicate polyline validator, and
the axis-aligned staple routing tiers, at the unit level (independent of
the full `layout_pseudoknot` orchestration -- see `test_pseudoknot_engine.py`).
"""

from __future__ import annotations

from rna_draw.geometry import Capsule, Disk, PrimitiveId
from rna_draw.layout.base import RoutedLine
from rna_draw.layout.pseudoknot import routing
from rna_draw.layout.pseudoknot.parsing import Stem
from rna_draw.layout.pseudoknot.placement import (
    MAX_PK_B_LENGTH,
    _PlacementState,
    _route_pair,
    _try_in_plane,
    place_crossings,
)
from rna_draw.layout.pseudoknot.routing import enclosing_circle
from rna_draw.layout.pseudoknot.validate import polyline_capsules, polyline_is_clean
from rna_draw.overlap import OverlapParams, build_primitives


class TestPlaceCrossingsPKB:
    def test_short_isolated_crossing_becomes_in_plane_connector(self) -> None:
        # Two nucleotides far apart from everything else: a short (length
        # 1) crossing pair between them should stay checker-clean as a
        # straight PK-B connector -- no reason to escalate to PK-A.
        x = [0.0, 1000.0]
        y = [0.0, 0.0]
        pair_map = [-1, -1]
        stem = Stem(i=0, j=1, length=1)
        result = place_crossings(x, y, pair_map, [stem], OverlapParams())
        assert result.pk_b_pairs == [(0, 1)]
        assert result.pk_a_lines == []
        assert result.unplaced == []

    def test_no_crossings_is_a_no_op(self) -> None:
        result = place_crossings([0.0, 10.0], [0.0, 0.0], [-1, -1], [], OverlapParams())
        assert result.pk_b_pairs == []
        assert result.pk_a_lines == []
        assert result.unplaced == []

    def test_stem_longer_than_max_pk_b_length_never_becomes_in_plane(self) -> None:
        # A long crossing stem laid out along two long, well-separated
        # straight tracks: the straight PK-B rungs WOULD be checker-clean,
        # but MUST-FIX #4 forbids keeping a stem this long in-plane.
        n = MAX_PK_B_LENGTH + 1
        x = [float(k) * 30 for k in range(n)] + [float(k) * 30 for k in range(n)]
        y = [0.0] * n + [1000.0] * n
        pair_map = [-1] * (2 * n)
        stem = Stem(i=0, j=2 * n - 1, length=n)
        result = place_crossings(x, y, pair_map, [stem], OverlapParams())
        assert result.pk_b_pairs == []


class TestPlaceCrossingsPKA:
    # Three co-linear nucleotides at ordinary backbone spacing (30 units,
    # comfortably more than `2 * node_r == 20`, so only the crossing pair
    # itself -- not the backbone -- is the thing under test): pairing the
    # two ends would draw a straight capsule directly through the middle
    # disk's own footprint, so PK-B must reject it, forcing PK-A.
    X = [0.0, 30.0, 60.0]
    Y = [0.0, 0.0, 0.0]
    PAIR_MAP = [-1, -1, -1]

    def test_blocked_straight_chord_escalates_to_routed_lines(self) -> None:
        stem = Stem(i=0, j=2, length=1)
        result = place_crossings(self.X, self.Y, self.PAIR_MAP, [stem], OverlapParams())
        assert result.pk_b_pairs == []
        assert len(result.pk_a_lines) == 1
        assert result.unplaced == []

    def test_routed_line_touches_only_its_own_endpoints(self) -> None:
        stem = Stem(i=0, j=2, length=1)
        result = place_crossings(self.X, self.Y, self.PAIR_MAP, [stem], OverlapParams())
        line = result.pk_a_lines[0]
        assert line.points[0] == (self.X[0], self.Y[0])
        assert line.points[-1] == (self.X[2], self.Y[2])


class TestDensePseudoknotAllPairsPlaced:
    """The staple quality ladder (direct/staple/offset-staple) placing
    EVERY crossing pair of a moderately dense pseudoknot -- the plan's
    core deliverable (UNPLACED -> ~0 for the reachable common case; see
    `benchmarks/pseudoknot_gate.py` for the honest corpus-scale numbers,
    which do have a residual -- documented -- for genuinely embedded
    real-structure endpoints).
    """

    def test_zero_unplaced(self) -> None:
        # 10 co-linear nucleotides (ordinary 30-unit backbone spacing)
        # with 4 nested crossing stems, each pair's straight chord
        # blocked by the nucleotides between it -- every one must
        # escalate to (and succeed at) a PK-A routed staple.
        n = 10
        x = [float(k) * 30 for k in range(n)]
        y = [0.0] * n
        pair_map = [-1] * n
        stems = [
            Stem(i=0, j=9, length=1),
            Stem(i=1, j=8, length=1),
            Stem(i=2, j=7, length=1),
            Stem(i=3, j=6, length=1),
        ]
        result = place_crossings(x, y, pair_map, stems, OverlapParams())
        assert result.unplaced == []
        assert len(result.pk_a_lines) == 4

    def test_every_routed_line_is_mutually_clean(self) -> None:
        n = 10
        x = [float(k) * 30 for k in range(n)]
        y = [0.0] * n
        pair_map = [-1] * n
        stems = [
            Stem(i=0, j=9, length=1),
            Stem(i=1, j=8, length=1),
            Stem(i=2, j=7, length=1),
            Stem(i=3, j=6, length=1),
        ]
        params = OverlapParams()
        result = place_crossings(x, y, pair_map, stems, params)
        base_primitives = build_primitives(x, y, pair_map, params)
        all_caps = [
            polyline_capsules(line.points, line.i, line.j, params.pair_half_width, uid)
            for uid, line in enumerate(result.pk_a_lines, start=1)
        ]
        for uid, segments in enumerate(all_caps):
            others = [s for ouid, caps in enumerate(all_caps) if ouid != uid for s in caps]
            assert polyline_is_clean(segments, base_primitives, others, pair_map, params.tol)


class TestTwoCrossingLinesRouteDisjoint:
    """Two crossing pairs whose INDEPENDENTLY-computed best routes would
    physically collide (same shared rail) must end up mutually clean once
    threaded through the SAME `committed_lines` -- the second escalates to
    a genuinely different, disjoint route rather than either colliding or
    going unplaced.
    """

    N = 10
    X = [float(k) * 30 for k in range(N)]
    Y = [0.0] * N
    PAIR_MAP = [-1] * N

    def _route(self, i: int, j: int, committed: list[Capsule]) -> RoutedLine | None:
        params = OverlapParams()
        base_primitives = build_primitives(self.X, self.Y, self.PAIR_MAP, params)
        center, _base_radius = enclosing_circle(self.X, self.Y, params)
        state = _PlacementState(
            x=self.X,
            y=self.Y,
            params=params,
            pair_map=self.PAIR_MAP,
            committed_lines=committed,
            center=center,
            extent=max(max(self.X) - min(self.X), params.node_r),
        )
        return _route_pair(state, i, j, base_primitives)

    def test_naive_independent_routes_would_collide(self) -> None:
        # Sanity: WITHOUT committed-line knowledge, pair (1, 8)'s own best
        # staple lands on the exact same shared rail (y) as pair (0, 9)'s
        # (both co-linear chords share a midpoint x, so both independently
        # pick the same cardinal direction and offset), and pair (1, 8)'s
        # rail x-span sits entirely WITHIN pair (0, 9)'s -- their "across"
        # segments would overlap outright -- proving this scenario
        # actually exercises disjoint-routing, not a vacuous pass.
        line_a = self._route(0, 9, [])
        line_b_naive = self._route(1, 8, [])
        assert line_a is not None
        assert line_b_naive is not None
        rail_a_y = line_a.points[1][1]
        rail_b_y = line_b_naive.points[1][1]
        assert rail_a_y == rail_b_y
        rail_a_x = sorted(p[0] for p in line_a.points[1:3])
        rail_b_x = sorted(p[0] for p in line_b_naive.points[1:3])
        assert rail_a_x[0] <= rail_b_x[0] and rail_b_x[1] <= rail_a_x[1]

    def test_informed_second_route_is_disjoint_and_clean(self) -> None:
        params = OverlapParams()
        base_primitives = build_primitives(self.X, self.Y, self.PAIR_MAP, params)
        line_a = self._route(0, 9, [])
        assert line_a is not None
        caps_a = polyline_capsules(line_a.points, 0, 9, params.pair_half_width, line_uid=1)

        line_b = self._route(1, 8, caps_a)
        assert line_b is not None
        caps_b = polyline_capsules(line_b.points, 1, 8, params.pair_half_width, line_uid=2)

        assert line_b.points != line_a.points
        assert polyline_is_clean(caps_b, base_primitives, caps_a, self.PAIR_MAP, params.tol)
        assert polyline_is_clean(caps_a, base_primitives, caps_b, self.PAIR_MAP, params.tol)


class TestPkBRejectedAgainstCommittedPkALine:
    """Regression for BUG 1: `_try_in_plane` used to validate a candidate
    PK-B rung only against `check_overlaps(candidate_pair_map)`, which is
    blind to already-committed PK-A routed lines (they never enter any
    `pair_map`). A rung that is clean against every OTHER primitive but
    slices straight through an earlier crossing stem's committed line must
    now be rejected (demoted to PK-A) instead of silently accepted.
    """

    def test_rung_clean_against_pair_map_but_crossing_committed_line_is_rejected(self) -> None:
        params = OverlapParams()
        x = [0.0, 100.0, 0.0, 100.0]
        y = [0.0, 0.0, 100.0, 100.0]
        pair_map = [-1, -1, -1, -1]
        # A committed PK-A line from an unrelated, already-placed crossing
        # stem: a vertical segment at x=50 that the candidate rung's
        # straight (0,0)-(100,0) chord must cross at (50, 0).
        committed_lines = polyline_capsules(
            [(50.0, -50.0), (50.0, 150.0)], 9, 10, params.pair_half_width, line_uid=1
        )
        state = _PlacementState(
            x=x, y=y, params=params, pair_map=pair_map, committed_lines=committed_lines
        )
        stem = Stem(i=0, j=1, length=1)
        assert _try_in_plane(state, stem) is False

    def test_rung_accepted_when_no_committed_line_is_in_the_way(self) -> None:
        params = OverlapParams()
        x = [0.0, 100.0, 0.0, 100.0]
        y = [0.0, 0.0, 100.0, 100.0]
        pair_map = [-1, -1, -1, -1]
        state = _PlacementState(x=x, y=y, params=params, pair_map=pair_map, committed_lines=[])
        stem = Stem(i=0, j=1, length=1)
        assert _try_in_plane(state, stem) is True


class TestEnclosingCircle:
    def test_encloses_every_point(self) -> None:
        x = [0.0, 50.0, -30.0, 10.0]
        y = [0.0, 20.0, -10.0, 40.0]
        center, radius = enclosing_circle(x, y, OverlapParams())
        for px, py in zip(x, y):
            dist = ((px - center[0]) ** 2 + (py - center[1]) ** 2) ** 0.5
            assert dist <= radius

    def test_single_point_layout_has_positive_radius(self) -> None:
        _center, radius = enclosing_circle([5.0], [5.0], OverlapParams())
        assert radius > 0.0


def _assert_axis_aligned(points: list[tuple[float, float]]) -> None:
    """Every consecutive segment is purely horizontal or purely vertical."""
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        assert x0 == x1 or y0 == y1, f"diagonal segment {(x0, y0)} -> {(x1, y1)}"


class TestStapleRoute:
    def test_succeeds_with_no_obstacles(self) -> None:
        # `routing.find_staple_route` (tier 1) exercised directly: with
        # empty `base_primitives`/`committed_lines`, the very first
        # (direction, offset) candidate must succeed immediately.
        params = OverlapParams()
        center, _radius = enclosing_circle([0.0, 20.0], [0.0, 0.0], params)
        line = routing.find_staple_route(
            0, 1, (0.0, 0.0), (20.0, 0.0), center, 20.0, params, [], [], [-1, -1]
        )
        assert line is not None
        assert line.i == 0
        assert line.j == 1
        assert line.points[0] == (0.0, 0.0)
        assert line.points[-1] == (20.0, 0.0)

    def test_is_a_3_segment_orthogonal_staple(self) -> None:
        params = OverlapParams()
        center, _radius = enclosing_circle([0.0, 20.0], [0.0, 0.0], params)
        line = routing.find_staple_route(
            0, 1, (0.0, 0.0), (20.0, 0.0), center, 20.0, params, [], [], [-1, -1]
        )
        assert line is not None
        assert len(line.points) == 4
        _assert_axis_aligned(line.points)

    def test_never_produces_a_diagonal_segment(self) -> None:
        # A pair whose chord is itself diagonal (different x AND y): the
        # staple must still route using only horizontal/vertical segments.
        params = OverlapParams()
        center, _radius = enclosing_circle([0.0, 40.0], [0.0, 30.0], params)
        line = routing.find_staple_route(
            0, 1, (0.0, 0.0), (40.0, 30.0), center, 50.0, params, [], [], [-1, -1]
        )
        assert line is not None
        _assert_axis_aligned(line.points)


class TestOffsetStapleRoute:
    def test_succeeds_with_no_obstacles(self) -> None:
        params = OverlapParams()
        center, _radius = enclosing_circle([0.0, 20.0], [0.0, 0.0], params)
        line = routing.find_offset_staple_route(
            0, 1, (0.0, 0.0), (20.0, 0.0), center, 20.0, params, [], [], [-1, -1]
        )
        assert line is not None
        _assert_axis_aligned(line.points)

    def test_never_produces_a_diagonal_segment(self) -> None:
        params = OverlapParams()
        center, _radius = enclosing_circle([0.0, 40.0], [0.0, 30.0], params)
        line = routing.find_offset_staple_route(
            0, 1, (0.0, 0.0), (40.0, 30.0), center, 50.0, params, [], [], [-1, -1]
        )
        assert line is not None
        _assert_axis_aligned(line.points)


class TestPolylineValidate:
    def test_clean_polyline_far_from_everything(self) -> None:
        points = [(0.0, 0.0), (500.0, 500.0), (1000.0, 0.0)]
        segments = polyline_capsules(points, 0, 1, half_width=7.5, line_uid=0)
        base_primitives = [
            Disk(pid=PrimitiveId("nt", 0), cx=0.0, cy=0.0, radius=10.0, ends=frozenset({0})),
            Disk(pid=PrimitiveId("nt", 1), cx=1000.0, cy=0.0, radius=10.0, ends=frozenset({1})),
        ]
        assert polyline_is_clean(segments, base_primitives, [], [-1, -1], 1e-6) is True

    def test_polyline_through_an_unrelated_disk_is_dirty(self) -> None:
        points = [(0.0, 0.0), (50.0, 0.0), (100.0, 0.0)]
        segments = polyline_capsules(points, 0, 2, half_width=7.5, line_uid=0)
        obstacle = Disk(pid=PrimitiveId("nt", 9), cx=50.0, cy=0.0, radius=10.0, ends=frozenset({9}))
        pair_map = [-1] * 3 + [-1] * 7
        assert polyline_is_clean(segments, [obstacle], [], pair_map, 1e-6) is False

    def test_committed_line_blocks_a_new_one(self) -> None:
        points = [(0.0, 0.0), (50.0, 0.0), (100.0, 0.0)]
        segments = polyline_capsules(points, 0, 2, half_width=7.5, line_uid=1)
        committed = polyline_capsules(
            [(50.0, -50.0), (50.0, 50.0)], 5, 6, half_width=7.5, line_uid=2
        )
        pair_map = [-1] * 10
        assert polyline_is_clean(segments, [], committed, pair_map, 1e-6) is False

    def test_build_primitives_smoke(self) -> None:
        # Sanity: `build_primitives` (frozen `overlap.py`) is usable
        # directly as `base_primitives` input, per the module contract.
        x, y = [0.0, 20.0], [0.0, 0.0]
        pair_map = [1, 0]
        primitives = build_primitives(x, y, pair_map, OverlapParams())
        assert len(primitives) > 0


def test_capsule_frozen_dataclass_hashable_ends() -> None:
    # `ends` must be a frozenset (not a plain set) for Capsule to remain
    # usable in the exclusion-check codepaths that treat it as such.
    cap = Capsule(
        pid=PrimitiveId("pkline", 0),
        x0=0.0,
        y0=0.0,
        x1=1.0,
        y1=1.0,
        half_width=1.0,
        ends=frozenset({0, 1}),
    )
    assert isinstance(cap.ends, frozenset)
