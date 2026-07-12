"""Tests for `rna_draw.layout.pseudoknot.placement`/`validate`/`routing`:
the PK-B/PK-A escalation, the geometry-predicate polyline validator, and
the two routing tiers, at the unit level (independent of the full
`layout_pseudoknot` orchestration -- see `test_pseudoknot_engine.py`).
"""

from __future__ import annotations

from rna_draw.geometry import Capsule, Disk, PrimitiveId
from rna_draw.layout.pseudoknot import routing
from rna_draw.layout.pseudoknot.parsing import Stem
from rna_draw.layout.pseudoknot.placement import MAX_PK_B_LENGTH, place_crossings
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


class TestRingRoute:
    def test_succeeds_with_no_obstacles(self) -> None:
        # `routing.find_ring_route` (tier 2) exercised directly: with an
        # empty `base_primitives`/`committed_lines`, the very first ring
        # scale must succeed immediately.
        params = OverlapParams()
        center, base_radius = enclosing_circle([0.0, 20.0], [0.0, 0.0], params)
        line = routing.find_ring_route(
            0, 1, (0.0, 0.0), (20.0, 0.0), center, base_radius, params, [], [], [-1, -1]
        )
        assert line is not None
        assert line.i == 0
        assert line.j == 1
        assert line.points[0] == (0.0, 0.0)
        assert line.points[-1] == (20.0, 0.0)


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
