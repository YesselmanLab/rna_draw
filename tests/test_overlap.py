"""Tests for `rna_draw.overlap`: building, exclusions, checking, brute-force
equivalence, and performance.
"""

from __future__ import annotations

import random
import time

import pytest

from rna_draw.geometry import Capsule, Disk, PrimitiveId
from rna_draw.overlap import (
    OverlapKind,
    OverlapParams,
    OverlapReport,
    Witness,
    build_backbone_capsules,
    build_disks,
    build_pair_capsules,
    build_primitives,
    check_overlaps,
    check_overlaps_bruteforce,
    is_excluded,
    rescale_coords,
)
from rna_draw.overlap import test_pair as _test_pair
from rna_draw.render_rna import RNARenderer, get_pairmap_from_secstruct
from rna_draw.spatial_hash import SpatialHash

NODE_R = 10
PRIMARY_SPACE = 20
PAIR_SPACE = 23


def _hairpin_coords(secstruct: str) -> tuple[list[float], list[float], list[int]]:
    """Lay out `secstruct` with the legacy renderer and return coords + pair_map."""
    renderer = RNARenderer()
    renderer.setup_tree(
        secstruct, NODE_R=NODE_R, PRIMARY_SPACE=PRIMARY_SPACE, PAIR_SPACE=PAIR_SPACE
    )
    x, y = list(renderer.xarray_), list(renderer.yarray_)
    pair_map = get_pairmap_from_secstruct(secstruct)
    return x, y, pair_map


def _random_case(n: int, seed: int) -> tuple[list[float], list[float], list[int]]:
    """Build scattered coordinates + a random valid symmetric pair_map.

    Coordinates are packed into a small box (relative to the default
    `node_r`) so overlaps are common, exercising both the geometry
    predicates and the exclusion rules.

    Args:
        n: Number of nucleotides.
        seed: RNG seed for reproducibility.

    Returns:
        `(x, y, pair_map)`.
    """
    rng = random.Random(seed)
    box = 4.0 * NODE_R
    x = [rng.uniform(0, box) for _ in range(n)]
    y = [rng.uniform(0, box) for _ in range(n)]
    pair_map = [-1] * n
    shuffled = list(range(n))
    rng.shuffle(shuffled)
    num_pairs = rng.randint(0, n // 2)
    for k in range(num_pairs):
        i, j = shuffled[2 * k], shuffled[2 * k + 1]
        pair_map[i] = j
        pair_map[j] = i
    return x, y, pair_map


def _find_witness(report: OverlapReport, id_a: PrimitiveId, id_b: PrimitiveId) -> Witness | None:
    """Locate a witness for an unordered pair of primitive ids, if present."""
    wanted = {id_a, id_b}
    for witness in report.witnesses:
        if {witness.id_a, witness.id_b} == wanted:
            return witness
    return None


class TestBuildPrimitives:
    def test_disk_count_and_ends(self) -> None:
        x, y = [0.0, 10.0, 20.0], [0.0, 0.0, 0.0]
        disks = build_disks(x, y, node_r=5.0)
        assert len(disks) == 3
        assert [d.pid for d in disks] == [PrimitiveId("nt", i) for i in range(3)]
        assert [d.ends for d in disks] == [frozenset({i}) for i in range(3)]

    def test_backbone_capsule_count_and_ends(self) -> None:
        x, y = [0.0, 10.0, 20.0], [0.0, 0.0, 0.0]
        capsules = build_backbone_capsules(x, y, half_width=3.0)
        assert len(capsules) == 2
        assert [c.ends for c in capsules] == [frozenset({0, 1}), frozenset({1, 2})]

    def test_pair_capsule_count_and_ends(self) -> None:
        x, y = [0.0, 10.0, 20.0, 30.0], [0.0, 0.0, 0.0, 0.0]
        pair_map = [3, -1, -1, 0]
        capsules = build_pair_capsules(x, y, pair_map, half_width=3.0)
        assert len(capsules) == 1
        assert capsules[0].pid == PrimitiveId("pair", 0)
        assert capsules[0].ends == frozenset({0, 3})

    def test_build_primitives_total_and_order(self) -> None:
        x, y = [0.0, 10.0, 20.0, 30.0], [0.0, 0.0, 0.0, 0.0]
        pair_map = [3, -1, -1, 0]
        params = OverlapParams()
        primitives = build_primitives(x, y, pair_map, params)
        assert len(primitives) == 4 + 3 + 1
        assert all(isinstance(p, Disk) for p in primitives[:4])
        assert all(isinstance(p, Capsule) for p in primitives[4:7])
        assert all(isinstance(p, Capsule) for p in primitives[7:])


class TestIsExcluded:
    def test_adjacent_disks_excluded(self) -> None:
        d0 = Disk(PrimitiveId("nt", 0), 0, 0, 10, frozenset({0}))
        d1 = Disk(PrimitiveId("nt", 1), 5, 0, 10, frozenset({1}))
        assert is_excluded(d0, d1, pair_map=[-1, -1]) is True

    def test_paired_disks_excluded(self) -> None:
        d0 = Disk(PrimitiveId("nt", 0), 0, 0, 10, frozenset({0}))
        d5 = Disk(PrimitiveId("nt", 5), 100, 0, 10, frozenset({5}))
        pair_map = [5, -1, -1, -1, -1, 0]
        assert is_excluded(d0, d5, pair_map) is True

    def test_disk_excluded_from_its_own_capsule(self) -> None:
        disk = Disk(PrimitiveId("nt", 0), 0, 0, 10, frozenset({0}))
        capsule = Capsule(PrimitiveId("bb", 0), 0, 0, 20, 0, 5, frozenset({0, 1}))
        assert is_excluded(disk, capsule, pair_map=[-1, -1]) is True

    def test_capsules_sharing_endpoint_excluded(self) -> None:
        c1 = Capsule(PrimitiveId("bb", 0), 0, 0, 20, 0, 5, frozenset({0, 1}))
        c2 = Capsule(PrimitiveId("bb", 1), 20, 0, 40, 0, 5, frozenset({1, 2}))
        assert is_excluded(c1, c2, pair_map=[-1, -1, -1]) is True

    def test_distant_disks_not_excluded(self) -> None:
        d0 = Disk(PrimitiveId("nt", 0), 0, 0, 10, frozenset({0}))
        d5 = Disk(PrimitiveId("nt", 5), 500, 0, 10, frozenset({5}))
        pair_map = [-1, -1, -1, -1, -1, -1]
        assert is_excluded(d0, d5, pair_map) is False


class TestPairOrdering:
    """`test_pair` must dispatch the same regardless of argument order."""

    def test_capsule_then_disk_matches_disk_then_capsule(self) -> None:
        disk = Disk(PrimitiveId("nt", 0), 5, 0, 10, frozenset({0}))
        capsule = Capsule(PrimitiveId("bb", 1), 0, 0, 10, 0, 5, frozenset({1, 2}))
        forward = _test_pair(disk, capsule, tol=1e-6)
        reversed_order = _test_pair(capsule, disk, tol=1e-6)
        assert forward is not None
        assert reversed_order is not None
        assert forward.kind == reversed_order.kind == OverlapKind.DISK_CAPSULE
        assert forward.overlap_depth == pytest.approx(reversed_order.overlap_depth)


class TestKnownGood:
    def test_clean_hairpin_passes(self) -> None:
        x, y, pair_map = _hairpin_coords("((((....))))")
        report = check_overlaps(x, y, pair_map)
        assert report.passed is True
        assert report.witnesses == []


class TestKnownBad:
    def test_identical_coords_far_in_index_flagged(self) -> None:
        # nucleotides 0 and 3 are neither backbone-adjacent nor paired, yet
        # sit at the exact same point.
        x, y = [0.0, 30.0, 60.0, 0.0], [0.0, 0.0, 0.0, 0.0]
        pair_map = [-1, -1, -1, -1]
        report = check_overlaps(x, y, pair_map)
        witness = _find_witness(report, PrimitiveId("nt", 0), PrimitiveId("nt", 3))
        assert witness is not None
        assert witness.kind == OverlapKind.DISK_DISK
        assert witness.id_a < witness.id_b

    def test_disk_on_foreign_pair_capsule_flagged(self) -> None:
        # Nucleotides 0 and 3 are paired; nucleotide 1 is unrelated but sits
        # exactly on the (0, 3) pair connector's axis.
        x, y = [0.0, 50.0, 75.0, 100.0], [0.0, 0.0, 50.0, 0.0]
        pair_map = [3, -1, -1, 0]
        report = check_overlaps(x, y, pair_map)
        witness = _find_witness(report, PrimitiveId("nt", 1), PrimitiveId("pair", 0))
        assert witness is not None
        assert witness.kind == OverlapKind.DISK_CAPSULE
        assert witness.id_a < witness.id_b

    def test_crossing_pair_capsules_flagged(self) -> None:
        # Corners of a square; pairs (0, 1) and (2, 3) are its two diagonals,
        # which cross at the center. Disks/backbone shrunk to isolate the
        # pair-capsule crossing. Pair-capsule pids are keyed by the lower
        # nucleotide index, so the two connectors are ("pair", 0) and
        # ("pair", 2).
        x, y = [0.0, 10.0, 10.0, 0.0], [0.0, 10.0, 0.0, 10.0]
        pair_map = [1, 0, 3, 2]
        params = OverlapParams(node_r=0.01, backbone_half_width=0.01, pair_half_width=1.0)
        report = check_overlaps(x, y, pair_map, params)
        witness = _find_witness(report, PrimitiveId("pair", 0), PrimitiveId("pair", 2))
        assert witness is not None
        assert witness.kind == OverlapKind.CAPSULE_CAPSULE
        assert witness.id_a < witness.id_b


class TestLoopLocalCloseness:
    """Risk R2: geometry, not index distance, decides overlap."""

    # A zig-zag of 5 nucleotides where every gap-2 pair (0,2), (1,3), (2,4)
    # sits `SEPARATION` apart -- closer than touching (`2 * node_r`) but not
    # backbone-adjacent or paired. Backbone/pair half-widths are zeroed so
    # the only signal under test is the disk-disk relationship.
    SEPARATION = 14.0

    def _coords(self) -> tuple[list[float], list[float], list[int]]:
        sep = self.SEPARATION
        x = [0.0, sep / 2, sep, sep + sep / 2, 2 * sep]
        y = [0.0, 200.0, 0.0, 200.0, 0.0]
        pair_map = [-1, -1, -1, -1, -1]
        return x, y, pair_map

    def test_flagged_at_default_node_r(self) -> None:
        x, y, pair_map = self._coords()
        params = OverlapParams(node_r=10.0, backbone_half_width=0.0, pair_half_width=0.0)
        report = check_overlaps(x, y, pair_map, params)
        assert report.passed is False
        for i, j in ((0, 2), (1, 3), (2, 4)):
            witness = _find_witness(report, PrimitiveId("nt", i), PrimitiveId("nt", j))
            assert witness is not None
            assert witness.kind == OverlapKind.DISK_DISK

    def test_clears_once_node_r_lowered_below_half_separation(self) -> None:
        x, y, pair_map = self._coords()
        small_node_r = self.SEPARATION / 2 - 1.0  # 2 * node_r < separation
        params = OverlapParams(node_r=small_node_r, backbone_half_width=0.0, pair_half_width=0.0)
        report = check_overlaps(x, y, pair_map, params)
        assert report.passed is True


class TestToleranceBoundary:
    def test_exactly_touching_disks_not_flagged(self) -> None:
        x, y = [0.0, 100.0, 20.0], [0.0, 100.0, 0.0]
        pair_map = [-1, -1, -1]
        params = OverlapParams(node_r=10.0, backbone_half_width=0.0, pair_half_width=0.0)
        report = check_overlaps(x, y, pair_map, params)
        assert report.passed is True


class TestHashEqualsBruteforce:
    @pytest.mark.parametrize("n", [5, 20, 50, 120])
    @pytest.mark.parametrize("seed", range(5))
    def test_witness_sets_match(self, n: int, seed: int) -> None:
        x, y, pair_map = _random_case(n, seed)
        hashed = check_overlaps(x, y, pair_map)
        brute = check_overlaps_bruteforce(x, y, pair_map)
        assert set(hashed.witnesses) == set(brute.witnesses)


class TestReportApi:
    def test_counts_by_kind_sum_to_num_overlaps(self) -> None:
        x, y, pair_map = _random_case(n=40, seed=1)
        report = check_overlaps(x, y, pair_map)
        assert sum(report.counts_by_kind.values()) == report.num_overlaps

    def test_passed_iff_no_witnesses(self) -> None:
        x, y, pair_map = _hairpin_coords("((((....))))")
        report = check_overlaps(x, y, pair_map)
        assert report.passed == (report.num_overlaps == 0)

    def test_invalid_lengths_raise(self) -> None:
        with pytest.raises(ValueError):
            check_overlaps([0.0, 1.0], [0.0], [-1, -1])

    def test_empty_input_raises(self) -> None:
        with pytest.raises(ValueError):
            check_overlaps([], [], [])

    def test_asymmetric_pair_map_raises(self) -> None:
        with pytest.raises(ValueError):
            check_overlaps([0.0, 1.0], [0.0, 1.0], [1, -1])


class TestRescaleCoords:
    def test_rescales_to_target_median_step(self) -> None:
        x, y = [0.0, 5.0, 10.0, 15.0], [0.0, 0.0, 0.0, 0.0]
        scaled_x, scaled_y = rescale_coords(x, y, target_backbone_step=20.0)
        assert scaled_x == pytest.approx([0.0, 20.0, 40.0, 60.0])
        assert scaled_y == pytest.approx([0.0, 0.0, 0.0, 0.0])

    def test_requires_at_least_two_points(self) -> None:
        with pytest.raises(ValueError):
            rescale_coords([0.0], [0.0], target_backbone_step=20.0)

    def test_zero_median_step_raises(self) -> None:
        # Every consecutive nucleotide sits at the same point, so the
        # median backbone step is zero and no scale factor can be derived.
        with pytest.raises(ValueError):
            rescale_coords([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], target_backbone_step=20.0)


class TestPerformance:
    def test_five_thousand_nt_completes_quickly(self) -> None:
        motif = "(((....)))" * 500  # 5000 nt
        x, y, pair_map = _hairpin_coords(motif)
        start = time.perf_counter()
        check_overlaps(x, y, pair_map)
        elapsed = time.perf_counter() - start
        assert elapsed < 5.0

    def test_candidate_count_much_less_than_brute_force(self) -> None:
        motif = "(((....)))" * 80  # ~800 nt
        x, y, pair_map = _hairpin_coords(motif)
        params = OverlapParams()
        primitives = build_primitives(x, y, pair_map, params)
        max_half_width = max(params.backbone_half_width, params.pair_half_width)
        cell_size = 2 * (params.node_r + max_half_width)
        grid = SpatialHash(cell_size)
        for index, primitive in enumerate(primitives):
            grid.insert(index, primitive.aabb)

        n = len(primitives)
        brute_total = n * (n - 1) // 2
        assert len(grid.candidate_pairs()) < 0.1 * brute_total
