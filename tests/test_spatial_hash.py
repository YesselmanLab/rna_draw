"""Tests for `rna_draw.spatial_hash`: bucketing and the superset property."""

from __future__ import annotations

import itertools
import random

from rna_draw.spatial_hash import SpatialHash

Aabb = tuple[float, float, float, float]


def _aabb_overlap(a: Aabb, b: Aabb) -> bool:
    """Reference AABB-overlap test used only to validate the superset claim."""
    ax0, ay0, ax1, ay1 = a
    bx0, by0, bx1, by1 = b
    return ax0 <= bx1 and bx0 <= ax1 and ay0 <= by1 and by0 <= ay1


class TestSpatialHashBasics:
    def test_multicell_primitive_lands_in_every_touched_cell(self) -> None:
        grid = SpatialHash(cell_size=10.0)
        grid.insert(0, (0.0, 0.0, 25.0, 5.0))  # spans cells (0,0), (1,0), (2,0)
        assert grid._cells[(0, 0)] == [0]
        assert grid._cells[(1, 0)] == [0]
        assert grid._cells[(2, 0)] == [0]

    def test_overlapping_aabbs_share_a_bucket(self) -> None:
        grid = SpatialHash(cell_size=10.0)
        grid.insert(0, (0.0, 0.0, 5.0, 5.0))
        grid.insert(1, (2.0, 2.0, 7.0, 7.0))
        assert (0, 1) in grid.candidate_pairs()

    def test_far_apart_aabbs_do_not_share_a_bucket(self) -> None:
        grid = SpatialHash(cell_size=10.0)
        grid.insert(0, (0.0, 0.0, 5.0, 5.0))
        grid.insert(1, (1000.0, 1000.0, 1005.0, 1005.0))
        assert grid.candidate_pairs() == set()

    def test_candidate_pairs_deduplicated_across_shared_cells(self) -> None:
        grid = SpatialHash(cell_size=5.0)
        # Both items span multiple cells and overlap in more than one bucket.
        grid.insert(0, (0.0, 0.0, 12.0, 12.0))
        grid.insert(1, (1.0, 1.0, 11.0, 11.0))
        pairs = grid.candidate_pairs()
        assert pairs == {(0, 1)}


class TestInsertSegment:
    """`insert_segment` tiles a long segment into small pieces instead of
    one huge AABB, but must stay a superset of the true footprint's cells.
    """

    def test_zero_length_segment_falls_back_to_a_point_aabb(self) -> None:
        grid = SpatialHash(cell_size=10.0)
        grid.insert_segment(0, 5.0, 5.0, 5.0, 5.0, pad=2.0)
        assert grid._cells == {(0, 0): [0]}

    def test_short_segment_matches_a_single_padded_aabb_insert(self) -> None:
        # A segment shorter than one cell should need exactly one piece,
        # landing in the same cells as a plain padded-AABB `insert`.
        via_segment = SpatialHash(cell_size=10.0)
        via_segment.insert_segment(0, 0.0, 0.0, 3.0, 4.0, pad=1.0)

        via_aabb = SpatialHash(cell_size=10.0)
        via_aabb.insert(0, (0.0 - 1.0, 0.0 - 1.0, 3.0 + 1.0, 4.0 + 1.0))

        assert via_segment._cells == via_aabb._cells

    def test_long_diagonal_touches_far_fewer_cells_than_one_big_aabb(self) -> None:
        # A near-diagonal segment: one whole-AABB insert touches ~(L/cell)^2
        # cells, but tiling touches ~(L/cell) -- the fix for the perf bug.
        cell_size = 10.0
        x0, y0, x1, y1 = 0.0, 0.0, 990.0, 995.0

        tiled = SpatialHash(cell_size)
        tiled.insert_segment(0, x0, y0, x1, y1, pad=2.0)

        whole = SpatialHash(cell_size)
        whole.insert(
            0, (min(x0, x1) - 2.0, min(y0, y1) - 2.0, max(x0, x1) + 2.0, max(y0, y1) + 2.0)
        )

        assert len(tiled._cells) < len(whole._cells) / 10

    def test_tiled_cells_are_a_subset_of_the_whole_padded_aabb(self) -> None:
        # Tiling must never claim cells the (correctness-proven) whole-AABB
        # insert wouldn't also cover.
        cell_size = 10.0
        x0, y0, x1, y1 = 3.0, -7.0, 123.0, 58.0
        pad = 4.0

        tiled = SpatialHash(cell_size)
        tiled.insert_segment(0, x0, y0, x1, y1, pad)

        whole = SpatialHash(cell_size)
        whole_aabb = (min(x0, x1) - pad, min(y0, y1) - pad, max(x0, x1) + pad, max(y0, y1) + pad)
        whole.insert(0, whole_aabb)

        assert set(tiled._cells) <= set(whole._cells)

    def test_candidate_pairs_still_superset_of_true_overlaps_for_segments(self) -> None:
        # A point near the middle of a long diagonal segment, well inside
        # its `pad` footprint but far from either endpoint -- exactly the
        # case a whole-endpoint-only AABB could miss if padding were wrong.
        cell_size = 10.0
        grid = SpatialHash(cell_size)
        grid.insert_segment(0, 0.0, 0.0, 500.0, 500.0, pad=3.0)
        # Point at t=0.5 on the segment, offset by < pad perpendicular.
        grid.insert(1, (250.0 - 1.0, 250.0 - 1.0, 250.0 + 1.0, 250.0 + 1.0))
        assert (0, 1) in grid.candidate_pairs()


def _random_aabbs(n: int, box: float, size: float, seed: int) -> list[Aabb]:
    rng = random.Random(seed)
    aabbs = []
    for _ in range(n):
        cx, cy = rng.uniform(0, box), rng.uniform(0, box)
        half = rng.uniform(size * 0.5, size * 1.5)
        aabbs.append((cx - half, cy - half, cx + half, cy + half))
    return aabbs


class TestSupersetProperty:
    def test_candidate_pairs_is_superset_of_true_overlaps(self) -> None:
        aabbs = _random_aabbs(n=150, box=200.0, size=10.0, seed=7)
        true_overlaps = {
            (i, j)
            for i, j in itertools.combinations(range(len(aabbs)), 2)
            if _aabb_overlap(aabbs[i], aabbs[j])
        }

        for cell_size in (1.0, 5.0, 10.0, 25.0, 100.0):
            grid = SpatialHash(cell_size=cell_size)
            for index, aabb in enumerate(aabbs):
                grid.insert(index, aabb)
            candidates = grid.candidate_pairs()
            assert true_overlaps <= candidates

    def test_superset_property_many_seeds(self) -> None:
        for seed in range(10):
            aabbs = _random_aabbs(n=60, box=100.0, size=8.0, seed=seed)
            true_overlaps = {
                (i, j)
                for i, j in itertools.combinations(range(len(aabbs)), 2)
                if _aabb_overlap(aabbs[i], aabbs[j])
            }
            grid = SpatialHash(cell_size=15.0)
            for index, aabb in enumerate(aabbs):
                grid.insert(index, aabb)
            assert true_overlaps <= grid.candidate_pairs()
