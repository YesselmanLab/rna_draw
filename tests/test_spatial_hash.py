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
