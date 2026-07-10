"""Unit tests for `rna_draw.geometry` distance math and overlap predicates."""

from __future__ import annotations

import math

import pytest

from rna_draw.geometry import (
    Capsule,
    Disk,
    PrimitiveId,
    capsule_capsule_overlap,
    clamp,
    disk_capsule_overlap,
    disks_overlap,
    point_segment_distance,
    segment_segment_distance,
)


def _disk(index: int, cx: float, cy: float, radius: float) -> Disk:
    return Disk(pid=PrimitiveId("nt", index), cx=cx, cy=cy, radius=radius, ends=frozenset({index}))


def _capsule(index: int, x0: float, y0: float, x1: float, y1: float, half_width: float) -> Capsule:
    return Capsule(
        pid=PrimitiveId("bb", index),
        x0=x0,
        y0=y0,
        x1=x1,
        y1=y1,
        half_width=half_width,
        ends=frozenset({index, index + 1}),
    )


class TestClamp:
    def test_within_range_unchanged(self) -> None:
        assert clamp(0.5, 0.0, 1.0) == 0.5

    def test_below_lo_clamped(self) -> None:
        assert clamp(-1.0, 0.0, 1.0) == 0.0

    def test_above_hi_clamped(self) -> None:
        assert clamp(2.0, 0.0, 1.0) == 1.0


class TestPointSegmentDistance:
    def test_interior_projection(self) -> None:
        # Point (5, 3) projects onto (0,0)-(10,0) at (5, 0): distance 3.
        assert point_segment_distance(5, 3, 0, 0, 10, 0) == pytest.approx(3.0)

    def test_clamped_to_start_endpoint(self) -> None:
        # Point (-5, 4) projects before the segment start; clamp to (0, 0).
        assert point_segment_distance(-5, 4, 0, 0, 10, 0) == pytest.approx(math.hypot(5, 4))

    def test_clamped_to_end_endpoint(self) -> None:
        # Point (15, 4) projects past the segment end; clamp to (10, 0).
        assert point_segment_distance(15, 4, 0, 0, 10, 0) == pytest.approx(math.hypot(5, 4))

    def test_degenerate_segment_is_point_distance(self) -> None:
        assert point_segment_distance(3, 4, 0, 0, 0, 0) == pytest.approx(5.0)


class TestSegmentSegmentDistance:
    def test_crossing_segments_touch(self) -> None:
        dist = segment_segment_distance((0, 0), (10, 10), (0, 10), (10, 0))
        assert dist == pytest.approx(0.0, abs=1e-9)

    def test_parallel_segments(self) -> None:
        dist = segment_segment_distance((0, 0), (10, 0), (0, 5), (10, 5))
        assert dist == pytest.approx(5.0)

    def test_skew_segments(self) -> None:
        # Segment 1 along x-axis; segment 2 a short vertical segment offset
        # in x, floating above -- closest points are each segment's nearest
        # endpoint/projection.
        dist = segment_segment_distance((0, 0), (10, 0), (20, 3), (20, 8))
        assert dist == pytest.approx(math.hypot(10, 3))

    def test_shared_endpoint_touches(self) -> None:
        dist = segment_segment_distance((0, 0), (5, 5), (5, 5), (10, 0))
        assert dist == pytest.approx(0.0, abs=1e-9)

    def test_collinear_overlapping_segments_touch(self) -> None:
        dist = segment_segment_distance((0, 0), (10, 0), (5, 0), (15, 0))
        assert dist == pytest.approx(0.0, abs=1e-9)

    def test_first_segment_degenerate(self) -> None:
        # Segment 1 collapses to a point off the axis of segment 2.
        dist = segment_segment_distance((5, 4), (5, 4), (0, 0), (10, 0))
        assert dist == pytest.approx(4.0)

    def test_second_segment_degenerate(self) -> None:
        # Segment 2 collapses to a point off the axis of segment 1.
        dist = segment_segment_distance((0, 0), (10, 0), (5, 4), (5, 4))
        assert dist == pytest.approx(4.0)

    def test_both_segments_degenerate(self) -> None:
        dist = segment_segment_distance((0, 0), (0, 0), (3, 4), (3, 4))
        assert dist == pytest.approx(5.0)


class TestDisksOverlap:
    def test_clearly_overlapping(self) -> None:
        d1, d2 = _disk(0, 0, 0, 10), _disk(1, 5, 0, 10)
        depth = disks_overlap(d1, d2, tol=1e-6)
        assert depth is not None
        assert depth == pytest.approx(15.0)

    def test_exactly_touching_not_flagged(self) -> None:
        d1, d2 = _disk(0, 0, 0, 10), _disk(1, 20, 0, 10)
        assert disks_overlap(d1, d2, tol=1e-6) is None

    def test_clearly_separate(self) -> None:
        d1, d2 = _disk(0, 0, 0, 10), _disk(1, 100, 0, 10)
        assert disks_overlap(d1, d2, tol=1e-6) is None


class TestDiskCapsuleOverlap:
    def test_clearly_overlapping(self) -> None:
        disk = _disk(0, 5, 0, 10)
        capsule = _capsule(0, 0, 0, 10, 0, 5)
        depth = disk_capsule_overlap(disk, capsule, tol=1e-6)
        assert depth is not None
        assert depth == pytest.approx(15.0)

    def test_exactly_touching_not_flagged(self) -> None:
        disk = _disk(0, 5, 15, 10)
        capsule = _capsule(0, 0, 0, 10, 0, 5)
        assert disk_capsule_overlap(disk, capsule, tol=1e-6) is None

    def test_clearly_separate(self) -> None:
        disk = _disk(0, 5, 100, 10)
        capsule = _capsule(0, 0, 0, 10, 0, 5)
        assert disk_capsule_overlap(disk, capsule, tol=1e-6) is None


class TestCapsuleCapsuleOverlap:
    def test_clearly_overlapping(self) -> None:
        c1 = _capsule(0, 0, 0, 10, 0, 5)
        c2 = _capsule(1, 0, 3, 10, 3, 5)
        depth = capsule_capsule_overlap(c1, c2, tol=1e-6)
        assert depth is not None
        assert depth == pytest.approx(7.0)

    def test_exactly_touching_not_flagged(self) -> None:
        c1 = _capsule(0, 0, 0, 10, 0, 5)
        c2 = _capsule(1, 0, 10, 10, 10, 5)
        assert capsule_capsule_overlap(c1, c2, tol=1e-6) is None

    def test_clearly_separate(self) -> None:
        c1 = _capsule(0, 0, 0, 10, 0, 5)
        c2 = _capsule(1, 0, 100, 10, 100, 5)
        assert capsule_capsule_overlap(c1, c2, tol=1e-6) is None


class TestAabb:
    def test_disk_aabb(self) -> None:
        disk = _disk(0, 5, 5, 10)
        assert disk.aabb == (-5, -5, 15, 15)

    def test_capsule_aabb_inflated(self) -> None:
        capsule = _capsule(0, 0, 0, 10, 0, 5)
        assert capsule.aabb == (-5, -5, 15, 5)
