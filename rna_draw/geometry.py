"""Geometric primitives and distance math for the overlap checker.

Two shapes are drawn on an RNA layout: nucleotide "disks" (circles) and
"capsules" (thick line segments, i.e. a backbone or base-pair connector
inflated by a half-width). This module defines those primitives and the
pure distance/overlap math the checker needs. No engine, layout, or I/O
logic lives here -- only geometry.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import hypot


@dataclass(frozen=True, order=True)
class PrimitiveId:
    """Stable identity for a drawn primitive.

    Orderable (by ``kind`` then ``index``) so witnesses can report a
    deterministic ``id_a < id_b``.

    Args:
        kind: Primitive family, e.g. "nt" (nucleotide disk), "bb" (backbone
            capsule), or "pair" (base-pair capsule).
        index: Index within that family (typically the lower nucleotide
            index for capsules, or the nucleotide index for disks).
    """

    kind: str
    index: int


@dataclass(frozen=True)
class Disk:
    """A nucleotide glyph: a circle of a given radius.

    Args:
        pid: Stable identity of this disk.
        cx: Center x-coordinate.
        cy: Center y-coordinate.
        radius: Disk radius.
        ends: Nucleotide indices this disk represents; always a singleton
            ``{i}`` for the disk drawn at nucleotide ``i``.
    """

    pid: PrimitiveId
    cx: float
    cy: float
    radius: float
    ends: frozenset[int]

    @property
    def aabb(self) -> tuple[float, float, float, float]:
        """Axis-aligned bounding box as ``(minx, miny, maxx, maxy)``."""
        return (
            self.cx - self.radius,
            self.cy - self.radius,
            self.cx + self.radius,
            self.cy + self.radius,
        )


@dataclass(frozen=True)
class Capsule:
    """A thick segment: a backbone or base-pair connector.

    Args:
        pid: Stable identity of this capsule.
        x0: Start point x-coordinate.
        y0: Start point y-coordinate.
        x1: End point x-coordinate.
        y1: End point y-coordinate.
        half_width: Half the stroke width (inflation radius around the
            segment axis).
        ends: The two nucleotide indices this capsule connects.
    """

    pid: PrimitiveId
    x0: float
    y0: float
    x1: float
    y1: float
    half_width: float
    ends: frozenset[int]

    @property
    def aabb(self) -> tuple[float, float, float, float]:
        """Segment bounding box inflated by ``half_width`` on all sides."""
        w = self.half_width
        return (
            min(self.x0, self.x1) - w,
            min(self.y0, self.y1) - w,
            max(self.x0, self.x1) + w,
            max(self.y0, self.y1) + w,
        )


def clamp(value: float, lo: float, hi: float) -> float:
    """Clamp ``value`` to the closed interval ``[lo, hi]``.

    Args:
        value: Value to clamp.
        lo: Lower bound.
        hi: Upper bound.

    Returns:
        ``value`` restricted to ``[lo, hi]``.
    """
    if value < lo:
        return lo
    if value > hi:
        return hi
    return value


def point_segment_distance(
    px: float, py: float, ax: float, ay: float, bx: float, by: float
) -> float:
    """Shortest distance between a point and a segment.

    Projects the point onto the segment's supporting line, clamps the
    projection parameter to ``[0, 1]`` so it stays on the segment, then
    returns the Euclidean distance to that clamped closest point.

    Args:
        px: Point x-coordinate.
        py: Point y-coordinate.
        ax: Segment start x-coordinate.
        ay: Segment start y-coordinate.
        bx: Segment end x-coordinate.
        by: Segment end y-coordinate.

    Returns:
        Distance from ``(px, py)`` to the segment ``(a, b)``.
    """
    abx, aby = bx - ax, by - ay
    denom = abx * abx + aby * aby
    if denom == 0.0:
        # Degenerate segment (a == b): distance reduces to point-point.
        return hypot(px - ax, py - ay)
    t = ((px - ax) * abx + (py - ay) * aby) / denom
    t = clamp(t, 0.0, 1.0)
    return hypot(px - (ax + t * abx), py - (ay + t * aby))


def _closest_params(
    ax: float, ay: float, abx: float, aby: float, cx: float, cy: float, cdx: float, cdy: float
) -> tuple[float, float]:
    """Clamped parametric solution for the closest points of two segments.

    Segment 1 is ``a + s*ab`` for ``s in [0, 1]``; segment 2 is
    ``c + t*cd`` for ``t in [0, 1]``. Follows the clamped-parametric
    approach in Ericson, *Real-Time Collision Detection* (2005),
    section 5.1.9, "ClosestPtSegmentSegment".

    Args:
        ax: Segment 1 start x.
        ay: Segment 1 start y.
        abx: Segment 1 direction x (``b - a``).
        aby: Segment 1 direction y (``b - a``).
        cx: Segment 2 start x.
        cy: Segment 2 start y.
        cdx: Segment 2 direction x (``d - c``).
        cdy: Segment 2 direction y (``d - c``).

    Returns:
        Tuple ``(s, t)`` of clamped parameters in ``[0, 1]``.
    """
    r_x, r_y = ax - cx, ay - cy
    a_dot_a = abx * abx + aby * aby
    c_dot_c = cdx * cdx + cdy * cdy
    c_dot_r = cdx * r_x + cdy * r_y
    a_dot_r = abx * r_x + aby * r_y
    a_dot_c = abx * cdx + aby * cdy

    denom = a_dot_a * c_dot_c - a_dot_c * a_dot_c
    s = 0.0
    if denom != 0.0:
        s = clamp((a_dot_c * c_dot_r - a_dot_r * c_dot_c) / denom, 0.0, 1.0)
    t = (a_dot_c * s + c_dot_r) / c_dot_c if c_dot_c != 0.0 else 0.0
    if t < 0.0:
        t = 0.0
        s = clamp(-a_dot_r / a_dot_a, 0.0, 1.0) if a_dot_a != 0.0 else 0.0
    elif t > 1.0:
        t = 1.0
        s = clamp((a_dot_c - a_dot_r) / a_dot_a, 0.0, 1.0) if a_dot_a != 0.0 else 0.0
    return s, t


def segment_segment_distance(
    p1: tuple[float, float],
    q1: tuple[float, float],
    p2: tuple[float, float],
    q2: tuple[float, float],
) -> float:
    """Shortest distance between two segments ``p1-q1`` and ``p2-q2``.

    Reference: Ericson, *Real-Time Collision Detection* (2005), section
    5.1.9, "ClosestPtSegmentSegment". Degenerate (zero-length) segments are
    handled as special cases up front so the general clamped-parametric
    solve only ever runs on two genuine segments.

    Args:
        p1: Start of segment 1.
        q1: End of segment 1.
        p2: Start of segment 2.
        q2: End of segment 2.

    Returns:
        Distance between the closest points of the two segments.
    """
    abx, aby = q1[0] - p1[0], q1[1] - p1[1]
    cdx, cdy = q2[0] - p2[0], q2[1] - p2[1]
    seg1_degenerate = abx == 0.0 and aby == 0.0
    seg2_degenerate = cdx == 0.0 and cdy == 0.0

    if seg1_degenerate and seg2_degenerate:
        return hypot(p1[0] - p2[0], p1[1] - p2[1])
    if seg1_degenerate:
        return point_segment_distance(p1[0], p1[1], p2[0], p2[1], q2[0], q2[1])
    if seg2_degenerate:
        return point_segment_distance(p2[0], p2[1], p1[0], p1[1], q1[0], q1[1])

    s, t = _closest_params(p1[0], p1[1], abx, aby, p2[0], p2[1], cdx, cdy)
    closest1 = (p1[0] + s * abx, p1[1] + s * aby)
    closest2 = (p2[0] + t * cdx, p2[1] + t * cdy)
    return hypot(closest1[0] - closest2[0], closest1[1] - closest2[1])


def disks_overlap(d1: Disk, d2: Disk, tol: float) -> float | None:
    """Test whether two disks overlap.

    Args:
        d1: First disk.
        d2: Second disk.
        tol: Tolerance subtracted from the required clearance; exact
            touching (``dist == r1 + r2``) is not flagged.

    Returns:
        Penetration depth ``(r1 + r2 - dist)`` if the disks overlap beyond
        ``tol``, otherwise ``None``.
    """
    dist = hypot(d1.cx - d2.cx, d1.cy - d2.cy)
    required = d1.radius + d2.radius
    if dist < required - tol:
        return required - dist
    return None


def disk_capsule_overlap(d: Disk, c: Capsule, tol: float) -> float | None:
    """Test whether a disk overlaps a capsule.

    Args:
        d: The disk.
        c: The capsule.
        tol: Tolerance subtracted from the required clearance.

    Returns:
        Penetration depth if the disk overlaps the capsule beyond ``tol``,
        otherwise ``None``.
    """
    dist = point_segment_distance(d.cx, d.cy, c.x0, c.y0, c.x1, c.y1)
    required = d.radius + c.half_width
    if dist < required - tol:
        return required - dist
    return None


def capsule_capsule_overlap(c1: Capsule, c2: Capsule, tol: float) -> float | None:
    """Test whether two capsules overlap.

    Args:
        c1: First capsule.
        c2: Second capsule.
        tol: Tolerance subtracted from the required clearance.

    Returns:
        Penetration depth if the capsules overlap beyond ``tol``, otherwise
        ``None``.
    """
    dist = segment_segment_distance((c1.x0, c1.y0), (c1.x1, c1.y1), (c2.x0, c2.y0), (c2.x1, c2.y1))
    required = c1.half_width + c2.half_width
    if dist < required - tol:
        return required - dist
    return None
