"""Pure geometry for the constructive layout engine.

Two primitives compose every layout `ConstructiveEngine` builds: a
straight-helix "ladder" (`place_stem`) and a round loop's exact circular
member packing (`pack_loop_angles`). Both are pure math -- no checker or
structure-tree import -- so they can be unit-tested in isolation from the
domain layer (`rna_draw.layout.constructive.engine`), which is the only
place that verifies a composed layout against the frozen checker.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

Point = tuple[float, float]


def rotate90_ccw(vec: Point) -> Point:
    """Rotate a 2D vector 90 degrees counterclockwise.

    Args:
        vec: The vector to rotate.

    Returns:
        `vec` rotated 90 degrees counterclockwise.
    """
    return (-vec[1], vec[0])


def rotate(vec: Point, theta: float) -> Point:
    """Rotate a 2D vector by `theta` radians, counterclockwise.

    Args:
        vec: The vector to rotate.
        theta: Rotation angle in radians.

    Returns:
        `vec` rotated by `theta`.
    """
    c, s = math.cos(theta), math.sin(theta)
    return (vec[0] * c - vec[1] * s, vec[0] * s + vec[1] * c)


def midpoint(a: Point, b: Point) -> Point:
    """The midpoint of two points.

    Args:
        a: First point.
        b: Second point.

    Returns:
        `(a + b) / 2`.
    """
    return ((a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0)


def translate(point: Point, direction: Point, distance: float) -> Point:
    """Move `point` by `distance` along `direction`.

    Args:
        point: Starting point.
        direction: Unit vector to move along (not normalized here -- callers
            pass already-unit vectors, e.g. `place_stem`'s `axis_dir`).
        distance: Signed distance to move.

    Returns:
        `point + distance * direction`.
    """
    return (point[0] + distance * direction[0], point[1] + distance * direction[1])


@dataclass(frozen=True)
class StemLadder:
    """Rigid straight-ladder coordinates for one stem's paired nucleotides.

    Args:
        strand_a: Coordinates of the "near" strand (the 5' side of the
            stem's outermost pair): rung 0 (outermost, at the stem's
            `base_point`) first, rung `depth - 1` (innermost, closes the
            loop this stem leads to) last.
        strand_b: Coordinates of the "far" strand (the 3' side), same rung
            order and count; `strand_b[k]` is `strand_a[k]`'s base-pair
            partner, `pair_space` apart across the rung.
        axis_dir: Unit vector the ladder rises along, rung 0 toward rung
            `depth - 1` (base toward tip).
    """

    strand_a: list[Point]
    strand_b: list[Point]
    axis_dir: Point


def place_stem(
    depth: int,
    base_point: Point,
    axis_dir: Point,
    primary_space: float,
    pair_space: float,
) -> StemLadder:
    """Place a `depth`-pair stem as a rigid straight ladder.

    Both strands run parallel to `axis_dir`, `pair_space` apart across each
    rung; consecutive rungs are `primary_space` apart along `axis_dir` --
    the exact backbone step a renderer expects between stacked pairs.

    Args:
        depth: Number of stacked base pairs (rungs) in this stem.
        base_point: Center of rung 0 (the outermost pair), i.e. the point
            where this stem attaches to its parent loop.
        axis_dir: Unit vector the ladder rises along, from `base_point`
            toward the loop this stem closes.
        primary_space: Backbone step between consecutive rungs along
            `axis_dir`.
        pair_space: Center-to-center distance across each rung.

    Returns:
        A `StemLadder` with `depth` positions per strand.

    Raises:
        ValueError: If `depth < 1`.
    """
    if depth < 1:
        raise ValueError(f"place_stem requires depth >= 1, got {depth}")
    perp = rotate90_ccw(axis_dir)
    half_pair = pair_space / 2.0
    strand_a: list[Point] = []
    strand_b: list[Point] = []
    for k in range(depth):
        rung_center = translate(base_point, axis_dir, k * primary_space)
        strand_a.append(translate(rung_center, perp, -half_pair))
        strand_b.append(translate(rung_center, perp, half_pair))
    return StemLadder(strand_a=strand_a, strand_b=strand_b, axis_dir=axis_dir)


@dataclass(frozen=True)
class LoopPacking:
    """Result of `pack_loop_angles`: the fitted radius and each slot's angle.

    Args:
        radius: The loop-circle radius that made every slot fit.
        angles: One angle (radians, measured counterclockwise from the
            reserved sector's positive edge) per input half-width, in the
            same order, strictly increasing.
    """

    radius: float
    angles: list[float]


def _angular_half_width(half_width: float, radius: float) -> float:
    """Angle subtended by a chord-clearance `half_width` at `radius`.

    Args:
        half_width: Required clearance, as a straight-line (chord) distance
            from the loop center's radial line through the slot's anchor.
        radius: Candidate loop-circle radius.

    Returns:
        `asin(half_width / radius)`, clamped to `pi / 2` if `half_width >=
        radius` (this radius cannot fit the slot yet; the caller's
        monotone search keeps growing the radius until it can).
    """
    if radius <= 0.0:
        return math.pi
    ratio = half_width / radius
    if ratio >= 1.0:
        return math.pi / 2.0
    return math.asin(ratio)


def _total_required_angle(
    half_widths: list[float], reserved_half_width: float, radius: float
) -> float:
    """Sum of every slot's (and the reserved sector's) full angular width.

    Args:
        half_widths: Each slot's chord-clearance half-width.
        reserved_half_width: The reserved sector's chord-clearance half-width.
        radius: Candidate loop-circle radius.

    Returns:
        The total angle a full trip around the circle must accommodate at
        `radius`; the circle fits iff this is `<= 2 * pi`.
    """
    total = 2.0 * _angular_half_width(reserved_half_width, radius)
    for half_width in half_widths:
        total += 2.0 * _angular_half_width(half_width, radius)
    return total


def _angles_at(half_widths: list[float], reserved_half_width: float, radius: float) -> list[float]:
    """Slot center angles at a `radius` already known to fit.

    Args:
        half_widths: Each slot's chord-clearance half-width, in fixed order.
        reserved_half_width: The reserved sector's chord-clearance half-width.
        radius: A radius for which `_total_required_angle(...) <= 2 * pi`.

    Returns:
        One center angle per slot, counterclockwise from the reserved
        sector's positive edge, each interval touching (never overlapping)
        its neighbors.
    """
    angles = []
    cursor = _angular_half_width(reserved_half_width, radius)
    for half_width in half_widths:
        half_angle = _angular_half_width(half_width, radius)
        cursor += half_angle
        angles.append(cursor)
        cursor += half_angle
    return angles


def pack_loop_angles(
    half_widths: list[float],
    reserved_half_width: float,
    radius_floor: float,
    radius_step: float,
    max_steps: int = 100_000,
) -> LoopPacking:
    """Exact 1D circular angular-interval packing with monotone radius inflation.

    Places `len(half_widths)` slots (in fixed order) plus one reserved
    sector around a circle: each item's chord-clearance `half_width` (or
    `reserved_half_width`) subtends a full angular width `2 * asin(half_width
    / radius)`; a radius fits iff these widths sum to at most `2 * pi`.
    Searches the fixed grid `radius_floor + n * radius_step` (`n = 0, 1, 2,
    ...`) and returns the first (smallest) radius that fits -- a
    deterministic, monotone search: growing `radius` strictly shrinks every
    `asin` term, so once a radius fits, every larger radius on the grid
    fits too (this is what makes the search a plain increasing scan rather
    than an iterative solver).

    Args:
        half_widths: Each slot's required chord-clearance half-width, in
            fixed order (e.g. a loop's backbone member order).
        reserved_half_width: Chord-clearance half-width of the sector
            reserved for the loop's own closing-pair stem (`0.0` for a loop
            with no closing pair, i.e. no reservation needed).
        radius_floor: Smallest radius the grid search tries.
        radius_step: Grid step size; must be `> 0`.
        max_steps: Safety cap on grid steps (unreachable in practice: as
            `radius -> infinity`, every `asin` term `-> 0`, so the sum
            always eventually fits within `2 * pi`).

    Returns:
        `LoopPacking(radius, angles)`.

    Raises:
        ValueError: If `radius_step <= 0`.
        RuntimeError: If no radius within `max_steps` grid points fits
            (should be unreachable; see `max_steps`).
    """
    if radius_step <= 0.0:
        raise ValueError(f"pack_loop_angles requires radius_step > 0, got {radius_step}")
    for step in range(max_steps):
        radius = radius_floor + step * radius_step
        if _total_required_angle(half_widths, reserved_half_width, radius) <= 2.0 * math.pi:
            return LoopPacking(radius, _angles_at(half_widths, reserved_half_width, radius))
    raise RuntimeError(
        f"pack_loop_angles: no fitting radius found within {max_steps} steps from "
        f"{radius_floor} (step {radius_step})"
    )


def loop_member_point(center: Point, zero_dir: Point, radius: float, angle: float) -> Point:
    """A point on a loop's circle at `angle` counterclockwise from `zero_dir`.

    Args:
        center: Loop circle center.
        zero_dir: Unit vector defining angle 0 (the reserved sector's seam,
            i.e. the direction back toward the stem that closes this loop).
        radius: Circle radius.
        angle: Angle in radians, counterclockwise from `zero_dir`.

    Returns:
        `center + radius * rotate(zero_dir, angle)`.
    """
    direction = rotate(zero_dir, angle)
    return translate(center, direction, radius)


__all__ = [
    "Point",
    "StemLadder",
    "LoopPacking",
    "rotate90_ccw",
    "rotate",
    "midpoint",
    "translate",
    "place_stem",
    "pack_loop_angles",
    "loop_member_point",
]
