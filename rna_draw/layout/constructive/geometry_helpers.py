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


def stem_base_for_attachment(attachment: Point, axis_dir: Point, pair_space: float) -> Point:
    """The `place_stem` `base_point` that puts the near strand's rung-0
    point (`strand_a[0]`) exactly at `attachment`.

    `place_stem` centers rung 0 ON its own `base_point`, offsetting each
    strand `pair_space / 2` to either side perpendicular to `axis_dir` --
    so `strand_a[0]` naturally lands `pair_space / 2` off `base_point`, not
    on it. The constructive engine needs the OPPOSITE: a branch's
    attachment nucleotide (`closing_pair[0]`, what a sibling's backbone
    bond actually connects to) must land EXACTLY at the point its parent
    packer assigned (`envelope.py`'s soundness argument measures every
    subtree's reach from that same exact point -- an uncompensated
    `pair_space / 2` offset lets the OTHER strand of this rung, which ends
    up on the near side instead, sit close enough to a sibling's
    connecting backbone edge to be clipped by it; this is what makes that
    offset load-bearing rather than cosmetic). This computes the
    `base_point` that cancels `place_stem`'s own offset.

    Args:
        attachment: Where `closing_pair[0]`'s own coordinate must land.
        axis_dir: The `axis_dir` `place_stem` will be called with.
        pair_space: The `pair_space` `place_stem` will be called with.

    Returns:
        The `base_point` to pass to `place_stem`.
    """
    perp = rotate90_ccw(axis_dir)
    return translate(attachment, perp, pair_space / 2.0)


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

    Any leftover slack (`2 * pi` minus the angle every slot plus the
    reserved sector actually needs) is split evenly BEFORE the first slot
    and AFTER the last one, rather than dumped entirely after the last slot
    (which the naive "start right at the reserved edge" placement would
    do). This keeps the disjointness proof identical (a uniform rotation of
    every slot changes no pairwise angular gap), but matters for a
    consuming engine that also draws a straight backbone edge from OUTSIDE
    this ring (e.g. to the ring's own closing pair) to a slot's anchor: with
    all the slack on one side, a lopsided packing (one dominant slot among
    few) can leave that slot's own near edge closer, chord-wise, to the
    reserved sector via the "slack" side than the packed order suggests,
    letting such an edge approach nearly tangent to -- and graze -- the
    slot's own paired-partner disk. Splitting the slack pushes a lopsided
    packing's far slot toward sitting radially oppposite the reserved
    sector instead, so an edge reaching it arrives closer to radially
    (moving away from a tangential neighbor) rather than tangentially.

    Args:
        half_widths: Each slot's chord-clearance half-width, in fixed order.
        reserved_half_width: The reserved sector's chord-clearance half-width.
        radius: A radius for which `_total_required_angle(...) <= 2 * pi`.

    Returns:
        One center angle per slot, counterclockwise from the reserved
        sector's positive edge, each interval touching (never overlapping)
        its neighbors.
    """
    slack = 2.0 * math.pi - _total_required_angle(half_widths, reserved_half_width, radius)
    angles = []
    cursor = _angular_half_width(reserved_half_width, radius) + slack / 2.0
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


@dataclass(frozen=True)
class BulgePacking:
    """Straight-line placement plan for a single-child bulge/interior loop.

    A loop with exactly one child branch (a bulge -- unpaired bases on one
    strand -- or an interior loop -- unpaired bases on both strands) is
    placed as a near-straight continuation of the parent stem's axis
    instead of a full circular envelope (see `pack_loop_angles`): the
    reach of a CHAIN of such loops then grows linearly with chain length
    instead of compounding a fresh circular envelope (and its
    `radius_floor`-driven blow-up) at every nesting level.

    Args:
        steps: Axial step count (in `primary_space` units) from the
            closing rung to the child branch's own attachment rung --
            `max(n_before, n_after) + 1`, so BOTH rails' bulge content sits
            strictly closer to the closing rung than the child's own
            attachment, leaving it clear of the child's own subtree.
        near_axials: Axial step index (`1..n_before`) for each
            `unpaired_before` member, in structure order.
        far_axials: Axial step index for each `unpaired_after` member, in
            structure order -- descending (`n_after` down to `1`), since
            the first `unpaired_after` member (structurally adjacent to the
            child) sits closest to the child's attachment, and the last
            (structurally adjacent to the loop's far closing-pair index)
            sits closest to the closing rung.
    """

    steps: int
    near_axials: list[int]
    far_axials: list[int]


def pack_bulge_linear(n_before: int, n_after: int) -> BulgePacking:
    """Compute a single-child bulge/interior loop's straight-line axial plan.

    Args:
        n_before: Unpaired member count on the near (5') strand, before the
            child branch's own start.
        n_after: Unpaired member count on the far (3') strand, after the
            child branch's own end.

    Returns:
        The `BulgePacking` `place_bulge_geometry` places from.
    """
    steps = max(n_before, n_after) + 1
    near_axials = list(range(1, n_before + 1))
    far_axials = [n_after - k for k in range(n_after)]
    return BulgePacking(steps=steps, near_axials=near_axials, far_axials=far_axials)


def place_bulge_geometry(
    tip_a: Point,
    tip_b: Point,
    axis_dir: Point,
    packing: BulgePacking,
    primary_space: float,
    pair_space: float,
) -> tuple[list[Point], list[Point], Point]:
    """Place a bulge/interior loop's unpaired members + child attachment.

    Both rails run parallel to `axis_dir`, `pair_space` apart across the
    (virtual) rung at each axial step -- the SAME near/far offset
    convention `place_stem` uses (near = `-pair_space/2`, far =
    `+pair_space/2`), so a bulge chain composes with an ordinary stem
    ladder as one continuous pair of parallel rails.

    Args:
        tip_a: The loop's own closing pair near-strand coordinate (the
            parent stem's innermost rung, near side).
        tip_b: The loop's own closing pair far-strand coordinate.
        axis_dir: Unit vector the bulge continues along (unchanged from
            the parent stem's own axis -- the near-straight-continuation
            approach this function implements).
        packing: The axial plan from `pack_bulge_linear`.
        primary_space: Backbone step per axial unit.
        pair_space: Center-to-center rail separation.

    Returns:
        `(near_points, far_points, attachment)`: one point per
        `packing.near_axials` entry, one per `packing.far_axials` entry,
        and the point the child branch's own attachment
        (`closing_pair[0]`) must land at exactly.
    """
    perp = rotate90_ccw(axis_dir)
    half_pair = pair_space / 2.0
    center0 = midpoint(tip_a, tip_b)

    def axial_point(step: int, side: float) -> Point:
        rung_center = translate(center0, axis_dir, step * primary_space)
        return translate(rung_center, perp, side * half_pair)

    near_points = [axial_point(step, -1.0) for step in packing.near_axials]
    far_points = [axial_point(step, 1.0) for step in packing.far_axials]
    attachment = axial_point(packing.steps, -1.0)
    return near_points, far_points, attachment


def pack_line_positions(extents: list[float], primary_space: float) -> list[float]:
    """Place `len(extents)` slots along an open line, disjoint by construction.

    Unlike `pack_loop_angles` (a closed circle, used for a loop's interior),
    the exterior loop has an OPEN boundary (dangling 5'/3' tails, no
    closing pair to reserve a seam for) -- see the constructive engine's
    `_place_exterior`. Slot `k`'s center sits `max(primary_space, extents[k
    - 1] + extents[k])` past slot `k - 1`'s, strictly increasing.

    This is disjoint for EVERY pair, not just adjacent ones: for `i < j`,
    the gap `positions[j] - positions[i]` is the sum of the consecutive
    gaps between them, and that sum already includes (as two of its
    non-negative terms) the single gap `>= extents[i] + extents[j-1] >=
    ...` -- concretely, the very first term in the sum is `>= extents[i] +
    extents[i+1]` and the very last is `>= extents[j-1] + extents[j]`, so
    the total is at least `extents[i] + extents[j]` (the rest only adds).
    So a disk of radius `extents[k]` centered at `positions[k]` never
    overlaps any other slot's disk.

    Args:
        extents: Each slot's required half-width -- a disk radius, centered
            at the slot's own position, bounding everything (e.g. a whole
            child subtree's envelope) that must not collide with a
            neighboring slot's own disk.
        primary_space: Minimum consecutive-slot spacing floor, so two tiny
            extents (e.g. two bare unpaired nts) still get a normal
            backbone step rather than crowding tighter.

    Returns:
        One x position per extent, `positions[0] == 0.0`, strictly
        increasing (empty list if `extents` is empty).
    """
    positions: list[float] = []
    x = 0.0
    prev_extent = 0.0
    for index, extent in enumerate(extents):
        if index > 0:
            x += max(primary_space, prev_extent + extent)
        positions.append(x)
        prev_extent = extent
    return positions


__all__ = [
    "Point",
    "StemLadder",
    "LoopPacking",
    "BulgePacking",
    "rotate90_ccw",
    "rotate",
    "midpoint",
    "translate",
    "stem_base_for_attachment",
    "place_stem",
    "pack_loop_angles",
    "pack_bulge_linear",
    "place_bulge_geometry",
    "pack_line_positions",
    "loop_member_point",
]
