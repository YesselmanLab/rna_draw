"""Phase 3 PK-A candidate-route geometry: `placement._route_pair`'s quality
ladder (see that module) for routing a non-overlapping, ORTHOGONAL-ONLY
line from nucleotide `i` to `j` when a straight in-plane connector (PK-B)
is rejected.

STYLE (locked, user-directed): every crossing-bond segment is STRAIGHT and
AXIS-ALIGNED (horizontal or vertical) -- never a diagonal, never a curve.
The standard depiction here is a 3-segment ORTHOGONAL "STAPLE": from `i`
go straight OUT (perpendicular to the structure, along one axis), straight
ACROSS (the other axis) to `j`'s own out-leg, then straight back IN to
`j` -- a bracket/staple bridging the two nucleotides, e.g.:

    i x---+           +---x j
          |           |
          +-----------+

Tier 0 (`find_direct_route`): the bare `[p_i, p_j]` chord, but ONLY when
it is ALREADY axis-aligned (`p_i`/`p_j` share an x or a y) -- a diagonal
direct chord is never accepted, however clean.

Tier 1 (`find_staple_route`): the 3-segment staple above, escalating OFFSET
(how far the shared rail sits from the structure, as a fraction of the
layout's own extent) and DIRECTION (up/down/left/right, i.e. which axis is
the out/in leg and which side of the chord the rail sits on) -- cheap and
sufficient for the vast majority of crossing pairs. Processing pairs
innermost-first with a consistent direction preference (`_direction_order`,
nearest cardinal to "away from center" tried first) naturally nests a
multi-bp crossing stem's staples into a tidy parallel family, like nested
brackets -- see `placement._route_stem`'s docstring.

Tier 2 (`find_offset_staple_route`): a 5-segment "staggered" staple -- each
endpoint FIRST takes a short axis-aligned local jog (perpendicular to its
own out-leg) before heading to the shared rail. Catches an endpoint whose
immediate neighbors block every plain (unjogged) out-leg at every offset,
without ever leaving the axis-aligned grid.

Every tier generates AND checks candidates (via `validate`), returning the
first clean `RoutedLine` or `None`; `placement.py` only orchestrates which
tier to try, in what order, and what to do once every tier fails (drop the
pair, flagged, never a silent or drawn overlap).
"""

from __future__ import annotations

from math import hypot

from rna_draw.geometry import Capsule
from rna_draw.layout.base import RoutedLine
from rna_draw.overlap import OverlapParams, Primitive

from .validate import polyline_capsules, polyline_is_clean

Point = tuple[float, float]
Direction = str  # one of "up", "down", "left", "right"

# Tier 1/2 offset ladder: how far the shared rail sits from the chord's own
# span, as a fraction of the layout's bounding-box diagonal (`extent`).
# Mirrors the pre-existing bow-magnitude ladder's range (0.05..12.8):
# small steps cover the common short/local crossing cheaply, the large
# steps are a near-guaranteed escape for a pair embedded deep inside a big
# structure. MEASURED: a much denser ladder (15 steps up to 20.0) was
# tried and made no difference on the real corpus -- the bottleneck for
# the residual unplaced pairs is having only 4 candidate directions
# (axis-aligned, vs. the old continuous bearing fan), not ladder
# resolution, so the simpler 9-step ladder is kept.
_OFFSET_FACTORS: tuple[float, ...] = (0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8)

# Tier 2 local-jog ladder (`find_offset_staple_route`): short, in units of
# `extent`, tried both signs (`_JOG_SIGNS`) independently per endpoint.
_JOG_FACTORS: tuple[float, ...] = (0.05, 0.15)
_JOG_SIGNS: tuple[float, ...] = (1.0, -1.0)

_DIRECTIONS: tuple[Direction, ...] = ("up", "down", "left", "right")
_OUTER_MARGIN_FACTOR = 1.5
_AXIS_ALIGN_TOL = 1e-6


def enclosing_circle(x: list[float], y: list[float], params: OverlapParams) -> tuple[Point, float]:
    """A circle, centered on the layout's bbox center, that encloses every
    primitive the layout can draw -- used only for `center` (the staple
    ladder's shared reference point for direction preference).
    """
    center = ((max(x) + min(x)) / 2, (max(y) + min(y)) / 2)
    half_diag = hypot(max(x) - center[0], max(y) - center[1])
    margin = _OUTER_MARGIN_FACTOR * (
        params.node_r + max(params.backbone_half_width, params.pair_half_width)
    )
    return center, max(half_diag + margin, params.node_r)


def find_direct_route(
    i: int,
    j: int,
    p_i: Point,
    p_j: Point,
    params: OverlapParams,
    base_primitives: list[Primitive],
    committed_lines: list[Capsule],
    pair_map: list[int],
) -> RoutedLine | None:
    """Tier 0 (best case): the single straight `i -> j` segment -- ONLY if
    it is already axis-aligned (`p_i`/`p_j` share an x or a y). A diagonal
    direct chord is never accepted, however checker-clean, per the
    no-diagonals style contract.
    """
    if abs(p_i[0] - p_j[0]) > _AXIS_ALIGN_TOL and abs(p_i[1] - p_j[1]) > _AXIS_ALIGN_TOL:
        return None
    return _accept_if_clean(i, j, [p_i, p_j], params, base_primitives, committed_lines, pair_map)


def find_staple_route(
    i: int,
    j: int,
    p_i: Point,
    p_j: Point,
    center: Point,
    extent: float,
    params: OverlapParams,
    base_primitives: list[Primitive],
    committed_lines: list[Capsule],
    pair_map: list[int],
) -> RoutedLine | None:
    """Tier 1: the 3-segment orthogonal staple, escalating (direction, offset).

    Directions are tried in `_direction_order`'s consistent preference
    (nearest cardinal to "away from center" first) so sibling pairs of the
    same crossing stem tend to pick the SAME direction, nesting into a
    nested-bracket family instead of scattering.
    """
    for direction in _direction_order(center, p_i, p_j):
        for factor in _OFFSET_FACTORS:
            points = _staple_points(p_i, p_j, direction, factor * extent)
            line = _accept_if_clean(
                i, j, points, params, base_primitives, committed_lines, pair_map
            )
            if line is not None:
                return line
    return None


def find_offset_staple_route(
    i: int,
    j: int,
    p_i: Point,
    p_j: Point,
    center: Point,
    extent: float,
    params: OverlapParams,
    base_primitives: list[Primitive],
    committed_lines: list[Capsule],
    pair_map: list[int],
) -> RoutedLine | None:
    """Tier 2: a 5-segment staggered staple -- each endpoint independently
    jogs sideways (perpendicular to its own out-leg) before heading to the
    shared rail, catching an endpoint whose plain out-leg (tier 1) is
    blocked by an immediate neighbor at every offset.
    """
    jogs = [sign * factor * extent for factor in _JOG_FACTORS for sign in _JOG_SIGNS]
    for direction in _direction_order(center, p_i, p_j):
        for factor in _OFFSET_FACTORS:
            offset = factor * extent
            for jog_i in jogs:
                for jog_j in jogs:
                    points = _offset_staple_points(p_i, p_j, direction, offset, jog_i, jog_j)
                    line = _accept_if_clean(
                        i, j, points, params, base_primitives, committed_lines, pair_map
                    )
                    if line is not None:
                        return line
    return None


def _direction_order(center: Point, p_i: Point, p_j: Point) -> list[Direction]:
    """Cardinal directions, nearest to "away from `center`" first.

    Tried in this consistent order for every pair, so sibling pairs of the
    same crossing stem tend to converge on the SAME direction (nesting
    into a parallel bracket family) rather than picking directions
    independently at random.
    """
    mid = ((p_i[0] + p_j[0]) / 2, (p_i[1] + p_j[1]) / 2)
    dx, dy = mid[0] - center[0], mid[1] - center[1]
    primary = _nearest_cardinal(dx, dy)
    return [primary] + [d for d in _DIRECTIONS if d != primary]


def _nearest_cardinal(dx: float, dy: float) -> Direction:
    """The cardinal direction nearest to the `(dx, dy)` bearing."""
    if dx == 0.0 and dy == 0.0:
        return "up"
    if abs(dx) >= abs(dy):
        return "right" if dx >= 0 else "left"
    return "up" if dy >= 0 else "down"


def _staple_points(p_i: Point, p_j: Point, direction: Direction, offset: float) -> list[Point]:
    """The 4 vertices (3 segments) of an orthogonal staple `p_i -> p_j`.

    `direction` picks BOTH which axis the out/in legs travel along and
    which side of the chord the shared rail sits on (see the module
    docstring's ASCII diagram): "up"/"down" -> vertical out/in legs, a
    horizontal rail; "left"/"right" -> horizontal out/in legs, a vertical
    rail.
    """
    if direction in ("up", "down"):
        rail = _rail(p_i[1], p_j[1], direction == "up", offset)
        return [p_i, (p_i[0], rail), (p_j[0], rail), p_j]
    rail = _rail(p_i[0], p_j[0], direction == "right", offset)
    return [p_i, (rail, p_i[1]), (rail, p_j[1]), p_j]


def _offset_staple_points(
    p_i: Point, p_j: Point, direction: Direction, offset: float, jog_i: float, jog_j: float
) -> list[Point]:
    """The 6 vertices (5 segments) of a staggered staple: each endpoint
    jogs sideways by `jog_i`/`jog_j` (perpendicular to its own out-leg)
    before joining the shared rail -- see `find_offset_staple_route`.
    """
    if direction in ("up", "down"):
        rail = _rail(p_i[1], p_j[1], direction == "up", offset)
        start_i, start_j = (p_i[0] + jog_i, p_i[1]), (p_j[0] + jog_j, p_j[1])
        return [p_i, start_i, (start_i[0], rail), (start_j[0], rail), start_j, p_j]
    rail = _rail(p_i[0], p_j[0], direction == "right", offset)
    start_i, start_j = (p_i[0], p_i[1] + jog_i), (p_j[0], p_j[1] + jog_j)
    return [p_i, start_i, (rail, start_i[1]), (rail, start_j[1]), start_j, p_j]


def _rail(coord_i: float, coord_j: float, positive_side: bool, offset: float) -> float:
    """The shared rail coordinate `offset` past the far side of `(coord_i, coord_j)`."""
    if positive_side:
        return max(coord_i, coord_j) + offset
    return min(coord_i, coord_j) - offset


def _accept_if_clean(
    i: int,
    j: int,
    points: list[Point],
    params: OverlapParams,
    base_primitives: list[Primitive],
    committed_lines: list[Capsule],
    pair_map: list[int],
) -> RoutedLine | None:
    """Build `points`' capsules and return a `RoutedLine` iff they're clean."""
    segments = polyline_capsules(points, i, j, params.pair_half_width, line_uid=0)
    if polyline_is_clean(segments, base_primitives, committed_lines, pair_map, params.tol):
        return RoutedLine(i=i, j=j, points=points)
    return None


__all__ = [
    "enclosing_circle",
    "find_direct_route",
    "find_staple_route",
    "find_offset_staple_route",
]
