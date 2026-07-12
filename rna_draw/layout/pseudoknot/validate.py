"""PK-A polyline validation (Phase 3).

Builds `Capsule` primitives for a candidate routed line and tests them
against the current layout's own primitives + any already-committed PK-A
lines, reusing the frozen checker's own exclusion rule and geometry
predicates UNCHANGED (`overlap.is_excluded`, `geometry.disk_capsule_
overlap`, `geometry.capsule_capsule_overlap`). This module never edits
`overlap.py`/`geometry.py` semantics; it only calls their public building
blocks -- see the plan's "How the frozen checker validates crossing
elements".
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

from rna_draw.geometry import (
    Capsule,
    Disk,
    PrimitiveId,
    capsule_capsule_overlap,
    disk_capsule_overlap,
)
from rna_draw.overlap import Primitive, is_excluded

Point = tuple[float, float]

# Synthetic capsule "ends" ids for a polyline's interior vertices, offset
# well clear of any real nucleotide index (structures top out at
# `constructive.engine._MAX_NUCLEOTIDES` == 4000) and namespaced per line
# via `line_uid`, so `is_excluded`'s shared-endpoint rule (reused
# UNCHANGED) correctly treats two consecutive segments of the SAME
# polyline as "supposed to touch" at their shared joint -- exactly like two
# backbone capsules sharing a nucleotide.
_SYNTHETIC_ID_BASE = 1_000_000


def polyline_capsules(
    points: Sequence[Point], i: int, j: int, half_width: float, line_uid: int
) -> list[Capsule]:
    """Build one `Capsule` per consecutive pair of a candidate PK-A polyline.

    Args:
        points: Sampled polyline points; `points[0]` sits at nucleotide
            `i`'s position and `points[-1]` at nucleotide `j`'s.
        i: 5' nucleotide index this line connects.
        j: 3' nucleotide index this line connects.
        half_width: Capsule half-width (matches the render stroke).
        line_uid: A value unique across every polyline built during one
            `layout_pseudoknot` call -- only used to keep synthetic
            interior-vertex ids from colliding across different lines.

    Returns:
        `len(points) - 1` capsules; the first carries `i` as one of its
        `ends`, the last carries `j`, and every interior joint carries a
        synthetic id shared by its two adjacent segments.
    """
    ids = [_vertex_id(k, len(points), i, j, line_uid) for k in range(len(points))]
    return [
        Capsule(
            pid=PrimitiveId("pkline", line_uid * 1000 + k),
            x0=points[k][0],
            y0=points[k][1],
            x1=points[k + 1][0],
            y1=points[k + 1][1],
            half_width=half_width,
            ends=frozenset({ids[k], ids[k + 1]}),
        )
        for k in range(len(points) - 1)
    ]


def _vertex_id(k: int, n_points: int, i: int, j: int, line_uid: int) -> int:
    """The real nucleotide id at a polyline's two endpoints, else synthetic."""
    if k == 0:
        return i
    if k == n_points - 1:
        return j
    return -(_SYNTHETIC_ID_BASE + line_uid * 1000 + k)


def polyline_is_clean(
    segments: Sequence[Capsule],
    base_primitives: Sequence[Primitive],
    committed_lines: Sequence[Capsule],
    pair_map: Sequence[int],
    tol: float,
) -> bool:
    """Whether every segment of a candidate polyline clears the layout.

    Tests each segment against every base primitive (disks, backbone and
    pair capsules from the current nested+PK-B layout), every already
    committed PK-A line's segments, and its own non-adjacent sibling
    segments (a self-crossing guard) -- via `is_excluded` +
    `disk_capsule_overlap`/`capsule_capsule_overlap`, unchanged.

    Args:
        segments: The candidate polyline's own capsules (from
            `polyline_capsules`).
        base_primitives: Every disk/backbone-capsule/pair-capsule the
            current layout draws (`overlap.build_primitives`'s output).
        committed_lines: Segments of every PK-A line already accepted.
        pair_map: The current augmented pair_map; only reachable by
            `is_excluded`'s disk-disk branch, which never fires here since
            one side of every tested pair is always a `Capsule`, but
            required by its signature.
        tol: Tolerance passed to the geometry predicates.

    Returns:
        True iff no segment overlaps anything beyond `tol`.
    """
    others: list[Primitive] = [*base_primitives, *committed_lines]
    for index, segment in enumerate(segments):
        siblings = [s for pos, s in enumerate(segments) if pos != index]
        if not capsule_is_clean(segment, others, pair_map, tol):
            return False
        if not capsule_is_clean(segment, siblings, pair_map, tol):
            return False
    return True


def capsule_is_clean(
    segment: Capsule, candidates: Iterable[Primitive], pair_map: Sequence[int], tol: float
) -> bool:
    """Whether `segment` overlaps none of `candidates` (excluded pairs skipped).

    Exposed (not just `polyline_is_clean`-internal) so callers can also
    probe a SINGLE candidate capsule cheaply -- e.g. `placement._route_pair`
    testing one local "exit leg" before assembling a full candidate route.

    Args:
        segment: The capsule to test.
        candidates: Primitives to test `segment` against.
        pair_map: Passed through to `is_excluded`; see `polyline_is_clean`.
        tol: Tolerance passed to the geometry predicates.

    Returns:
        True iff `segment` overlaps none of `candidates` beyond `tol`.
    """
    for other in candidates:
        if is_excluded(segment, other, pair_map):
            continue
        if _overlap_depth(segment, other, tol) is not None:
            return False
    return True


def _overlap_depth(segment: Capsule, other: Primitive, tol: float) -> float | None:
    """Dispatch to the right frozen geometry predicate for `other`'s type."""
    if isinstance(other, Disk):
        return disk_capsule_overlap(other, segment, tol)
    return capsule_capsule_overlap(segment, other, tol)


__all__ = ["polyline_capsules", "polyline_is_clean", "capsule_is_clean"]
