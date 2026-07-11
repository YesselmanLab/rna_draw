"""Engine-agnostic overlap checker for a drawn RNA secondary-structure layout.

Given nucleotide coordinates, a pairing map, and geometry parameters, this
module builds the primitives a renderer actually draws (nucleotide disks,
backbone capsules, base-pair capsules), excludes the pairs that are
*supposed* to touch (backbone neighbors, base-pair partners, a disk at a
capsule's own endpoint), and reports every remaining overlap. Overlap is
judged purely by geometry: two non-adjacent, non-paired primitives that
happen to sit close together are a real, reportable overlap (see the M2
plan's Risk R2) -- there is no index-distance exclusion.
"""

from __future__ import annotations

import itertools
import statistics
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from enum import Enum
from math import hypot
from typing import Union

from .geometry import (
    Capsule,
    Disk,
    PrimitiveId,
    capsule_capsule_overlap,
    disk_capsule_overlap,
    disks_overlap,
)
from .spatial_hash import SpatialHash

Primitive = Union[Disk, Capsule]


class OverlapKind(str, Enum):
    """The three kinds of primitive-pair overlap the checker distinguishes."""

    DISK_DISK = "disk_disk"
    DISK_CAPSULE = "disk_capsule"
    CAPSULE_CAPSULE = "capsule_capsule"


@dataclass(frozen=True)
class Witness:
    """Evidence of a single overlapping primitive pair.

    Args:
        kind: Which pair of primitive types overlapped.
        id_a: Identity of the lower-ordered primitive (``id_a < id_b``).
        id_b: Identity of the higher-ordered primitive.
        separation: Actual center-to-center (disks) or axis-to-axis
            (capsules) distance between the two primitives.
        overlap_depth: How far the primitives interpenetrate, i.e. the
            required clearance minus the actual separation. Positive.
    """

    kind: OverlapKind
    id_a: PrimitiveId
    id_b: PrimitiveId
    separation: float
    overlap_depth: float


@dataclass
class OverlapReport:
    """Result of an overlap check: every witnessed overlap plus a tally.

    Args:
        witnesses: All overlaps found, sorted by ``(kind, id_a, id_b)`` for
            deterministic output.
        counts_by_kind: Number of witnesses of each `OverlapKind`.
    """

    witnesses: list[Witness]
    counts_by_kind: dict[OverlapKind, int]

    @property
    def passed(self) -> bool:
        """Whether the layout is free of overlaps."""
        return not self.witnesses

    @property
    def num_overlaps(self) -> int:
        """Total number of witnessed overlaps."""
        return len(self.witnesses)


@dataclass
class OverlapParams:
    """Geometry knobs the overlap checker judges a layout against.

    These are parameters of a *given* rendering, not a calibration to any
    particular layout engine (that is a later milestone's job) -- see the
    M2 plan's Risk R1.

    Args:
        node_r: Nucleotide disk radius, in layout (data) units. Matches
            the renderer's ``Circle(radius=NODE_R)`` (`DrawParameters`).
        backbone_half_width: Half-width of a backbone connector capsule.
            Default `0.75 * node_r` (full width `1.5 * node_r`): a
            principled layout-unit stand-in for the renderer's
            `linewidth=15` pair connectors relative to radius-10 disks --
            there is no exact points-to-data-units conversion (the two
            are drawn in different matplotlib unit systems), so this is
            stated as a parameter default, not a calibrated value.
        pair_half_width: Half-width of a base-pair connector capsule.
            Same rationale and default as `backbone_half_width`.
        tol: Tolerance subtracted from the required clearance before
            flagging an overlap; primitives that exactly touch
            (`separation == required clearance`) are not flagged.
    """

    node_r: float = 10.0
    backbone_half_width: float = 7.5
    pair_half_width: float = 7.5
    tol: float = 1e-6


def build_disks(x: Sequence[float], y: Sequence[float], node_r: float) -> list[Disk]:
    """Build one disk primitive per nucleotide.

    Args:
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.
        node_r: Disk radius shared by every nucleotide.

    Returns:
        One `Disk` per nucleotide, in nucleotide-index order.
    """
    return [
        Disk(pid=PrimitiveId("nt", i), cx=x[i], cy=y[i], radius=node_r, ends=frozenset({i}))
        for i in range(len(x))
    ]


def build_backbone_capsules(
    x: Sequence[float], y: Sequence[float], half_width: float
) -> list[Capsule]:
    """Build one capsule per consecutive backbone connection.

    Args:
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.
        half_width: Half-width of every backbone capsule.

    Returns:
        One `Capsule` per consecutive nucleotide pair ``(i, i+1)``.
    """
    return [
        Capsule(
            pid=PrimitiveId("bb", i),
            x0=x[i],
            y0=y[i],
            x1=x[i + 1],
            y1=y[i + 1],
            half_width=half_width,
            ends=frozenset({i, i + 1}),
        )
        for i in range(len(x) - 1)
    ]


def build_pair_capsules(
    x: Sequence[float], y: Sequence[float], pair_map: Sequence[int], half_width: float
) -> list[Capsule]:
    """Build one capsule per base pair.

    Args:
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.
        half_width: Half-width of every pair capsule.

    Returns:
        One `Capsule` per base pair `(i, j)` with `i < j`.
    """
    return [
        Capsule(
            pid=PrimitiveId("pair", i),
            x0=x[i],
            y0=y[i],
            x1=x[j],
            y1=y[j],
            half_width=half_width,
            ends=frozenset({i, j}),
        )
        for i, j in enumerate(pair_map)
        if j != -1 and i < j
    ]


def build_primitives(
    x: Sequence[float],
    y: Sequence[float],
    pair_map: Sequence[int],
    params: OverlapParams,
) -> list[Primitive]:
    """Build every primitive a renderer draws for this layout.

    Args:
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.
        params: Geometry parameters (disk radius, capsule half-widths).

    Returns:
        All disks, then all backbone capsules, then all pair capsules, in
        that stable order so a primitive's list index is reproducible.
    """
    disks: list[Primitive] = list(build_disks(x, y, params.node_r))
    backbone: list[Primitive] = list(build_backbone_capsules(x, y, params.backbone_half_width))
    pairs: list[Primitive] = list(build_pair_capsules(x, y, pair_map, params.pair_half_width))
    return disks + backbone + pairs


def is_excluded(a: Primitive, b: Primitive, pair_map: Sequence[int]) -> bool:
    """Decide whether a primitive pair is *supposed* to touch.

    Disk-disk pairs are excluded when backbone-adjacent (`|i - j| == 1`) or
    base-paired partners (`pair_map[i] == j`); no other disk-disk pair is
    excluded by index distance. Any pair involving at least one capsule is
    excluded exactly when the two primitives share a nucleotide endpoint
    (a disk at a capsule's own end, or two capsules meeting at a shared
    nucleotide) -- captured generically via their `ends` sets.

    Args:
        a: First primitive.
        b: Second primitive.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.

    Returns:
        True if this pair is expected to touch and should not be reported.
    """
    if isinstance(a, Disk) and isinstance(b, Disk):
        i, j = next(iter(a.ends)), next(iter(b.ends))
        return abs(i - j) == 1 or pair_map[i] == j
    return bool(a.ends & b.ends)


def _ordered(a: Primitive, b: Primitive) -> tuple[Primitive, Primitive]:
    """Normalize a mixed pair so the `Disk` comes first.

    Args:
        a: First primitive.
        b: Second primitive.

    Returns:
        `(a, b)` unless `a` is a `Capsule` and `b` is a `Disk`, in which
        case `(b, a)`.
    """
    if isinstance(a, Capsule) and isinstance(b, Disk):
        return b, a
    return a, b


def test_pair(a: Primitive, b: Primitive, tol: float) -> Witness | None:
    """Test one primitive pair for overlap and build a `Witness` if found.

    Args:
        a: First primitive.
        b: Second primitive.
        tol: Tolerance passed to the geometry overlap predicates.

    Returns:
        A `Witness` describing the overlap, or `None` if the primitives do
        not overlap beyond `tol`.

    Raises:
        TypeError: If the pair is not one of disk/disk, disk/capsule, or
            capsule/capsule (should be unreachable given `Primitive`).
    """
    first, second = _ordered(a, b)
    if isinstance(first, Disk) and isinstance(second, Disk):
        kind = OverlapKind.DISK_DISK
        required = first.radius + second.radius
        depth = disks_overlap(first, second, tol)
    elif isinstance(first, Disk) and isinstance(second, Capsule):
        kind = OverlapKind.DISK_CAPSULE
        required = first.radius + second.half_width
        depth = disk_capsule_overlap(first, second, tol)
    elif isinstance(first, Capsule) and isinstance(second, Capsule):
        kind = OverlapKind.CAPSULE_CAPSULE
        required = first.half_width + second.half_width
        depth = capsule_capsule_overlap(first, second, tol)
    else:
        raise TypeError(f"Unsupported primitive pair: {type(first)}, {type(second)}")

    if depth is None:
        return None
    id_a, id_b = sorted((a.pid, b.pid))
    return Witness(
        kind=kind, id_a=id_a, id_b=id_b, separation=required - depth, overlap_depth=depth
    )


def _validate_inputs(x: Sequence[float], y: Sequence[float], pair_map: Sequence[int]) -> None:
    """Validate coordinate/pair-map shapes before building primitives.

    Args:
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.

    Raises:
        ValueError: If lengths mismatch, inputs are empty, or `pair_map`
            is not symmetric.
    """
    if not (len(x) == len(y) == len(pair_map)):
        raise ValueError("x, y, and pair_map must all have equal length")
    if len(x) == 0:
        raise ValueError("x, y, and pair_map must be non-empty")
    for i, j in enumerate(pair_map):
        if j != -1 and pair_map[j] != i:
            raise ValueError(f"pair_map is not symmetric at index {i} (partner {j})")


def _collect_witnesses(
    primitives: Sequence[Primitive],
    candidate_pairs: Iterable[tuple[int, int]],
    pair_map: Sequence[int],
    tol: float,
) -> list[Witness]:
    """Test every candidate index pair, skipping excluded pairs.

    Args:
        primitives: All primitives, indexed as `build_primitives` produced
            them.
        candidate_pairs: Index pairs `(i, j)` to test (a superset of the
            true overlaps is sufficient -- exclusion and the geometry
            predicates filter out the rest).
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.
        tol: Tolerance passed to `test_pair`.

    Returns:
        A `Witness` for every genuinely overlapping, non-excluded pair.
    """
    witnesses = []
    for i, j in candidate_pairs:
        a, b = primitives[i], primitives[j]
        if is_excluded(a, b, pair_map):
            continue
        witness = test_pair(a, b, tol)
        if witness is not None:
            witnesses.append(witness)
    return witnesses


def _build_report(witnesses: list[Witness]) -> OverlapReport:
    """Sort witnesses deterministically and tally them by kind.

    Args:
        witnesses: Unsorted witnesses collected from a check.

    Returns:
        An `OverlapReport` with witnesses sorted by `(kind, id_a, id_b)`.
    """
    witnesses = sorted(witnesses, key=lambda w: (w.kind, w.id_a, w.id_b))
    counts: dict[OverlapKind, int] = {kind: 0 for kind in OverlapKind}
    for witness in witnesses:
        counts[witness.kind] += 1
    return OverlapReport(witnesses=witnesses, counts_by_kind=counts)


def _hash_cell_size(params: OverlapParams) -> float:
    """Pick a spatial-hash cell size sized to a typical primitive footprint.

    Args:
        params: Geometry parameters in use.

    Returns:
        `2 * (node_r + max_half_width)`, a performance knob only -- see
        `SpatialHash`'s correctness argument.
    """
    max_half_width = max(params.backbone_half_width, params.pair_half_width)
    size = 2 * (params.node_r + max_half_width)
    # Guard the degenerate all-zero-geometry call: cell_size must stay positive
    # so SpatialHash never divides by zero. Correctness is unaffected (the value
    # is a performance knob only), so any positive fallback is valid.
    return size if size > 0 else 1.0


def _insert_primitive(grid: SpatialHash, index: int, primitive: Primitive) -> None:
    """Insert one primitive into the spatial hash.

    Capsules are inserted as tiled segment pieces (`SpatialHash.insert_segment`)
    so a long diagonal chord (e.g. a base pair on a circle layout) doesn't
    blow up to one huge AABB covering `O((length / cell_size) ** 2)` cells;
    disks keep their single small AABB.

    Args:
        grid: Spatial hash to insert into.
        index: The primitive's index in the primitives list.
        primitive: The disk or capsule to insert.
    """
    if isinstance(primitive, Capsule):
        grid.insert_segment(
            index, primitive.x0, primitive.y0, primitive.x1, primitive.y1, primitive.half_width
        )
    else:
        grid.insert(index, primitive.aabb)


def check_overlaps(
    x: Sequence[float],
    y: Sequence[float],
    pair_map: Sequence[int],
    params: OverlapParams | None = None,
) -> OverlapReport:
    """Check a rendered layout for overlapping primitives (spatial-hash).

    Args:
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired (see `render_rna.get_pairmap_from_secstruct`).
        params: Geometry parameters; defaults to `OverlapParams()`.

    Returns:
        An `OverlapReport` listing every overlap between disks, backbone
        capsules, and pair capsules, excluding pairs supposed to touch.
    """
    params = params or OverlapParams()
    _validate_inputs(x, y, pair_map)
    primitives = build_primitives(x, y, pair_map, params)

    grid = SpatialHash(_hash_cell_size(params))
    for index, primitive in enumerate(primitives):
        _insert_primitive(grid, index, primitive)

    witnesses = _collect_witnesses(primitives, grid.candidate_pairs(), pair_map, params.tol)
    return _build_report(witnesses)


def check_overlaps_bruteforce(
    x: Sequence[float],
    y: Sequence[float],
    pair_map: Sequence[int],
    params: OverlapParams | None = None,
) -> OverlapReport:
    """Reference/oracle overlap check: O(n^2), no spatial hash.

    Identical overlap logic to `check_overlaps`; used in tests to confirm
    the spatial-hash version drops no true overlaps.

    Args:
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.
        params: Geometry parameters; defaults to `OverlapParams()`.

    Returns:
        An `OverlapReport` equal (as a set of witnesses) to
        `check_overlaps`'s report for the same inputs.
    """
    params = params or OverlapParams()
    _validate_inputs(x, y, pair_map)
    primitives = build_primitives(x, y, pair_map, params)

    all_pairs = itertools.combinations(range(len(primitives)), 2)
    witnesses = _collect_witnesses(primitives, all_pairs, pair_map, params.tol)
    return _build_report(witnesses)


def rescale_coords(
    x: Sequence[float], y: Sequence[float], target_backbone_step: float
) -> tuple[list[float], list[float]]:
    """Uniformly rescale coordinates to a target backbone step.

    Scales `x` and `y` so the *median* consecutive-nucleotide distance
    equals `target_backbone_step`. Useful for comparing layouts produced
    in different units (e.g. a ViennaRNA puzzler layout vs. this
    renderer's `PRIMARY_SPACE`).

    Args:
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.
        target_backbone_step: Desired median consecutive-nucleotide
            distance after rescaling.

    Returns:
        `(scaled_x, scaled_y)`.

    Raises:
        ValueError: If fewer than 2 points are given, or the median
            consecutive-nucleotide step is zero (cannot derive a scale).
    """
    n = len(x)
    if n < 2:
        raise ValueError("rescale_coords requires at least 2 points")
    steps = [hypot(x[i + 1] - x[i], y[i + 1] - y[i]) for i in range(n - 1)]
    median_step = statistics.median(steps)
    if median_step == 0:
        raise ValueError("median consecutive-nucleotide step is zero; cannot rescale")
    scale = target_backbone_step / median_step
    return [xi * scale for xi in x], [yi * scale for yi in y]


__all__ = [
    "OverlapKind",
    "Witness",
    "OverlapReport",
    "OverlapParams",
    "build_disks",
    "build_backbone_capsules",
    "build_pair_capsules",
    "build_primitives",
    "is_excluded",
    "test_pair",
    "check_overlaps",
    "check_overlaps_bruteforce",
    "rescale_coords",
]
