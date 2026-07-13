"""`SafeFallbackEngine`: a circle layout scaled analytically toward clean.

Places nucleotides at equal angles on a circle and scales the radius so
`check_overlaps` passes -- the guaranteed-terminating backstop
`rna_draw.layout.pipeline.layout_guaranteed` falls back to when a real
engine's layout fails the M2 checker, so the pipeline can honestly promise
"checker-clean, or flagged" -- never a silent overlap.

Why a single check pass suffices to size the radius: the layout is a pure
uniform scaling of a fixed angular configuration, so *every* pairwise
primitive separation is proportional to `radius`, while the required
clearances (`node_r + half_width`) are fixed. From one `check_overlaps`
pass, a witness with actual separation `s` and interpenetration
`overlap_depth d` needs its separation grown to the required clearance
`s + d`; because separation scales linearly with radius, multiplying the
radius by `(s + d) / s == 1 + d / s` clears *that* witness, and taking the
maximum of `1 + d / s` over all witnesses clears them all at once. So a
single analytic scale -- not an unbounded doubling loop -- reaches a clean
radius. We iterate at most `max_passes` times only to absorb floating-point
residue on the tightest witness.

The honest catch (why this is not "always clean"): a clean circle for a
large, densely paired structure (tiny hairpin loops sitting almost on top
of their own base-pair chord) can demand an astronomically large radius,
and the frozen `check_overlaps` cost grows with the radius (long chords
tile into `O(chord_length / cell_size)` pieces -- the checker's separate
O(n^2) follow-up). Rather than freeze on such a case, the radius is capped
at a cost budget (`_MAX_SEGMENT_PIECES`) that keeps every `check_overlaps`
pass affordable. Small/medium structures reach their (cheap) clean radius
under the cap and come back clean; a giant structure is capped and returns
the largest affordable radius, still *checker-run* and left for the caller
to flag -- never silently claimed clean. `SafeFallbackEngine` therefore
returns coordinates; `pipeline._fallback_result` runs the frozen checker on
them and reports `flagged=True`, so the never-silent contract holds
regardless of whether the capped circle happens to be clean.
"""

from __future__ import annotations

import math

from rna_draw.overlap import OverlapParams, check_overlaps
from rna_draw.render_rna import get_pairmap_from_secstruct

from .base import is_pseudoknot_free

# Cost budget: cap the radius so the frozen (O(chord_length/cell_size))
# checker stays affordable on every pass. A pair chord of length L tiles
# into ~L/cell_size spatial-hash pieces; summed over all pairs this is the
# dominant cost, so we bound the total tiled-piece count. Sized so a single
# `check_overlaps` pass on a several-thousand-nt circle stays around a
# couple of seconds; well above what any structure whose clean radius is
# genuinely cheap needs, so those still reach clean.
_MAX_SEGMENT_PIECES = 450_000.0

# Multiplicative safety margin on the analytic scale so the grown radius
# clears the required clearance strictly (past the checker's `tol`), rather
# than landing exactly on it.
_SCALE_MARGIN = 1.0 + 1e-6


class SafeFallbackEngine:
    """Lays out any structure as a circle, scaled analytically toward clean."""

    name = "fallback"

    def __init__(self, params: OverlapParams | None = None, max_passes: int = 3) -> None:
        """Store the geometry parameters to gate the scaling against.

        Args:
            params: Disk/capsule geometry to satisfy; defaults to
                `OverlapParams()`.
            max_passes: Bounded number of `check_overlaps` passes (each
                followed by an analytic radius solve). A single pass sizes
                the radius; the small remainder absorbs floating-point
                residue and a capped best-effort. Never an unbounded loop.
        """
        self._params = params or OverlapParams()
        self._max_passes = max(1, max_passes)

    def layout(self, secstruct: str) -> tuple[list[float], list[float]]:
        """Lay out `secstruct` as a circle, scaled toward overlap-free.

        Args:
            secstruct: Dot-bracket secondary structure.

        Returns:
            `(x, y)`, one entry per nucleotide in sequence order.
        """
        n = len(secstruct)
        if n == 0:
            return [], []
        pair_map = (
            get_pairmap_from_secstruct(secstruct) if is_pseudoknot_free(secstruct) else ([-1] * n)
        )
        unit = _unit_circle_positions(n)
        return self._scale_until_clean(unit, pair_map)

    def _scale_until_clean(
        self, unit: list[tuple[float, float]], pair_map: list[int]
    ) -> tuple[list[float], list[float]]:
        """Analytically solve the radius from bounded `check_overlaps` passes.

        Each pass runs the frozen checker once; if it passes, we are done.
        Otherwise the radius is scaled by `max(1 + overlap_depth / separation)`
        over the witnesses -- the exact factor that grows every overlapping
        separation to its required clearance in one step -- then capped at
        the cost budget so the next pass stays affordable. If the cap stops
        the radius from growing (a genuinely huge clean radius the frozen
        checker can't afford), we stop early and return the largest
        affordable radius; the pipeline runs the checker on it and flags it.

        Args:
            unit: Unit-circle `(cos, sin)` positions, one per nucleotide.
            pair_map: Entry `i` holds the partner index of nucleotide `i`,
                or `-1` if unpaired (pseudoknot rungs already zeroed).

        Returns:
            `(x, y)` at the first clean radius, or (when the clean radius
            exceeds the cost budget) at the largest affordable radius.
        """
        radius_cap = _radius_cap(unit, pair_map, self._params)
        radius = min(_initial_radius(len(unit), self._params), radius_cap)
        x, y = _scaled_positions(unit, radius)
        for _ in range(self._max_passes):
            report = check_overlaps(x, y, pair_map, self._params)
            if report.passed:
                return x, y
            scale = _analytic_scale(report.witnesses)
            next_radius = min(radius * scale, radius_cap)
            if next_radius <= radius:
                # The cost budget (or a degenerate scale) blocks further
                # growth; the current radius is the best affordable one.
                break
            radius = next_radius
            x, y = _scaled_positions(unit, radius)
        return x, y


def _analytic_scale(witnesses: list) -> float:
    """The one-step radius factor that clears every overlap witness.

    For a witness with separation `s` and interpenetration `overlap_depth d`,
    the required clearance is `s + d`; since separations scale linearly with
    the circle radius, scaling by `(s + d) / s` grows this pair to exactly
    its clearance. The maximum over all witnesses clears them all at once.

    Args:
        witnesses: The overlap witnesses from one `check_overlaps` pass.

    Returns:
        `max(1 + overlap_depth / separation)` over witnesses with a positive
        separation, times a small safety margin; `2.0` if no witness has a
        usable (positive) separation, so a degenerate pass still makes
        progress rather than stalling.
    """
    factors = [
        (w.separation + w.overlap_depth) / w.separation
        for w in witnesses
        if w.separation > 0.0
    ]
    if not factors:
        return 2.0
    return max(factors) * _SCALE_MARGIN


def _radius_cap(
    unit: list[tuple[float, float]], pair_map: list[int], params: OverlapParams
) -> float:
    """Largest radius whose `check_overlaps` pass stays within the cost budget.

    The checker tiles each pair chord into ~`chord_length / cell_size`
    spatial-hash pieces; summed over pairs this dominates its cost. At
    `radius = 1` the total chord length is `total_unit_chord`, and it scales
    linearly with the radius, so we cap the radius where the total piece
    count would reach `_MAX_SEGMENT_PIECES`.

    Args:
        unit: Unit-circle `(cos, sin)` positions, one per nucleotide.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.
        params: Geometry parameters (their `cell_size` sets the tiling).

    Returns:
        A positive radius cap; `math.inf` when there are no pair chords (the
        piece count then does not grow with the radius, so no cap is needed).
    """
    total_unit_chord = 0.0
    for i, j in enumerate(pair_map):
        if j != -1 and i < j:
            total_unit_chord += math.hypot(unit[i][0] - unit[j][0], unit[i][1] - unit[j][1])
    if total_unit_chord <= 0.0:
        return math.inf
    cell_size = 2 * (params.node_r + max(params.backbone_half_width, params.pair_half_width))
    cell_size = cell_size if cell_size > 0 else 1.0
    return _MAX_SEGMENT_PIECES * cell_size / total_unit_chord


def _unit_circle_positions(n: int) -> list[tuple[float, float]]:
    """Place `n` points at equal angles on the unit circle.

    Args:
        n: Number of nucleotides.

    Returns:
        `[(cos(2*pi*k/n), sin(2*pi*k/n)) for k in range(n)]` -- distinct
        angles give distinct points, and consecutive indices give
        consecutive angles, so the backbone traces a convex-polygon
        boundary.
    """
    return [(math.cos(2 * math.pi * k / n), math.sin(2 * math.pi * k / n)) for k in range(n)]


def _scaled_positions(
    unit: list[tuple[float, float]], radius: float
) -> tuple[list[float], list[float]]:
    """Scale unit-circle positions by `radius`.

    Args:
        unit: Unit-circle `(cos, sin)` positions.
        radius: Uniform scale factor.

    Returns:
        `(x, y)` coordinate lists.
    """
    return [radius * ux for ux, _ in unit], [radius * uy for _, uy in unit]


def _initial_radius(n: int, params: OverlapParams) -> float:
    """Pick a positive seed radius; the analytic solve corrects it.

    Args:
        n: Number of nucleotides.
        params: Geometry parameters in use.

    Returns:
        A positive seed radius scaled with `n` and the largest clearance
        requirement, so the backbone polygon starts near non-overlapping and
        the first analytic pass has few, well-separated witnesses.
    """
    clearance = params.node_r + max(params.backbone_half_width, params.pair_half_width)
    return max(2.0, 2 * clearance) * n


__all__ = ["SafeFallbackEngine"]
