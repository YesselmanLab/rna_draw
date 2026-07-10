"""`SafeFallbackEngine`: an always-clean-by-construction circle layout.

Places nucleotides at equal angles on a circle and scales the radius until
`check_overlaps` passes. This is the guaranteed-clean backstop
`rna_draw.layout.pipeline.layout_guaranteed` falls back to when a real
engine's layout fails the M2 checker, so the pipeline can honestly promise
"checker-clean, or flagged" -- never a silent overlap.

Why this terminates and is clean for any pseudoknot-free structure: the
layout is a pure uniform scaling of a fixed angular configuration, so
*every* pairwise primitive separation equals `radius * (a strictly
positive constant)`:

- no two nucleotides coincide (distinct angles);
- non-excluded pair chords never cross (well-nested pairs => non-crossing
  chords) and never cross the boundary backbone edges; backbone edges are
  convex-polygon sides (also non-crossing);
- therefore no non-excluded axes intersect and no disk center lies on a
  non-incident segment, so the minimum non-excluded separation at
  `radius=1` is a strictly positive constant `d*`.

Since the required clearances (`node_r + half_width`) are fixed, any
`radius > (max clearance) / d*` clears everything; doubling from any
positive seed reaches such a radius in `O(log)` steps. `max_doublings=40`
is astronomically safe. (For pseudoknots the crossing `[]{}` rungs are
zeroed out of `pair_map` before layout, so the loop still converges on the
primitives actually drawn; the pipeline reports `flagged=True` regardless
-- a pseudoknot is never claimed clean.)
"""

from __future__ import annotations

import math

from rna_draw.overlap import OverlapParams, check_overlaps
from rna_draw.render_rna import get_pairmap_from_secstruct

from .base import is_pseudoknot_free


class SafeFallbackEngine:
    """Lays out any structure as a circle, scaled until checker-clean."""

    name = "fallback"

    def __init__(self, params: OverlapParams | None = None, max_doublings: int = 40) -> None:
        """Store the geometry parameters to gate the scaling loop against.

        Args:
            params: Disk/capsule geometry to satisfy; defaults to
                `OverlapParams()`.
            max_doublings: Safety cap on radius-doubling iterations.
        """
        self._params = params or OverlapParams()
        self._max_doublings = max_doublings

    def layout(self, secstruct: str) -> tuple[list[float], list[float]]:
        """Lay out `secstruct` as a circle, scaled until overlap-free.

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
        """Double a seed radius until `check_overlaps` reports clean.

        Args:
            unit: Unit-circle `(cos, sin)` positions, one per nucleotide.
            pair_map: Entry `i` holds the partner index of nucleotide `i`,
                or `-1` if unpaired (pseudoknot rungs already zeroed).

        Returns:
            `(x, y)` at the first radius (or the last tried, as a
            best-effort fallback) that clears the checker.
        """
        radius = _initial_radius(len(unit), self._params)
        x, y = _scaled_positions(unit, radius)
        for _ in range(self._max_doublings):
            if check_overlaps(x, y, pair_map, self._params).passed:
                return x, y
            radius *= 2
            x, y = _scaled_positions(unit, radius)
        return x, y


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
    """Pick a positive seed radius; the doubling loop corrects it.

    Args:
        n: Number of nucleotides.
        params: Geometry parameters in use.

    Returns:
        A positive seed radius scaled with `n` and the largest clearance
        requirement, so typical structures need few doublings.
    """
    clearance = params.node_r + max(params.backbone_half_width, params.pair_half_width)
    return max(2.0, 2 * clearance) * n


__all__ = ["SafeFallbackEngine"]
