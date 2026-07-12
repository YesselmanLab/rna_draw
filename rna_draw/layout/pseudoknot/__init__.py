"""Pseudoknot layout (M3): multi-bracket parsing, max-nested extraction, and
checker-gated crossing-stem placement. Public entry point:
`layout_pseudoknot`, wired into `rna_draw.layout.pipeline.layout_guaranteed`
as the tier tried when `secstruct` is not pseudoknot-free.
"""

from __future__ import annotations

from .engine import layout_pseudoknot
from .extraction import max_nested_subset, stems_cross
from .parsing import Stem, group_stems, nested_secstruct, parse_all_pairs, stem_pairs
from .placement import PlacementResult, place_crossings

__all__ = [
    "layout_pseudoknot",
    "Stem",
    "parse_all_pairs",
    "group_stems",
    "stem_pairs",
    "nested_secstruct",
    "max_nested_subset",
    "stems_cross",
    "place_crossings",
    "PlacementResult",
]
