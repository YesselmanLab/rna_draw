"""Pluggable layout engines + the checker-gated `layout_guaranteed` pipeline.

Public API:

- `LayoutEngine`: the minimal `Protocol` every backend implements.
- `LayoutResult`: what `layout_guaranteed` returns.
- `EngineError` / `EngineUnavailableError`: engine failure signals.
- `LegacyEngine` / `PuzzlerEngine` / `SafeFallbackEngine`: the original
  three backends.
- `ViennaPuzzlerEngine` / `ViennaTurtleEngine` (M5.1): in-process
  ViennaRNA layout bindings behind the same `LayoutEngine` seam.
  `PuzzlerEngine` (subprocess `RNAplot`) stays the parity oracle.
- `production_engine` (M5 productionized): in-process clearance escalation
  + a wall-clock-bounded local overlap post-pass -- the composition
  `default_engine`/`-engine auto` now prefers when the compiled ViennaRNA
  binding is available (see `rna_draw.layout.production`).
- `layout_guaranteed` / `default_engine` / `resolve_engine`: engine
  selection and the guaranteed-clean-or-flagged pipeline.
"""

from __future__ import annotations

from .base import EngineError, EngineUnavailableError, LayoutEngine, LayoutResult, RoutedLine
from .fallback import SafeFallbackEngine
from .legacy import LegacyEngine
from .pipeline import default_engine, layout_guaranteed, resolve_engine
from .production import production_engine
from .pseudoknot import layout_pseudoknot
from .puzzler import PuzzlerEngine
from .vienna import ViennaPuzzlerEngine, ViennaTurtleEngine

__all__ = [
    "LayoutEngine",
    "LayoutResult",
    "RoutedLine",
    "EngineError",
    "EngineUnavailableError",
    "LegacyEngine",
    "PuzzlerEngine",
    "SafeFallbackEngine",
    "ViennaPuzzlerEngine",
    "ViennaTurtleEngine",
    "production_engine",
    "layout_pseudoknot",
    "layout_guaranteed",
    "default_engine",
    "resolve_engine",
]
