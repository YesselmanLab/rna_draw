"""Pluggable layout engines + the checker-gated `layout_guaranteed` pipeline.

Public API:

- `LayoutEngine`: the minimal `Protocol` every backend implements.
- `LayoutResult`: what `layout_guaranteed` returns.
- `EngineError` / `EngineUnavailableError`: engine failure signals.
- `LegacyEngine` / `PuzzlerEngine` / `SafeFallbackEngine`: the original
  three backends.
- `NativePuzzlerEngine` / `NativeTurtleEngine` (Milestone A step 10): the
  owned modern-C++ `rna_layout` port behind the same `LayoutEngine` seam --
  `production_engine`'s primary since this milestone (no ViennaRNA
  dependency).
- `ViennaPuzzlerEngine` / `ViennaTurtleEngine` (M5.1): in-process ViennaRNA
  layout bindings, RETIRED to oracle/parity-only status (Milestone A step
  11) -- built only under the `RNA_DRAW_BUILD_ORACLE` CMake option;
  constructing them raises `EngineUnavailableError` in the default build.
  `PuzzlerEngine` (subprocess `RNAplot`) stays the parity oracle for
  `LegacyEngine`.
- `production_engine` (M5 productionized, native since Milestone A step
  10): in-process clearance escalation + a wall-clock-bounded local
  overlap post-pass over the native `rna_layout` engine -- the composition
  `default_engine`/`-engine auto` prefers when `_layout_core` is available
  (see `rna_draw.layout.production`).
- `layout_guaranteed` / `default_engine` / `resolve_engine`: engine
  selection and the guaranteed-clean-or-flagged pipeline.
"""

from __future__ import annotations

from .base import EngineError, EngineUnavailableError, LayoutEngine, LayoutResult, RoutedLine
from .fallback import SafeFallbackEngine
from .legacy import LegacyEngine
from .native import NativePuzzlerEngine, NativeTurtleEngine
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
    "NativePuzzlerEngine",
    "NativeTurtleEngine",
    "ViennaPuzzlerEngine",
    "ViennaTurtleEngine",
    "production_engine",
    "layout_pseudoknot",
    "layout_guaranteed",
    "default_engine",
    "resolve_engine",
]
