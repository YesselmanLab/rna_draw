"""Pluggable layout engines + the checker-gated `layout_guaranteed` pipeline.

Public API:

- `LayoutEngine`: the minimal `Protocol` every backend implements.
- `LayoutResult`: what `layout_guaranteed` returns.
- `EngineError` / `EngineUnavailableError`: engine failure signals.
- `LegacyEngine` / `PuzzlerEngine` / `SafeFallbackEngine`: the original
  three backends.
- `ViennaPuzzlerEngine` / `ViennaNaviewEngine` / `ViennaTurtleEngine`
  (M5.1): in-process ViennaRNA bindings behind the same `LayoutEngine`
  seam. `PuzzlerEngine` (subprocess `RNAplot`) stays the parity oracle and
  the public default; these are additive engine choices only.
- `layout_guaranteed` / `default_engine` / `resolve_engine`: engine
  selection and the guaranteed-clean-or-flagged pipeline.
"""

from __future__ import annotations

from .base import EngineError, EngineUnavailableError, LayoutEngine, LayoutResult
from .fallback import SafeFallbackEngine
from .legacy import LegacyEngine
from .pipeline import default_engine, layout_guaranteed, resolve_engine
from .puzzler import PuzzlerEngine
from .vienna import ViennaNaviewEngine, ViennaPuzzlerEngine, ViennaTurtleEngine

__all__ = [
    "LayoutEngine",
    "LayoutResult",
    "EngineError",
    "EngineUnavailableError",
    "LegacyEngine",
    "PuzzlerEngine",
    "SafeFallbackEngine",
    "ViennaPuzzlerEngine",
    "ViennaNaviewEngine",
    "ViennaTurtleEngine",
    "layout_guaranteed",
    "default_engine",
    "resolve_engine",
]
