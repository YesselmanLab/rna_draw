"""Pluggable layout engines + the checker-gated `layout_guaranteed` pipeline.

Public API:

- `LayoutEngine`: the minimal `Protocol` every backend implements.
- `LayoutResult`: what `layout_guaranteed` returns.
- `EngineError` / `EngineUnavailableError`: engine failure signals.
- `LegacyEngine` / `PuzzlerEngine` / `SafeFallbackEngine`: the three
  backends this milestone ships.
- `layout_guaranteed` / `default_engine` / `resolve_engine`: engine
  selection and the guaranteed-clean-or-flagged pipeline.
"""

from __future__ import annotations

from .base import EngineError, EngineUnavailableError, LayoutEngine, LayoutResult
from .fallback import SafeFallbackEngine
from .legacy import LegacyEngine
from .pipeline import default_engine, layout_guaranteed, resolve_engine
from .puzzler import PuzzlerEngine

__all__ = [
    "LayoutEngine",
    "LayoutResult",
    "EngineError",
    "EngineUnavailableError",
    "LegacyEngine",
    "PuzzlerEngine",
    "SafeFallbackEngine",
    "layout_guaranteed",
    "default_engine",
    "resolve_engine",
]
