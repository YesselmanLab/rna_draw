"""Standalone native desktop RNA-structure editor (PySide6/Qt).

A thin Qt view over `rna_draw.gui.model.EditorModel`. Launch with::

    python -m rna_draw.gui.desktop

See `rna_draw.gui.desktop.app.main` for the entry point.
"""

from __future__ import annotations

__all__ = ["main", "build_app"]


def __getattr__(name):  # pragma: no cover - lazy Qt import
    if name in ("main", "build_app"):
        from .app import build_app, main

        return {"main": main, "build_app": build_app}[name]
    raise AttributeError(name)
