"""QApplication bootstrap + theme + entry point for the desktop editor.

`build_app` constructs (or reuses) the `QApplication` and a `MainWindow`
WITHOUT calling `exec()`, so tests can build the whole app headlessly under
``QT_QPA_PLATFORM=offscreen``. `main` is the console-script / ``python -m``
entry point that shows the window and runs the event loop.
"""

from __future__ import annotations

import sys

from PySide6 import QtWidgets

from rna_draw.gui.model import DEMO_SS

from .main_window import MainWindow


def build_app(
    ss: str = DEMO_SS, seq: str | None = None, argv: list[str] | None = None
) -> tuple[QtWidgets.QApplication, MainWindow]:
    """Build the `QApplication` and `MainWindow` without running the loop.

    Args:
        ss: Dot-bracket structure to draw on launch (defaults to the demo
            three-way junction).
        seq: Optional sequence for per-nt labels.
        argv: Argument vector for `QApplication` (defaults to `sys.argv`).

    Returns:
        ``(app, window)``; the window is NOT shown and the event loop is
        NOT started (the caller decides).
    """
    app = QtWidgets.QApplication.instance()
    if app is None:
        app = QtWidgets.QApplication(argv if argv is not None else sys.argv)
    window = MainWindow(ss=ss, seq=seq)
    # Default to the polished DARK theme (the mockup's default); this also
    # applies the QSS to the QApplication and pushes canvas colors to the view.
    window.set_theme("dark")
    return app, window


def main(argv: list[str] | None = None) -> int:
    """Launch the desktop editor (shows the window, runs the event loop)."""
    app, window = build_app(argv=argv)
    window.show()
    return app.exec()


__all__ = ["build_app", "main"]
