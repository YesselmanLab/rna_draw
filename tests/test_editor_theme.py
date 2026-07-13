"""Headless tests for the desktop editor's visual theme overhaul.

Covers the dark/light QSS themes + live toggle, the dot-grid canvas
``drawBackground``, and the clearly-visible teal box-select rectangle. All run
under ``QT_QPA_PLATFORM=offscreen``; nothing here exercises editor behavior
(selection/rotate/undo/save) -- those live in the other editor tests.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest

PySide6 = pytest.importorskip("PySide6")

from PySide6 import QtCore, QtGui, QtWidgets  # noqa: E402

from rna_draw.gui.desktop.app import build_app  # noqa: E402
from rna_draw.gui.desktop.theme import canvas_colors, qss_for  # noqa: E402
from rna_draw.gui.model import DEMO_SS  # noqa: E402


@pytest.fixture
def win():
    _, w = build_app(ss=DEMO_SS, seq=None)
    yield w
    w.close()


# -- (a) themes apply, switching is error-free ------------------------------


def test_default_theme_is_dark(win):
    assert win.current_theme() == "dark"


def test_set_theme_applies_non_empty_stylesheet(win):
    app = QtWidgets.QApplication.instance()
    win.set_theme("dark")
    assert app.styleSheet().strip()
    win.set_theme("light")
    assert app.styleSheet().strip()
    # switching back and forth must not raise
    win.set_theme("dark")
    win.set_theme("light")
    assert win.current_theme() == "light"


def test_toggle_action_switches_theme(win):
    win.set_theme("dark")
    win._theme_act.setChecked(True)  # checked == light
    assert win.current_theme() == "light"
    win._theme_act.setChecked(False)
    assert win.current_theme() == "dark"


# -- (b) the two themes differ ----------------------------------------------


def test_dark_and_light_stylesheets_differ(win):
    app = QtWidgets.QApplication.instance()
    win.set_theme("dark")
    dark = app.styleSheet()
    win.set_theme("light")
    light = app.styleSheet()
    assert dark != light
    assert qss_for("dark") != qss_for("light")
    # the canvas base color also flips between themes
    assert canvas_colors("dark")[0] != canvas_colors("light")[0]


# -- (c) drawBackground (dot grid) renders offscreen in each theme ----------


def test_drawbackground_renders_in_both_themes(win):
    view = win._view
    for theme in ("dark", "light"):
        win.set_theme(theme)
        img = QtGui.QImage(320, 240, QtGui.QImage.Format.Format_ARGB32)
        img.fill(0)
        painter = QtGui.QPainter(img)
        view.render(painter)  # exercises RnaGraphicsView.drawBackground
        painter.end()
        # the themed canvas base color must have been painted
        expected = QtGui.QColor(canvas_colors(theme)[0])
        assert img.pixelColor(3, 3) == expected


# -- (d) the box-select rectangle is visible during a drag, gone after ------


def test_box_select_rect_is_visible_then_removed(win):
    view = win._view
    view._options.mode = "select"
    assert view._sel_rect is None
    view._box_begin(QtCore.QPoint(20, 20), additive=False)
    assert view._sel_rect is not None
    assert view._sel_rect.scene() is not None  # actually in the scene (visible)
    # teal, translucent fill + a teal border == clearly visible
    assert view._sel_rect.brush().color().alpha() > 0
    assert view._sel_rect.pen().color().alpha() > 0
    view._box_update(QtCore.QPoint(240, 200))
    assert view._sel_rect.rect().width() > 0
    assert view._sel_rect.rect().height() > 0
    view._box_finish()
    assert view._sel_rect is None  # removed after the drag settles
