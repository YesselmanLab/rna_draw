"""Headless tests for the floating selection inspector popover.

Covers the contextual `FloatingSelectionPanel` that appears near a live
selection with inline styling actions:
  * show on a non-empty selection with the right readout; hide when cleared;
  * its Color action recolors the selection (scene fill) and is undoable;
  * its Highlight action adds a highlight halo;
  * its Clear action drops color + highlight;
  * `place_near` keeps the panel inside the viewport bounds;
  * MainWindow builds with the popover and a theme toggle restyles it.

A real mouse cannot run headlessly, so selection is driven through the model /
`_on_selection` and the panel's action handlers are invoked directly. All Qt
runs under ``QT_QPA_PLATFORM=offscreen``.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest

from rna_draw.gui.model import DEMO_SEQ, DEMO_SS

PySide6 = pytest.importorskip("PySide6")

from PySide6 import QtCore  # noqa: E402

from rna_draw.gui.desktop.app import build_app  # noqa: E402
from rna_draw.gui.desktop.floating_selection import FloatingSelectionPanel  # noqa: E402


def _fills(model) -> dict[int, str]:
    return {nt["id"]: nt["fill"] for nt in model.scene()["nucleotides"]}


# -- (a) show near a selection with a readout; hide when cleared ------------


def test_selection_shows_floating_panel_with_readout():
    _, win = build_app(ss=DEMO_SS, seq=DEMO_SEQ, argv=[])
    try:
        # isVisible() needs a shown ancestor (the window is never shown in the
        # headless tests), so assert the panel's own shown/hidden state.
        assert win._floating.isHidden()
        win._model.select_indices([0, 1, 2, 3, 4])
        win._on_selection((0, 4))
        assert not win._floating.isHidden()
        assert "5 nt" in win._floating._readout.text()
    finally:
        win.close()


def test_clearing_selection_hides_floating_panel():
    _, win = build_app(ss=DEMO_SS, seq=DEMO_SEQ, argv=[])
    try:
        win._model.select_indices([0, 1, 2])
        win._on_selection((0, 2))
        assert not win._floating.isHidden()
        win._model.deselect()
        win._on_selection(None)
        assert win._floating.isHidden()
    finally:
        win.close()


def test_helix_readout_names_the_helix():
    text = FloatingSelectionPanel._readout_text("helix", [6, 7, 8, 9, 10, 11])
    assert text.startswith("Helix 7")  # 1-based
    assert "6 nt" in text


# -- (b) the Color action recolors the selection and is undoable ------------


def test_floating_color_recolors_and_is_undoable():
    _, win = build_app(ss=DEMO_SS, seq=DEMO_SEQ, argv=[])
    try:
        win._model.select_indices([0, 1, 2, 3])
        win._on_selection((0, 3))
        baseline = _fills(win._model)
        win._floating._color_btn.set_color("#ff0000")
        # emit the panel's action exactly as the button click would
        win._floating.colorRequested.emit(win._floating._color_btn.color())
        win._model.deselect()
        fills = _fills(win._model)
        assert all(fills[i] == "#ff0000" for i in (0, 1, 2, 3))
        # undoable: reverts to the pre-color scheme fills
        assert win._model.undo()
        assert _fills(win._model) == baseline
    finally:
        win.close()


# -- (c) the Highlight action adds a highlight halo -------------------------


def test_floating_highlight_adds_highlight():
    _, win = build_app(ss=DEMO_SS, seq=DEMO_SEQ, argv=[])
    try:
        win._model.select_indices([8, 9, 10, 11])
        win._on_selection((8, 11))
        win._floating._highlight_btn.set_color("#ffd300")
        win._floating.highlightRequested.emit(win._floating._highlight_btn.color())
        assert win._model.scene()["highlights"] == [
            {"indices": [8, 9, 10, 11], "color": "#ffd300"}
        ]
        # and undoable
        assert win._model.undo()
        assert "highlights" not in win._model.scene()
    finally:
        win.close()


def test_floating_clear_drops_color_and_highlight():
    _, win = build_app(ss=DEMO_SS, seq=DEMO_SEQ, argv=[])
    try:
        win._model.select_indices([0, 1, 2])
        win._on_selection((0, 2))
        win._floating.colorRequested.emit("#00ff00")
        win._floating.highlightRequested.emit("#00ffff")
        assert win._model.color_overrides
        assert win._model.highlights
        # re-select the same set and clear
        win._model.select_indices([0, 1, 2])
        win._floating.clearRequested.emit()
        assert win._model.color_overrides == {}
        assert win._model.highlights == []
    finally:
        win.close()


# -- (d) place_near keeps the panel inside the viewport ---------------------


def test_place_near_stays_within_viewport():
    panel = FloatingSelectionPanel()
    viewport = QtCore.QRect(0, 0, 800, 600)
    for anchor in (
        QtCore.QRect(400, 300, 60, 40),   # center
        QtCore.QRect(2, 2, 30, 30),       # top-left corner
        QtCore.QRect(770, 570, 40, 40),   # bottom-right corner
        QtCore.QRect(390, 0, 40, 20),     # flush against the top edge
    ):
        panel.place_near(anchor, viewport)
        geo = panel.geometry()
        assert geo.left() >= viewport.left()
        assert geo.top() >= viewport.top()
        assert geo.right() <= viewport.right()
        assert geo.bottom() <= viewport.bottom()
    panel.deleteLater()


def test_view_selection_rect_maps_to_viewport():
    _, win = build_app(ss=DEMO_SS, seq=DEMO_SEQ, argv=[])
    try:
        win._model.select_indices([0, 1, 2, 3])
        win._view.restyle()
        rect = win._view.selection_view_rect()
        assert rect is not None
        assert rect.width() >= 0 and rect.height() >= 0
        # no selection -> no rect
        win._model.deselect()
        assert win._view.selection_view_rect() is None
    finally:
        win.close()


# -- (e) MainWindow builds with the popover; theme toggle restyles it -------


def test_window_builds_with_floating_panel():
    _, win = build_app(ss=DEMO_SS, seq=DEMO_SEQ, argv=[])
    try:
        assert isinstance(win._floating, FloatingSelectionPanel)
        # child of the view's viewport so it floats over the canvas
        assert win._floating.parent() is win._view.viewport()
    finally:
        win.close()


def test_theme_toggle_restyles_floating_panel():
    _, win = build_app(ss=DEMO_SS, seq=DEMO_SEQ, argv=[])
    try:
        win.set_theme("dark")
        dark_qss = win._floating.styleSheet()
        win.set_theme("light")
        light_qss = win._floating.styleSheet()
        assert dark_qss and light_qss
        assert dark_qss != light_qss
    finally:
        win.close()
