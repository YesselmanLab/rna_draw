"""Headless tests for the desktop style editor (`StylePanel` + live restyle).

All run under ``QT_QPA_PLATFORM=offscreen``: the panel mutates one
`StylePreset` in place, emits, and `MainWindow` applies the change to the
`EditorModel` (visual = no relayout) or re-runs layout (geometry). We assert
the model's scene dict reflects each change and that a preset round-trips.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest

PySide6 = pytest.importorskip("PySide6")

from PySide6 import QtCore, QtGui  # noqa: E402

from rna_draw.gui.desktop.app import build_app  # noqa: E402
from rna_draw.gui.desktop.style_panel import StylePanel  # noqa: E402
from rna_draw.gui.model import DEMO_SS, view_style  # noqa: E402
from rna_draw.style import StylePreset, default_preset  # noqa: E402

_SEQ = "GGGGAAACCCCAAAGGGGAAACCCCAAAAGGGGAAACCCCA"[: len(DEMO_SS)]


@pytest.fixture
def app_window():
    app, win = build_app(ss=DEMO_SS, seq=_SEQ)
    yield app, win
    win.close()


def test_panel_present_and_scene_has_style(app_window):
    _, win = app_window
    assert isinstance(win._panel, StylePanel)
    scene = win._model.scene()
    assert "style" in scene
    for key in ("pair_color", "background", "show_letters", "default_fill"):
        assert key in scene["style"]


def test_default_fill_change_updates_scene_fills(app_window):
    _, win = app_window
    panel, model = win._panel, win._model
    # render_type none so fills come from the default color
    panel.set_render_type("none")
    win._on_style_visual()
    panel.set_default_fill("#123456")
    win._on_style_visual()
    fills = {nt["fill"] for nt in model.scene()["nucleotides"]}
    assert fills == {"#123456"}


def test_render_type_recomputes_colors_via_colorer(app_window):
    _, win = app_window
    panel, model = win._panel, win._model
    panel.set_render_type("res_type")
    win._on_style_visual()
    scene = model.scene()
    # G -> COLORS["r"] (#ff6666), C -> COLORS["g"] (#71bc78) per the colorer
    fill_by_base = {nt["label"]: nt["fill"] for nt in scene["nucleotides"]}
    assert fill_by_base["G"] == "#ff6666"
    assert fill_by_base["C"] == "#71bc78"
    assert fill_by_base["A"] == "#ffd300"


def test_palette_edit_recolors_scheme(app_window):
    _, win = app_window
    panel, model = win._panel, win._model
    panel.set_render_type("res_type")
    win._on_style_visual()
    # recolor the 'r' palette slot (used for G in res_type) to pure black
    panel.set_palette_color("r", "#000000")
    win._on_style_visual()
    fills = {nt["label"]: nt["fill"] for nt in model.scene()["nucleotides"]}
    assert fills["G"] == "#000000"


def test_toggle_letters_updates_scene_style(app_window):
    _, win = app_window
    panel, model = win._panel, win._model
    panel.set_show_letters(False)
    win._on_style_visual()
    assert model.scene()["style"]["show_letters"] is False
    panel.set_show_letters(True)
    win._on_style_visual()
    assert model.scene()["style"]["show_letters"] is True


def test_connector_and_background_are_live_visual(app_window):
    _, win = app_window
    panel, model = win._panel, win._model
    panel._set_connector("nested_pair", "#abcdef")
    panel._set_view("background", "#202020")
    win._on_style_visual()
    style = model.scene()["style"]
    assert style["pair_color"] == "#abcdef"
    assert style["background"] == "#202020"


def test_visual_restyle_never_moves_geometry(app_window):
    _, win = app_window
    panel, model = win._panel, win._model
    before = [(nt["x"], nt["y"]) for nt in model.scene()["nucleotides"]]
    flagged_before = model.flagged
    panel.set_default_fill("#ff0000")
    panel.set_render_type("paired")
    win._on_style_visual()
    after = [(nt["x"], nt["y"]) for nt in model.scene()["nucleotides"]]
    assert before == after  # visual restyle cannot move a coordinate
    assert model.flagged == flagged_before


def test_relayout_at_new_node_r_is_checked(app_window):
    _, win = app_window
    panel, model = win._panel, win._model
    panel.set_node_r(14.0)
    win._on_style_relayout()
    scene = model.scene()
    assert scene["node_r"] > 0
    # never-silent: relayout adopts the checker's honest verdict
    assert isinstance(model.flagged, bool)
    assert len(scene["nucleotides"]) == len(DEMO_SS)


def test_reset_restores_defaults(app_window):
    _, win = app_window
    panel = win._panel
    panel.set_default_fill("#010203")
    panel._set_view("background", "#000000")
    panel.reset()
    style = view_style(panel.preset())
    default_style = view_style(default_preset())
    assert style["background"] == default_style["background"]
    # default_color reverts to the palette key 'e' (grey)
    assert panel.preset().default_color == default_preset().default_color


def test_preset_save_load_round_trip(app_window, tmp_path):
    _, win = app_window
    panel = win._panel
    panel.set_default_fill("#0a0b0c")
    panel._set_connector("nested_pair", "#111213")
    panel._set_view("background", "#141516")
    panel._set_view("pair_width", 3.5)
    panel.set_show_letters(False)
    path = tmp_path / "custom.rnastyle.json"
    panel.save_preset_file(str(path))

    reloaded = StylePreset.load(str(path))
    rstyle = view_style(reloaded)
    assert rstyle["default_fill"] == "#0a0b0c"
    assert rstyle["pair_color"] == "#111213"
    assert rstyle["background"] == "#141516"
    assert rstyle["pair_width"] == 3.5
    assert reloaded.layout_defaults.render_in_letters is False

    # round-trips back into a fresh panel too
    fresh = StylePanel()
    fresh.load_preset(reloaded)
    assert view_style(fresh.preset())["background"] == "#141516"


def test_render_restyled_scene_to_png(app_window, tmp_path):
    """Visual bonus: render a custom-palette, letters-on scene to a PNG."""
    _, win = app_window
    panel, view = win._panel, win._view
    panel.set_render_type("res_type")
    panel.set_palette_color("r", "#8e44ad")
    panel.set_show_letters(True)
    panel._set_view("background", "#101418")
    win._on_style_visual()

    scene = view._scene
    rect = scene.itemsBoundingRect().adjusted(-20, -20, 20, 20)
    img = QtGui.QImage(700, 560, QtGui.QImage.Format.Format_ARGB32)
    img.fill(QtGui.QColor("#101418"))
    painter = QtGui.QPainter(img)
    scene.render(painter, QtCore.QRectF(img.rect()), rect)
    painter.end()
    out = tmp_path / "restyled.png"
    assert img.save(str(out))
    assert out.stat().st_size > 0
