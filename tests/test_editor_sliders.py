"""Headless tests for the percentage size sliders + "Show spheres" toggle.

All run under ``QT_QPA_PLATFORM=offscreen``. The style panel now expresses the
nucleotide/connector sizes as PERCENTAGES of a base default (100% = the old
look) and carries a "Show spheres" switch; both are display-only, so they never
move a coordinate or the layout ``node_r`` the overlap checker gates on. We
assert the percentages map through ``view_style`` into the scene style block,
that the drawn disk radius scales (and disappears) accordingly, and that
everything round-trips through ``to_document`` / ``from_document``.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest

PySide6 = pytest.importorskip("PySide6")

from PySide6 import QtCore, QtGui, QtWidgets  # noqa: E402

from rna_draw.gui.desktop.app import build_app  # noqa: E402
from rna_draw.gui.desktop.scene_view import RnaScene, _NT_INDEX_KEY  # noqa: E402
from rna_draw.gui.model import (  # noqa: E402
    DEMO_SS,
    EditorModel,
    view_style,
)
from rna_draw.style import default_preset  # noqa: E402

_SEQ = "GGGGAAACCCCAAAGGGGAAACCCCAAAAGGGGAAACCCCA"[: len(DEMO_SS)]


@pytest.fixture
def app_window():
    app, win = build_app(ss=DEMO_SS, seq=_SEQ)
    yield app, win
    win.close()


def _disk_items(scene: RnaScene) -> list[QtWidgets.QGraphicsEllipseItem]:
    """Every nucleotide-disk ellipse item (they carry the nt index as data)."""
    return [
        it
        for it in scene.items()
        if isinstance(it, QtWidgets.QGraphicsEllipseItem)
        and it.data(_NT_INDEX_KEY) is not None
    ]


# -- (a) sphere % is a display scale: changes drawn radius, not geometry -------


def test_sphere_pct_scales_drawn_radius_without_moving_geometry(app_window):
    _, win = app_window
    panel, model, view = win._panel, win._model, win._view
    before_xy = list(zip(model._x, model._y))
    before_node_r = model._node_r

    panel._pct_sliders["sphere_pct"].setValue(60)
    win._on_style_visual()

    style = model.scene()["style"]
    assert style["sphere_scale"] == pytest.approx(0.60)
    # geometry + gate radius untouched (display-only, never a relayout)
    assert list(zip(model._x, model._y)) == before_xy
    assert model._node_r == before_node_r

    # the drawn disks are actually smaller than the layout radius
    disks = _disk_items(view._scene)
    assert disks
    drawn_r = disks[0].rect().width() / 2.0
    assert drawn_r == pytest.approx(model._node_r * 0.60)


def test_sphere_pct_larger_than_layout(app_window):
    _, win = app_window
    panel, model, view = win._panel, win._model, win._view
    panel._pct_sliders["sphere_pct"].setValue(220)
    win._on_style_visual()
    disks = _disk_items(view._scene)
    drawn_r = disks[0].rect().width() / 2.0
    assert drawn_r == pytest.approx(model._node_r * 2.20)


# -- (b) letter % and width % map through to the style block -------------------


def test_letter_and_width_pcts_map_to_style(app_window):
    _, win = app_window
    panel, model = win._panel, win._model
    panel._pct_sliders["letter_pct"].setValue(150)
    panel._pct_sliders["edge_pct"].setValue(200)
    panel._pct_sliders["pair_pct"].setValue(50)
    panel._pct_sliders["backbone_pct"].setValue(0)
    panel._pct_sliders["routed_pct"].setValue(300)
    win._on_style_visual()
    style = model.scene()["style"]
    assert style["letter_scale"] == pytest.approx(1.50)
    assert style["nt_edge_width"] == pytest.approx(1.0 * 2.0)
    assert style["pair_width"] == pytest.approx(1.6 * 0.5)
    assert style["backbone_width"] == pytest.approx(0.0)
    assert style["routed_width"] == pytest.approx(1.6 * 3.0)


# -- (c) "Show spheres" off: no disk items, but backbone/pairs/letters remain --


def test_show_spheres_off_hides_disks_only(app_window):
    _, win = app_window
    panel, model, view = win._panel, win._model, win._view
    assert model.scene()["style"]["show_spheres"] is True
    assert _disk_items(view._scene)  # on by default

    panel._spheres_cb.setChecked(False)
    win._on_style_visual()
    assert model.scene()["style"]["show_spheres"] is False

    scene = view._scene
    assert _disk_items(scene) == []  # no nucleotide disks
    # backbone (path) and pair lines still render
    paths = [it for it in scene.items() if isinstance(it, QtWidgets.QGraphicsPathItem)]
    lines = [it for it in scene.items() if isinstance(it, QtWidgets.QGraphicsLineItem)]
    texts = [
        it for it in scene.items() if isinstance(it, QtWidgets.QGraphicsSimpleTextItem)
    ]
    assert paths  # backbone
    assert lines  # base pairs
    assert texts  # letters still visible


# -- (d) round-trip through to_document / from_document ------------------------


def test_slider_view_round_trips_through_document():
    model = EditorModel.from_ss(DEMO_SS, seq=_SEQ)
    preset = default_preset()
    preset.extra["view"] = {
        "sphere_pct": 75,
        "letter_pct": 130,
        "edge_pct": 250,
        "pair_pct": 40,
        "backbone_pct": 10,
        "routed_pct": 180,
        "show_spheres": False,
    }
    model.set_style(preset)

    doc = model.to_document()
    restored = EditorModel.from_document(doc)
    style = view_style(restored.preset)
    assert style["sphere_scale"] == pytest.approx(0.75)
    assert style["letter_scale"] == pytest.approx(1.30)
    assert style["nt_edge_width"] == pytest.approx(1.0 * 2.5)
    assert style["pair_width"] == pytest.approx(1.6 * 0.4)
    assert style["backbone_width"] == pytest.approx(2.0 * 0.1)
    assert style["routed_width"] == pytest.approx(1.6 * 1.8)
    assert style["show_spheres"] is False


# -- (e) MainWindow builds with sliders/toggle; a slider drag restyles ---------


def test_mainwindow_has_sliders_and_toggle_live(app_window):
    _, win = app_window
    panel = win._panel
    for key in ("sphere_pct", "letter_pct", "edge_pct", "pair_pct", "backbone_pct", "routed_pct"):
        assert key in panel._pct_sliders
    assert panel._spheres_cb is not None

    fired = []
    panel.visualChanged.connect(lambda: fired.append(True))
    panel._pct_sliders["sphere_pct"].setValue(180)
    assert fired  # moving a slider emits visualChanged
    # and the model restyles through the wired handler
    win._on_style_visual()
    assert win._model.scene()["style"]["sphere_scale"] == pytest.approx(1.80)


# -- PNGs: spheres OFF (letters + backbone) and 60% sphere size ----------------


def _render_png(view, out_path: str, bg: str = "#ffffff") -> None:
    scene = view._scene
    rect = scene.itemsBoundingRect().adjusted(-20, -20, 20, 20)
    img = QtGui.QImage(700, 560, QtGui.QImage.Format.Format_ARGB32)
    img.fill(QtGui.QColor(bg))
    painter = QtGui.QPainter(img)
    scene.render(painter, QtCore.QRectF(img.rect()), rect)
    painter.end()
    assert img.save(out_path)


def test_render_pngs(app_window):
    _, win = app_window
    panel, view = win._panel, win._view
    scratch = os.environ.get(
        "CLAUDE_SCRATCH",
        "/private/tmp/claude-503/-Users-jyesselman2-local-code-python-developing-rna-draw/"
        "641187e9-926c-4fa2-a1cd-2ce8215f508a/scratchpad",
    )
    os.makedirs(scratch, exist_ok=True)

    panel.set_show_letters(True)
    panel._spheres_cb.setChecked(False)
    win._on_style_visual()
    _render_png(view, os.path.join(scratch, "spheres_off.png"))

    panel._spheres_cb.setChecked(True)
    panel._pct_sliders["sphere_pct"].setValue(60)
    win._on_style_visual()
    _render_png(view, os.path.join(scratch, "sphere_60pct.png"))
