"""Headless tests for the render-style switch + the RFview-like "Letters" mode.

All run under ``QT_QPA_PLATFORM=offscreen``. The editor now carries an overall
RENDER STYLE (``extra["view"]["render_style"]``): "spheres" (the classic disk
look) or "letters" (RFview-style colored letters, a gapped backbone, and typed
base-pair symbols keyed off each pair's ``bond``). This is display-only: no
coordinate or the layout ``node_r`` the overlap checker gates on ever moves, so
the never-silent contract is untouched. We assert the model tags each pair with
the right ``bond``, that Letters mode draws letters (no disks) + pair symbols,
that the switch round-trips through the document, and that MainWindow exposes a
live Render-style combo.
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
    _bond_type,
    view_style,
)
from rna_draw.style import default_preset  # noqa: E402

_DEMO_SEQ = "GGGGAAACCCCAAAGGGGAAACCCCAAAAGGGGAAACCCCA"[: len(DEMO_SS)]

# A hairpin whose four pairs are, from outside in, GC / AU / GU / non-canonical.
_BOND_SS = "((((....))))"
#            0123456789...
#  pair 0-11 = G-C (gc), 1-10 = A-U (au), 2-9 = G-U (gu), 3-8 = G-A (other)
_BOND_SEQ = "GAGGAAAAAUUC"


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    """A single QApplication for the module (RnaScene needs one to exist)."""
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    yield app


def _scratch() -> str:
    scratch = os.environ.get(
        "CLAUDE_SCRATCH",
        "/private/tmp/claude-503/-Users-jyesselman2-local-code-python-developing-rna-draw/"
        "641187e9-926c-4fa2-a1cd-2ce8215f508a/scratchpad",
    )
    os.makedirs(scratch, exist_ok=True)
    return scratch


def _disk_items(scene: RnaScene) -> list[QtWidgets.QGraphicsEllipseItem]:
    """Every nucleotide-disk ellipse item (they carry the nt index as data)."""
    return [
        it
        for it in scene.items()
        if isinstance(it, QtWidgets.QGraphicsEllipseItem)
        and it.data(_NT_INDEX_KEY) is not None
    ]


def _text_items(scene: RnaScene) -> list[QtWidgets.QGraphicsSimpleTextItem]:
    return [it for it in scene.items() if isinstance(it, QtWidgets.QGraphicsSimpleTextItem)]


def _line_items(scene: RnaScene) -> list[QtWidgets.QGraphicsLineItem]:
    return [it for it in scene.items() if isinstance(it, QtWidgets.QGraphicsLineItem)]


def _bond_of(model: EditorModel, i: int, j: int) -> str:
    for p in model.scene()["pairs"]:
        if {p["i"], p["j"]} == {i, j}:
            return p["bond"]
    raise AssertionError(f"no pair {(i, j)} in scene")


# -- (a) model tags each pair with the correct bond type ----------------------


def test_bond_type_classifier():
    assert _bond_type("G", "C") == "gc"
    assert _bond_type("C", "G") == "gc"
    assert _bond_type("A", "U") == "au"
    assert _bond_type("U", "A") == "au"
    assert _bond_type("A", "T") == "au"  # T normalized to U
    assert _bond_type("G", "U") == "gu"
    assert _bond_type("U", "G") == "gu"
    assert _bond_type("G", "A") == "other"
    assert _bond_type("C", "C") == "other"
    assert _bond_type("", "C") == "other"  # blank base
    assert _bond_type(" ", "G") == "other"


def test_scene_tags_pairs_with_bond():
    model = EditorModel.from_ss(_BOND_SS, seq=_BOND_SEQ)
    assert _bond_of(model, 0, 11) == "gc"
    assert _bond_of(model, 1, 10) == "au"
    assert _bond_of(model, 2, 9) == "gu"
    assert _bond_of(model, 3, 8) == "other"


def test_pairs_carry_bond_when_seq_blank():
    model = EditorModel.from_ss(_BOND_SS)  # no sequence
    for p in model.scene()["pairs"]:
        assert p["bond"] == "other"


# -- (b) Letters mode: no disks, but letters + pair symbols -------------------


def _letters_scene(ss: str, seq: str) -> tuple[EditorModel, RnaScene]:
    model = EditorModel.from_ss(ss, seq=seq)
    preset = default_preset()
    preset.extra["view"] = {"render_style": "letters"}
    model.set_style(preset)
    assert model.scene()["style"]["render_style"] == "letters"
    scene = RnaScene()
    scene.build(model.scene(), selection=None)
    return model, scene


def test_letters_mode_draws_letters_not_disks():
    _, scene = _letters_scene(_BOND_SS, _BOND_SEQ)
    assert _disk_items(scene) == []  # no nucleotide disks
    texts = _text_items(scene)
    assert len(texts) == len(_BOND_SEQ)  # one letter per nt
    assert _line_items(scene)  # gapped backbone + pair rungs


def test_letters_mode_pair_symbols_present():
    # gu (wobble) and other (LW placeholder) each draw a circle marker; those
    # are ellipse items WITHOUT the nt-index key, so they are not disks.
    _, scene = _letters_scene(_BOND_SS, _BOND_SEQ)
    circles = [
        it
        for it in scene.items()
        if isinstance(it, QtWidgets.QGraphicsEllipseItem)
        and it.data(_NT_INDEX_KEY) is None
    ]
    # gu + other = two marker circles (seq has no blank positions -> no dots).
    assert len(circles) == 2


def test_spheres_mode_still_draws_disks():
    model = EditorModel.from_ss(_BOND_SS, seq=_BOND_SEQ)
    assert model.scene()["style"]["render_style"] == "spheres"
    scene = RnaScene()
    scene.build(model.scene(), selection=None)
    assert _disk_items(scene)  # classic disks present


def test_letters_mode_selection_and_overlap_tint():
    # With a selection, the selected letters are tinted teal (not disks); the
    # overlap set would tint red the same way. We assert the render path runs
    # and still produces letters when a selection set is present.
    model = EditorModel.from_ss(_BOND_SS, seq=_BOND_SEQ)
    preset = default_preset()
    preset.extra["view"] = {"render_style": "letters"}
    model.set_style(preset)
    data = model.scene()
    data["selected"] = [0, 1]
    scene = RnaScene()
    scene.build(data, selection=None)
    assert _text_items(scene)
    assert _disk_items(scene) == []


# -- (c) round-trip through to_document / from_document -----------------------


def test_render_style_round_trips_through_document():
    model = EditorModel.from_ss(DEMO_SS, seq=_DEMO_SEQ)
    preset = default_preset()
    preset.extra["view"] = {"render_style": "letters"}
    model.set_style(preset)

    doc = model.to_document()
    restored = EditorModel.from_document(doc)
    assert view_style(restored.preset)["render_style"] == "letters"
    assert restored.scene()["style"]["render_style"] == "letters"


def test_default_render_style_is_spheres():
    model = EditorModel.from_ss(DEMO_SS, seq=_DEMO_SEQ)
    assert view_style(model.preset)["render_style"] == "spheres"


# -- (d) MainWindow: the Render-style combo restyles live ---------------------


@pytest.fixture
def app_window():
    app, win = build_app(ss=DEMO_SS, seq=_DEMO_SEQ)
    yield app, win
    win.close()


def test_mainwindow_has_render_style_combo(app_window):
    _, win = app_window
    panel = win._panel
    assert panel._render_style_cb is not None
    labels = [panel._render_style_cb.itemText(i) for i in range(panel._render_style_cb.count())]
    assert labels == ["Spheres", "Letters"]
    assert panel._render_style_cb.currentText() == "Spheres"


def test_render_style_combo_switches_live(app_window):
    _, win = app_window
    panel, model, view = win._panel, win._model, win._view
    assert model.scene()["style"]["render_style"] == "spheres"
    assert _disk_items(view._scene)  # spheres by default

    fired = []
    panel.visualChanged.connect(lambda: fired.append(True))
    panel._render_style_cb.setCurrentText("Letters")
    assert fired  # switching the combo emits visualChanged

    win._on_style_visual()
    assert model.scene()["style"]["render_style"] == "letters"
    assert _disk_items(view._scene) == []  # letters mode: no disks
    assert _text_items(view._scene)  # letters drawn


# -- PNGs: demo in Letters mode, and a GU + non-canonical structure -----------


def _render_png(view, out_path: str, bg: str = "#ffffff") -> None:
    scene = view._scene
    rect = scene.itemsBoundingRect().adjusted(-20, -20, 20, 20)
    img = QtGui.QImage(760, 600, QtGui.QImage.Format.Format_ARGB32)
    img.fill(QtGui.QColor(bg))
    painter = QtGui.QPainter(img)
    scene.render(painter, QtCore.QRectF(img.rect()), rect)
    painter.end()
    assert img.save(out_path)


def test_render_pngs(app_window):
    _, win = app_window
    panel, view = win._panel, win._view
    scratch = _scratch()

    # 1) The demo in Letters mode (colored letters + =/- pair symbols + thin
    #    gapped backbone).
    panel._render_style_cb.setCurrentText("Letters")
    win._on_style_visual()
    _render_png(view, os.path.join(scratch, "letters_demo.png"))

    # 2) A structure with a GU (wobble) and a non-canonical (LW) pair.
    win.load_ss(_BOND_SS, _BOND_SEQ)
    panel._render_style_cb.setCurrentText("Letters")
    win._on_style_visual()
    view.fit()
    _render_png(view, os.path.join(scratch, "letters_gu_lw.png"))
