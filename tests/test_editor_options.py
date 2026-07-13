"""Headless tests for the editor's sequence input + behavior options.

Covers the three additive features:
  * `EditorModel.set_sequence` relabels/recolors without moving coordinates
    (and pads/truncates a mismatched-length sequence, returning a note);
  * the "Highlight overlaps" master switch (default OFF) gates the red tint
    while the never-silent overlap COUNT stays honest;
  * overlaps are never tinted DURING an active drag.

All Qt runs under ``QT_QPA_PLATFORM=offscreen``.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest

PySide6 = pytest.importorskip("PySide6")

from PySide6 import QtGui  # noqa: E402

from rna_draw.gui.desktop.app import build_app  # noqa: E402
from rna_draw.gui.desktop.options_panel import EditorOptions, OptionsPanel  # noqa: E402
from rna_draw.gui.desktop.scene_view import OVERLAP_RED, RnaGraphicsView, RnaScene  # noqa: E402
from rna_draw.gui.model import DEMO_SS, EditorModel  # noqa: E402

_SEQ = "GGGGAAACCCCAAAGGGGAAACCCCAAAAGGGGAAACCCCA"[: len(DEMO_SS)]
_SCRATCH = (
    "/private/tmp/claude-503/-Users-jyesselman2-local-code-python-developing-rna-draw/"
    "641187e9-926c-4fa2-a1cd-2ce8215f508a/scratchpad"
)


@pytest.fixture
def app():
    yield build_app(ss=DEMO_SS, seq=None)[0]


def _red_disk_count(qscene) -> int:
    """Number of nucleotide disks currently filled with the overlap red."""
    n = 0
    for item in qscene.items():
        brush = getattr(item, "brush", None)
        if brush is None:
            continue
        try:
            if item.brush().color().rgb() == OVERLAP_RED.rgb():
                n += 1
        except (TypeError, AttributeError):
            continue
    return n


# -- (a) sequence input -----------------------------------------------------


def test_set_sequence_updates_labels_without_moving(app):
    model = EditorModel.from_ss(DEMO_SS)
    before = list(zip(model._x, model._y))
    note = model.set_sequence(_SEQ)
    assert note == ""  # exact length, no adjustment
    labels = [nt["label"] for nt in model.scene()["nucleotides"]]
    assert "".join(labels) == _SEQ
    # coordinates are untouched -- pure relabel
    assert list(zip(model._x, model._y)) == before


def test_set_sequence_recolors_when_res_type_active(app):
    model = EditorModel.from_ss(DEMO_SS)
    model.preset.extra.setdefault("view", {})["render_type"] = "res_type"
    model.set_sequence(_SEQ)
    fills_before = [nt["fill"] for nt in model.scene()["nucleotides"]]
    model.set_sequence("A" * len(DEMO_SS))
    fills_after = [nt["fill"] for nt in model.scene()["nucleotides"]]
    # different sequence under res_type -> different per-nt fills
    assert fills_before != fills_after
    assert len(set(fills_after)) == 1  # all A -> one color


def test_set_sequence_mismatch_pads_and_notes(app):
    model = EditorModel.from_ss(DEMO_SS)
    note = model.set_sequence("ACGU")  # far shorter than the structure
    assert note  # non-empty note reported
    assert "padded" in note
    labels = [nt["label"] for nt in model.scene()["nucleotides"]]
    assert len(labels) == len(DEMO_SS)
    assert "".join(labels).startswith("ACGU")

    note2 = model.set_sequence("A" * (len(DEMO_SS) + 10))
    assert "truncated" in note2
    assert len(model._seq) == len(DEMO_SS)


def test_set_sequence_normalizes_case_and_whitespace(app):
    model = EditorModel.from_ss("((((....))))")
    model.set_sequence("ggg gaaaa cccc")  # 12 non-space chars, mixed spacing
    labels = "".join(nt["label"] for nt in model.scene()["nucleotides"])
    assert labels == "GGGGAAAACCCC"


# -- (b) highlight master switch OFF: never-silent count still honest -------


def test_highlight_off_no_red_but_count_reported(app):
    model = EditorModel.from_ss(DEMO_SS)
    model.select(8)
    result = model.rotate_selection(0.6)  # forces overlaps
    assert result.flagged and result.overlaps  # model honestly flags

    view = RnaGraphicsView()
    view.set_options(EditorOptions(highlight_overlaps=False))  # default OFF
    view.set_model(model)

    # NO disk is tinted red...
    assert _red_disk_count(view._scene) == 0
    # ...but the overlap count the status bar reads is still correct.
    assert len(model.scene()["overlaps"]) == len(result.overlaps)
    assert model.flagged is True


def test_status_bar_honest_when_highlight_off(app):
    _, win = build_app(ss=DEMO_SS, seq=None)
    # highlighting is OFF by default
    assert win._options_panel.options().highlight_overlaps is False
    # drive the model into an overlap and report it as the view would
    win._model.select(8)
    result = win._model.rotate_selection(0.6)
    win._on_overlap(result.flagged, len(result.overlaps))
    assert "overlap" in win._overlap_label.text().lower()
    win.close()


def test_highlight_on_and_idle_tints_red(app):
    model = EditorModel.from_ss(DEMO_SS)
    model.select(8)
    model.rotate_selection(0.6)

    view = RnaGraphicsView()
    view.set_options(EditorOptions(highlight_overlaps=True))
    view.set_model(model)
    assert _red_disk_count(view._scene) > 0  # highlighting on, not dragging


# -- (c) never tint during an active drag -----------------------------------


def test_no_red_during_active_drag(app):
    model = EditorModel.from_ss(DEMO_SS)
    model.select(8)
    model.rotate_selection(0.6)

    view = RnaGraphicsView()
    view.set_options(EditorOptions(highlight_overlaps=True, highlight_only_after_drag=True))
    view.set_model(model)
    assert _red_disk_count(view._scene) > 0  # idle: red shows

    # simulate an active drag
    view._rotating = True
    view._rebuild()
    assert _red_disk_count(view._scene) == 0  # suppressed mid-drag

    view._rotating = False
    view._rebuild()
    assert _red_disk_count(view._scene) > 0  # revealed once settled


def test_scene_build_show_overlaps_flag():
    model = EditorModel.from_ss(DEMO_SS)
    model.select(8)
    model.rotate_selection(0.6)
    scene_dict = model.scene()
    assert scene_dict["overlaps"]  # data still carries the offenders

    qscene = RnaScene()
    qscene.build(scene_dict, model.selection, show_overlaps=False)
    assert _red_disk_count(qscene) == 0
    qscene.build(scene_dict, model.selection, show_overlaps=True)
    assert _red_disk_count(qscene) > 0


# -- (d) window construction + PNG ------------------------------------------


def test_window_has_sequence_field_and_embedded_options(app):
    _, win = build_app(ss=DEMO_SS, seq=None)
    assert win._seq_edit is not None
    assert isinstance(win._options_panel, OptionsPanel)
    # options now live in the single LEFT dock (no separate right Options dock)
    assert not hasattr(win, "_options_dock")
    assert win._panel.section("Editing") is not None
    assert win._options_panel.options().highlight_overlaps is False
    # applying a sequence through the window relabels the disks
    win._seq_edit.setText(_SEQ)
    win._apply_sequence()
    labels = "".join(nt["label"] for nt in win._model.scene()["nucleotides"])
    assert labels == _SEQ
    win.close()


def test_render_png_letters_on_no_red(app):
    model = EditorModel.from_ss(DEMO_SS)
    model.set_sequence(_SEQ)
    model.select(8)
    model.rotate_selection(0.6)  # overlapping geometry

    view = RnaGraphicsView()
    view.set_options(EditorOptions(highlight_overlaps=False))  # no red
    view.set_model(model)
    assert _red_disk_count(view._scene) == 0

    img = QtGui.QImage(700, 700, QtGui.QImage.Format.Format_ARGB32)
    img.fill(QtGui.QColor("white"))
    painter = QtGui.QPainter(img)
    view._scene.render(painter)
    painter.end()
    out = os.path.join(_SCRATCH, "editor_options_letters_no_red.png")
    assert img.save(out)
    assert os.path.getsize(out) > 0
