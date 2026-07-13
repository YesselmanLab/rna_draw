"""Headless tests for the editor's interaction features.

Covers the three additions layered on `EditorModel` + the desktop view:
  * rotation CLAMP -- a selected helix cannot swing past `MAX_HELIX_ROTATION`
    from its laid-out direction (no unphysical flip); small rotations pass
    through unchanged;
  * selection GRANULARITY -- `select_at(index, granularity)` resolves a click
    to a single residue, the enclosing helix, or the enclosing motif;
  * a residue FREE-FORM move that re-runs the never-silent checker;
  * the MainWindow mode toolbar (Move / Select / Edit placeholder) + the
    granularity control, and that switching them changes selection behavior.

All Qt runs under ``QT_QPA_PLATFORM=offscreen``.
"""

from __future__ import annotations

import math
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest

from rna_draw.gui.model import MAX_HELIX_ROTATION, DEMO_SS, EditorModel, Selection

PySide6 = pytest.importorskip("PySide6")

from rna_draw.gui.desktop.app import build_app  # noqa: E402

_SS = DEMO_SS  # three-way junction; helix [8,17] is the second stem
_SCRATCH = (
    "/private/tmp/claude-503/-Users-jyesselman2-local-code-python-developing-rna-draw/"
    "641187e9-926c-4fa2-a1cd-2ce8215f508a/scratchpad"
)


def _helix_base_dir(model: EditorModel, start: int, end: int, pivot) -> float:
    """Direction from the pivot to the helix base midpoint (radians)."""
    cx, cy = pivot
    mx = (model._x[start] + model._x[end]) / 2.0
    my = (model._y[start] + model._y[end]) / 2.0
    return math.atan2(my - cy, mx - cx)


# -- (a) rotation clamp -----------------------------------------------------


def test_huge_rotation_is_clamped_no_flip():
    model = EditorModel.from_ss(_SS)
    start, end, pivot = model.select(8)
    ref = _helix_base_dir(model, start, end, pivot)

    model.rotate_selection(3.0)  # ~171 deg: would flip an unclamped helix

    dev = abs(math.atan2(math.sin(_helix_base_dir(model, start, end, pivot) - ref),
                         math.cos(_helix_base_dir(model, start, end, pivot) - ref)))
    assert dev <= MAX_HELIX_ROTATION + 1e-6, "helix flipped past its clamp"
    # and it landed AT (not short of) the bound, since the request exceeded it
    assert dev >= MAX_HELIX_ROTATION - math.radians(1.0)


def test_huge_negative_rotation_is_clamped():
    model = EditorModel.from_ss(_SS)
    start, end, pivot = model.select(8)
    ref = _helix_base_dir(model, start, end, pivot)
    model.rotate_selection(-3.0)
    dev = abs(math.atan2(math.sin(_helix_base_dir(model, start, end, pivot) - ref),
                         math.cos(_helix_base_dir(model, start, end, pivot) - ref)))
    assert dev <= MAX_HELIX_ROTATION + 1e-6


def test_small_rotation_unaffected_by_clamp():
    model = EditorModel.from_ss(_SS)
    start, end, pivot = model.select(8)
    assert model._clamp_rotation(0.05) == pytest.approx(0.05)
    ref = _helix_base_dir(model, start, end, pivot)
    model.rotate_selection(0.05)
    dev = _helix_base_dir(model, start, end, pivot) - ref
    assert dev == pytest.approx(0.05, abs=1e-6)


def test_incremental_rotations_accumulate_to_the_bound():
    # many small drags in one direction still cannot cross the clamp
    model = EditorModel.from_ss(_SS)
    start, end, pivot = model.select(8)
    ref = _helix_base_dir(model, start, end, pivot)
    for _ in range(60):
        model.rotate_selection(0.1)
    dev = abs(math.atan2(math.sin(_helix_base_dir(model, start, end, pivot) - ref),
                         math.cos(_helix_base_dir(model, start, end, pivot) - ref)))
    assert dev <= MAX_HELIX_ROTATION + 1e-6


# -- (b) selection granularity ----------------------------------------------


def test_select_at_residue_is_single_nt():
    model = EditorModel.from_ss(_SS)
    sel = model.select_at(5, "residue")
    assert isinstance(sel, Selection)
    assert sel.kind == "residue"
    assert sel.indices == [5]
    assert sel.pivot is None  # residues do not rotate
    assert model.sel_indices == [5]


def test_select_at_helix_is_branch_slice_with_pivot():
    model = EditorModel.from_ss(_SS)
    sel = model.select_at(8, "helix")
    assert sel.kind == "helix"
    assert (sel.start, sel.end) == (8, 17)
    assert sel.indices == list(range(8, 18))
    assert sel.pivot is not None
    # back-compat: select() still returns the tuple form + sets pivot
    assert model.select(8) == (8, 17, sel.pivot)


def test_select_at_motif_is_enclosing_loop_members():
    model = EditorModel.from_ss(_SS)
    sel = model.select_at(5, "motif")  # nt 5 is an unpaired loop nt
    assert sel.kind == "motif"
    # the enclosing loop's ring: contains the clicked nt and is non-contiguous
    assert 5 in sel.indices
    assert len(sel.indices) > 1
    assert sel.indices != list(range(sel.indices[0], sel.indices[-1] + 1))


def test_select_at_motif_on_stem_is_whole_helix():
    model = EditorModel.from_ss(_SS)
    sel = model.select_at(8, "motif")  # nt 8 is paired (a stem)
    assert sel.kind == "motif"
    # the WHOLE helix: both strands of the 4-rung stacked run (7,18)..(10,15),
    # NOT the (8,17) branch slice and NOT just the base pair.
    assert set(sel.indices) == {7, 8, 9, 10, 15, 16, 17, 18}


def test_scene_carries_selected_set():
    model = EditorModel.from_ss(_SS)
    model.select_at(5, "motif")
    scene = model.scene()
    assert scene["selection_kind"] == "motif"
    assert set(scene["selected"]) == set(model.sel_indices)


# -- (c) residue free-form move re-runs the checker -------------------------


def test_move_residue_updates_coord_and_rechecks():
    model = EditorModel.from_ss(_SS)
    model.select_at(5, "residue")
    result = model.move_residue(5, 3.0, 4.0)
    assert model._x[5] == pytest.approx(3.0)
    assert model._y[5] == pytest.approx(4.0)
    # verdict reflects reality: the returned flag matches the model's flag
    assert result.flagged == model.flagged
    assert set(result.overlaps) == set(model.scene()["overlaps"])


def test_move_residue_can_flag_overlap_honestly():
    # dropping a residue onto a neighbour must be flagged, never silent
    model = EditorModel.from_ss(_SS)
    result = model.move_residue(5, model._x[6], model._y[6])
    assert result.flagged is True
    assert result.overlaps


def test_move_residue_only_moves_one_nt():
    model = EditorModel.from_ss(_SS)
    before = list(zip(model._x, model._y))
    model.move_residue(5, model._x[5] + 2.0, model._y[5] + 2.0)
    for i, (bx, by) in enumerate(before):
        if i == 5:
            continue
        assert (model._x[i], model._y[i]) == (bx, by)


# -- (d) mode toolbar + granularity control ---------------------------------


@pytest.fixture
def win():
    _, w = build_app(ss=_SS, seq=None)
    yield w
    w.close()


def test_mode_toolbar_present_with_move_select_edit(win):
    assert win._move_act.isChecked() is True  # Move is the default
    assert win._select_act.isChecked() is False
    assert win._edit_act.isEnabled() is False  # placeholder, extensible
    # the three modes are mutually exclusive
    assert win._mode_group.isExclusive() is True


def test_granularity_control_defaults_to_helix(win):
    assert win._gran_combo.currentText() == "Helix"
    items = [win._gran_combo.itemText(i) for i in range(win._gran_combo.count())]
    assert items == ["Residue", "Helix", "Motif"]
    # default options preserve current (helix) behavior
    assert win._options_panel.options().mode == "move"
    assert win._options_panel.options().granularity == "helix"


def test_switching_mode_and_granularity_updates_options(win):
    win._select_act.trigger()
    assert win._options_panel.options().mode == "select"
    win._gran_combo.setCurrentText("Residue")
    assert win._options_panel.options().granularity == "residue"
    win._move_act.trigger()
    assert win._options_panel.options().mode == "move"


def test_select_mode_residue_click_changes_selection_behavior(win):
    # Move mode click -> helix; Select+Residue click -> single residue
    win._move_act.trigger()
    win._view._select_at(5)
    assert win._model.sel_kind == "helix"  # nt 5 resolves to enclosing helix

    win._select_act.trigger()
    win._gran_combo.setCurrentText("Residue")
    win._view._select_at(5)
    assert win._model.sel_kind == "residue"
    assert win._model.sel_indices == [5]
    # a residue selection carries no rotate handle
    assert win._view._handle is None


def test_select_mode_helix_gets_handle(win):
    win._select_act.trigger()
    win._gran_combo.setCurrentText("Helix")
    win._view._select_at(8)
    assert win._model.sel_kind == "helix"
    assert win._view._handle is not None


# -- (e) PNG renders (visual proof) -----------------------------------------


def _render_png(view, out: str) -> None:
    from PySide6 import QtGui

    img = QtGui.QImage(700, 700, QtGui.QImage.Format.Format_ARGB32)
    img.fill(QtGui.QColor("white"))
    painter = QtGui.QPainter(img)
    view._scene.render(painter)
    painter.end()
    assert img.save(out)
    assert os.path.getsize(out) > 0


def test_render_residue_selection_png(win):
    win._select_act.trigger()
    win._gran_combo.setCurrentText("Residue")
    win._view._select_at(5)
    _render_png(win._view, os.path.join(_SCRATCH, "interaction_residue_selected.png"))


def test_render_helix_clamped_png(win):
    win._move_act.trigger()
    win._view._select_at(8)
    win._model.rotate_selection(3.0)  # driven to the clamp limit
    win._view._rebuild()
    _render_png(win._view, os.path.join(_SCRATCH, "interaction_helix_clamped.png"))
