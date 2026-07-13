"""Headless tests for the editor's UNDO / REDO (full-state snapshot history).

Every editor STATE change is undoable: helix rotate, residue move, per-region
color override, highlight, sequence edit, global style change, relayout, and
version/document restore. Undo/redo are pure state restoration that RE-RUN the
never-silent native checker so a restored pose is never silently clean.

The real mouse drag/key event cannot be driven headlessly, so a gesture is
simulated exactly how the view drives it: `model.push_undo()` once, then the
mutating op(s). The MainWindow Edit>Undo/Redo actions are exercised directly.
All run under ``QT_QPA_PLATFORM=offscreen``.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import copy

import pytest

PySide6 = pytest.importorskip("PySide6")

from rna_draw.gui.desktop.app import build_app  # noqa: E402
from rna_draw.gui.model import DEMO_SEQ, DEMO_SS, EditorModel  # noqa: E402
from rna_draw.style import default_preset  # noqa: E402

_SS = DEMO_SS
_SEQ = DEMO_SEQ


@pytest.fixture
def win():
    _, w = build_app(ss=_SS, seq=_SEQ)
    yield w
    w.close()


def _fresh_model() -> EditorModel:
    return EditorModel.from_ss(_SS, seq=_SEQ)


# -- (a) rotate a helix -----------------------------------------------------


def test_rotate_undo_restores_exact_coords_redo_reapplies():
    m = _fresh_model()
    m.select(8)
    before_x, before_y = list(m._x), list(m._y)
    m.push_undo()  # simulate the view's per-gesture push at drag start
    m.rotate_selection(0.4)
    after_x, after_y = list(m._x), list(m._y)
    assert after_x != before_x or after_y != before_y  # the drag moved coords

    assert m.can_undo() is True
    assert m.undo() is True
    assert m._x == pytest.approx(before_x)
    assert m._y == pytest.approx(before_y)

    assert m.can_redo() is True
    assert m.redo() is True
    assert m._x == pytest.approx(after_x)
    assert m._y == pytest.approx(after_y)


def test_whole_drag_is_one_undo_entry():
    m = _fresh_model()
    m.select(8)
    before_x = list(m._x)
    m.push_undo()  # ONE push at drag start
    for _ in range(5):  # many per-frame mutations, NO extra pushes
        m.rotate_selection(0.05)
    assert len(m._undo) == 1
    m.undo()
    assert m._x == pytest.approx(before_x)
    assert m.can_undo() is False


# -- (b) per-region color override ------------------------------------------


def test_color_override_undo_removes_redo_restores():
    m = _fresh_model()
    m.select(8)
    assert m.color_overrides == {}
    m.push_undo()
    m.apply_color_to_selection("#123456")
    assert m.color_overrides  # override present
    # scene fill reflects the override
    fill = {nt["id"]: nt["fill"] for nt in m.scene()["nucleotides"]}
    assert fill[8] == "#123456"

    m.undo()
    assert m.color_overrides == {}
    fill = {nt["id"]: nt["fill"] for nt in m.scene()["nucleotides"]}
    assert fill[8] != "#123456"  # reverted to scheme fill

    m.redo()
    assert m.color_overrides.get(8) == "#123456"


# -- (c) highlight ----------------------------------------------------------


def test_highlight_undo_removes():
    m = _fresh_model()
    m.select(8)
    assert m.highlights == []
    m.push_undo()
    m.highlight_selection("#ffcc00")
    assert len(m.highlights) == 1

    m.undo()
    assert m.highlights == []
    m.redo()
    assert len(m.highlights) == 1


# -- (d) sequence edit ------------------------------------------------------


def test_set_sequence_undo_restores_prior_seq_and_labels():
    m = _fresh_model()
    prior_seq = m._seq
    prior_labels = [nt.get("label") for nt in m.scene()["nucleotides"]]
    new_seq = "A" * len(_SS)
    m.push_undo()
    m.set_sequence(new_seq)
    assert m._seq == new_seq
    assert [nt.get("label") for nt in m.scene()["nucleotides"]] != prior_labels

    m.undo()
    assert m._seq == prior_seq
    assert [nt.get("label") for nt in m.scene()["nucleotides"]] == prior_labels


# -- (e) global style change ------------------------------------------------


def test_set_style_undo_restores_prior_style():
    m = _fresh_model()
    prior_fill = m.scene()["style"]["default_fill"]
    new_preset = default_preset()
    new_preset.default_color = "#ff0000"
    m.push_undo()
    m.set_style(new_preset)
    assert m.scene()["style"]["default_fill"] != prior_fill

    m.undo()
    assert m.scene()["style"]["default_fill"] == prior_fill
    m.redo()
    assert m.scene()["style"]["default_fill"] != prior_fill


# -- (e2) relayout ----------------------------------------------------------


def test_relayout_undo_restores_prior_coords():
    m = _fresh_model()
    before_x, before_y = list(m._x), list(m._y)
    preset = default_preset()
    preset.layout_defaults.node_r = 16.0
    m.push_undo()
    m.relayout(preset)
    # relayout changes node_r (geometry) even if coords land similar
    assert m._node_r != pytest.approx(10.0) or m._x != pytest.approx(before_x)

    m.undo()
    assert m._x == pytest.approx(before_x)
    assert m._y == pytest.approx(before_y)
    assert m._node_r == pytest.approx(10.0)


# -- (f) mixed edits undone in reverse return to the exact start ------------


def test_mixed_edits_full_reverse_returns_to_start():
    m = _fresh_model()
    start = {
        "x": list(m._x),
        "y": list(m._y),
        "seq": m._seq,
        "overrides": dict(m.color_overrides),
        "highlights": copy.deepcopy(m.highlights),
        "fill": m.scene()["style"]["default_fill"],
    }

    # edit 1: rotate a helix
    m.select(8)
    m.push_undo()
    m.rotate_selection(0.3)
    # edit 2: color override
    m.select(23)
    m.push_undo()
    m.apply_color_to_selection("#abcdef")
    # edit 3: highlight
    m.push_undo()
    m.highlight_selection("#00ff00")
    # edit 4: sequence
    m.push_undo()
    m.set_sequence("G" * len(_SS))
    # edit 5: global style
    p = default_preset()
    p.default_color = "#010203"
    m.push_undo()
    m.set_style(p)

    # undo all five in reverse
    for _ in range(5):
        assert m.undo() is True

    assert m._x == pytest.approx(start["x"])
    assert m._y == pytest.approx(start["y"])
    assert m._seq == start["seq"]
    assert m.color_overrides == start["overrides"]
    assert m.highlights == start["highlights"]
    assert m.scene()["style"]["default_fill"] == start["fill"]
    assert m.can_undo() is False


# -- (g) redo stack clears after a new edit ---------------------------------


def test_redo_cleared_after_new_edit():
    m = _fresh_model()
    m.select(8)
    m.push_undo()
    m.rotate_selection(0.2)
    m.undo()
    assert m.can_redo() is True
    # a brand-new edit must invalidate the redo history
    m.push_undo()
    m.highlight_selection("#ffffff")
    assert m.can_redo() is False


# -- (h) empty-stack no-op + depth cap --------------------------------------


def test_undo_redo_empty_is_safe_noop():
    m = _fresh_model()
    assert m.can_undo() is False
    assert m.undo() is False  # safe no-op
    assert m.can_redo() is False
    assert m.redo() is False


def test_depth_cap_drops_oldest():
    m = _fresh_model()
    m._undo_cap = 5
    m.select(8)
    for _ in range(20):
        m.push_undo()
        m.rotate_selection(0.01)
    assert len(m._undo) == 5


def test_snapshot_is_independent_of_later_mutation():
    m = _fresh_model()
    m.select(8)
    m.push_undo()
    snap_x = list(m._undo[-1]["x"])
    m.rotate_selection(0.5)  # mutate AFTER snapshot
    assert m._undo[-1]["x"] == snap_x  # snapshot uncorrupted


# -- (i) MainWindow Edit>Undo/Redo actions ----------------------------------


def test_mainwindow_has_edit_menu_with_shortcuts(win):
    from PySide6 import QtGui

    assert win._undo_act.shortcut() == QtGui.QKeySequence(QtGui.QKeySequence.StandardKey.Undo)
    assert win._redo_act.shortcut() == QtGui.QKeySequence(QtGui.QKeySequence.StandardKey.Redo)
    # an Edit menu exists in the menu bar
    titles = [a.text() for a in win.menuBar().actions()]
    assert "&Edit" in titles


def test_mainwindow_undo_action_reverts_state_and_refreshes(win):
    m = win._model
    m.select(8)
    before_x = list(m._x)
    # simulate a committed drag gesture
    m.push_undo()
    m.rotate_selection(0.4)
    win._sync_edit_actions()
    assert win._undo_act.isEnabled() is True

    win._undo_act.trigger()
    assert m._x == pytest.approx(before_x)
    assert win._redo_act.isEnabled() is True

    win._redo_act.trigger()
    assert m._x != pytest.approx(before_x)


def test_mainwindow_color_selection_is_undoable(win):
    m = win._model
    m.select(8)
    win._on_color_selection()  # discrete op pushes undo
    assert win._undo_act.isEnabled() is True
    assert m.color_overrides
    win._undo_act.trigger()
    assert m.color_overrides == {}


def test_mainwindow_version_restore_is_undoable(win):
    m0 = win._model
    # rotate so the current pose differs from a snapshot we take now
    m0.select(8)
    coords_at_snapshot = list(m0._x)
    win._take_snapshot(name="snap")
    # now change the drawing
    m0.push_undo()
    m0.rotate_selection(0.5)
    changed = list(win._model._x)
    assert changed != pytest.approx(coords_at_snapshot)

    # restore the snapshot (undoable, in place)
    win._restore_version(0)
    assert win._model._x == pytest.approx(coords_at_snapshot)
    # undo the restore -> back to the changed pose
    win._on_undo()
    assert win._model._x == pytest.approx(changed)


# -- (j) restore re-runs the never-silent checker ---------------------------


def test_restore_reruns_checker_flag_is_honest():
    # Build a model, then drag a residue far away to CREATE an overlap so the
    # flagged state is True; undo back to the clean pose must re-run the checker
    # and report clean (not stuck on the stale flag), and redo must re-flag.
    m = _fresh_model()
    assert m.flagged is False
    m.push_undo()
    # collapse a residue onto its neighbour -> guaranteed disk overlap
    m.move_residue(5, m._x[6], m._y[6])
    assert m.flagged is True

    m.undo()
    assert m.flagged is False  # checker re-run on restore: honest clean
    assert m.scene()["overlaps"] == []

    m.redo()
    assert m.flagged is True  # checker re-run on restore: honest overlap
    assert m.scene()["overlaps"]
