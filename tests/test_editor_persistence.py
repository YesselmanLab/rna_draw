"""Headless tests for editor SAVE/LOAD documents + version history.

Covers `EditorModel.to_document`/`from_document` (live-coord round-trip, style
and sequence fidelity, and the never-silent-on-load contract), the desktop
`.rnadoc.json` save/open flow, and the in-editor Versions snapshot list --
including that a saved file's ``extra["versions"]`` round-trips. The Qt pieces
run under ``QT_QPA_PLATFORM=offscreen``; the pure model logic needs no Qt.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest

from rna_draw.document import Document
from rna_draw.gui.model import DEMO_SEQ, DEMO_SS, EditorModel

PySide6 = pytest.importorskip("PySide6")

from rna_draw.gui.desktop.app import build_app  # noqa: E402
from rna_draw.gui.desktop.collapsible import CollapsibleSection  # noqa: E402


def _edited_model() -> EditorModel:
    """A demo model with a small (clean) hand-edit applied to helix [8, 17]."""
    model = EditorModel.from_ss(DEMO_SS, seq=DEMO_SEQ)
    model.select(8)
    model.rotate_selection(0.05)
    return model


def _overlapping_model() -> EditorModel:
    """A demo model rotated far enough to flag an overlap (never clean)."""
    model = EditorModel.from_ss(DEMO_SS, seq=DEMO_SEQ)
    model.select(8)
    model.rotate_selection(0.6)
    assert model.flagged is True  # precondition: the edit really overlaps
    return model


# -- (a) to_document -> from_document round-trips coords/seq/style -----------


def test_to_from_document_round_trips_coords_exactly():
    model = _edited_model()
    doc = model.to_document(name="v1")
    restored = EditorModel.from_document(doc)
    # LIVE edited coords are captured and restored byte-identically (no
    # re-layout, no float drift).
    assert restored._x == model._x
    assert restored._y == model._y
    assert restored._pair_map == model._pair_map
    assert restored._node_r == model._node_r


def test_to_from_document_round_trips_seq_and_style():
    model = _edited_model()
    # tweak the style so we prove it is carried, not defaulted
    model.preset.default_color = "b"
    model.preset.extra.setdefault("view", {})["render_type"] = "res_type"
    doc = model.to_document()
    restored = EditorModel.from_document(doc)
    assert restored._seq == model._seq
    assert restored.preset.to_dict() == model.preset.to_dict()
    assert doc.extra.get("name") is None  # name omitted when not supplied


# -- (b) save to file, load, coords match -----------------------------------


def test_save_load_file_round_trips_coords(tmp_path):
    model = _edited_model()
    path = tmp_path / "drawing.rnadoc.json"
    model.to_document().save(path)
    restored = EditorModel.from_document(Document.load(path))
    assert restored._x == pytest.approx(model._x)
    assert restored._y == pytest.approx(model._y)
    assert restored._seq == model._seq


# -- (c) snapshot -> edit -> restore snapshot returns to snapshot coords -----


def test_snapshot_then_edit_then_restore():
    model = _edited_model()
    snapshot = model.to_document(name="snap")
    snap_x = list(model._x)
    # further edit moves the geometry away from the snapshot
    model.select(8)
    model.rotate_selection(0.08)
    assert model._x != snap_x
    # restoring the snapshot returns to exactly its coordinates
    restored = EditorModel.from_document(snapshot)
    assert restored._x == snap_x


# -- (d) loading overlapping stored coords surfaces flagged (never-silent) ---


def test_load_overlapping_layout_is_flagged_not_silent():
    overlapping = _overlapping_model()
    doc = overlapping.to_document()
    # Tamper the ADVISORY checker cache to falsely claim the layout is clean;
    # the live gate must override it.
    doc.derived.layout.checker = {"verdict": "passed", "node_r": doc.derived.layout.node_r}
    doc.derived.layout.flagged = False

    restored = EditorModel.from_document(doc)
    assert restored.flagged is True  # never presented as silently clean
    assert restored.load_note  # a human-readable surfacing note
    assert restored.scene()["flagged"] is True
    assert restored.scene()["overlaps"]  # offending nts marked for tinting


def test_load_clean_layout_has_no_note():
    doc = _edited_model().to_document()
    restored = EditorModel.from_document(doc)
    assert restored.flagged is False
    assert restored.load_note == ""


# -- (e) saved file's extra["versions"] round-trips the snapshot list --------


def test_saved_file_versions_round_trip(tmp_path):
    _, win = build_app(ss=DEMO_SS, seq=DEMO_SEQ, argv=[])
    try:
        win._take_snapshot(name="first")
        win._model.select(8)
        win._model.rotate_selection(0.05)
        win._take_snapshot(name="second")
        assert win._versions_list.count() == 2

        path = str(tmp_path / "hist.rnadoc.json")
        win._current_path = path
        win._write_document(path)

        # the file carries the snapshot list under extra["versions"]
        doc = Document.load(path)
        assert len(doc.extra.get("versions", [])) == 2

        # reopening rebuilds the in-editor history from the file
        win._versions = []
        win._versions_list.clear()
        win._open_document(path)
        names = [win._versions_list.item(i).text() for i in range(win._versions_list.count())]
        assert names == ["first", "second"]
    finally:
        win.close()


# -- (f) MainWindow constructs with File menu + Versions section -------------


def test_window_has_file_menu_and_versions_section():
    _, win = build_app(ss=DEMO_SS, seq=None, argv=[])
    try:
        menu_titles = [m.title() for m in win.menuBar().findChildren(PySide6.QtWidgets.QMenu)]
        assert "&File" in menu_titles
        section = win._panel.section("Versions")
        assert isinstance(section, CollapsibleSection)
    finally:
        win.close()


def test_open_document_through_window_restores_layout(tmp_path):
    _, win = build_app(ss=DEMO_SS, seq=DEMO_SEQ, argv=[])
    try:
        win._model.select(8)
        win._model.rotate_selection(0.05)
        edited_x = list(win._model._x)
        path = str(tmp_path / "doc.rnadoc.json")
        win._current_path = path
        win._write_document(path)

        # load a different structure, then reopen the document
        win.load_ss("(((...)))", "GGGAAACCC")
        assert len(win._model._x) == 9
        win._open_document(path)
        assert win._model._x == pytest.approx(edited_x)
    finally:
        win.close()
