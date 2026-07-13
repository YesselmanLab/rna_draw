"""Headless tests for the overhauled desktop editor UI.

Covers the UI-only changes: the single LEFT dock of collapsible sections (no
right Options dock), the folding behavior, the moved-out Structure & Sequence
inputs applied via Draw, the pre-filled demo sequence (letters on launch),
VARNA-style residue numbering, and that the render-type combo carries its full
item text (the cutoff fix). All run under ``QT_QPA_PLATFORM=offscreen``.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest

PySide6 = pytest.importorskip("PySide6")

from rna_draw.gui.desktop.app import build_app  # noqa: E402
from rna_draw.gui.desktop.collapsible import CollapsibleSection  # noqa: E402
from rna_draw.gui.model import DEMO_SEQ, DEMO_SS, EditorModel  # noqa: E402

_SECTION_TITLES = [
    "Structure & Sequence",
    "Nucleotides",
    "Connectors",
    "Colors / scheme",
    "Canvas",
    "Editing",
]


@pytest.fixture
def win():
    _, w = build_app(ss=DEMO_SS, seq=None)
    yield w
    w.close()


# -- (1) one left dock, no right Options dock -------------------------------


def test_single_left_dock_no_right_options_dock(win):
    # exactly one dock widget, docked LEFT
    from PySide6 import QtCore, QtWidgets

    docks = win.findChildren(QtWidgets.QDockWidget)
    assert len(docks) == 1
    assert win.dockWidgetArea(win._dock) == QtCore.Qt.DockWidgetArea.LeftDockWidgetArea
    # the old right Options dock is gone
    assert not hasattr(win, "_options_dock")


def test_all_sections_present(win):
    for title in _SECTION_TITLES:
        sec = win._panel.section(title)
        assert isinstance(sec, CollapsibleSection), f"missing section {title!r}"


def test_default_fold_state(win):
    # Structure & Sequence + Nucleotides open; the rest folded
    assert win._panel.section("Structure & Sequence").is_expanded() is True
    assert win._panel.section("Nucleotides").is_expanded() is True
    for title in ("Connectors", "Colors / scheme", "Canvas", "Editing"):
        assert win._panel.section(title).is_expanded() is False


# -- (2) folding hides/shows content ----------------------------------------


def test_collapsible_toggles_content_visibility(win):
    sec = win._panel.section("Nucleotides")
    content = sec.content_widget()
    assert sec.is_expanded() is True
    assert content.isHidden() is False
    sec.set_expanded(False)
    assert sec.is_expanded() is False
    assert content.isHidden() is True  # content explicitly hidden when folded
    sec.set_expanded(True)
    assert content.isHidden() is False


def test_multiple_sections_open_at_once(win):
    win._panel.section("Connectors").set_expanded(True)
    assert win._panel.section("Nucleotides").is_expanded() is True
    assert win._panel.section("Connectors").is_expanded() is True


# -- (2b) structure + sequence moved out of the toolbar, applied via Draw ----


def test_structure_and_sequence_apply_via_draw(win):
    ss = "((((....))))"
    seq = "GGGGAAAACCCC"
    win._ss_edit.setText(ss)
    win._seq_edit.setText(seq)
    win._draw_from_edit()
    scene = win._model.scene()
    assert len(scene["nucleotides"]) == len(ss)
    labels = "".join(nt["label"] for nt in scene["nucleotides"])
    assert labels == seq


def test_returnpressed_triggers_draw(win):
    win._ss_edit.setText("(((...)))")
    win._seq_edit.setText("GGGAAACCC")
    win._ss_edit.returnPressed.emit()
    assert len(win._model.scene()["nucleotides"]) == 9


def test_structure_fields_not_in_toolbar(win):
    from PySide6 import QtWidgets

    bars = win.findChildren(QtWidgets.QToolBar)
    toolbar_edits = []
    for bar in bars:
        toolbar_edits.extend(bar.findChildren(QtWidgets.QLineEdit))
    assert toolbar_edits == []  # ss/seq line edits left the toolbar


# -- (4) demo starts with a sequence (letters on launch) --------------------


def test_demo_starts_with_sequence(win):
    # the Sequence field is pre-filled and the disks carry the demo letters
    assert win._seq_edit.text() == DEMO_SEQ
    labels = "".join(nt["label"] for nt in win._model.scene()["nucleotides"])
    assert labels == DEMO_SEQ
    assert win._model.scene()["style"]["show_letters"] is True


# -- (3) render-type combo carries its full item text (cutoff fix) ----------


def test_render_combo_has_full_items(win):
    cb = win._panel._render_cb
    items = [cb.itemText(i) for i in range(cb.count())]
    assert items == ["none", "res_type", "paired"]
    # the dock is wide enough that the field can show the longest label
    assert win._dock.minimumWidth() >= 280


# -- (5) residue numbering --------------------------------------------------


def test_numbering_toggle_adds_numbers_at_expected_indices(win):
    win._numbers_cb.setChecked(True)
    win._number_interval_sp.setValue(10)
    scene = win._model.scene()
    assert "numbers" in scene
    nums = scene["numbers"]
    idxs = sorted(n["index"] for n in nums)
    n = len(DEMO_SS)  # 41
    expected = sorted({0, n - 1} | {i for i in range(n) if (i + 1) % 10 == 0})
    assert idxs == expected  # 1, 10, 20, 30, 40, and the last (41)
    texts = {n["index"]: n["text"] for n in nums}
    assert texts[0] == "1"
    assert texts[9] == "10"
    assert texts[n - 1] == str(n)


def test_numbering_off_emits_no_numbers(win):
    win._numbers_cb.setChecked(False)
    scene = win._model.scene()
    assert scene.get("numbers", []) == []


def test_numbering_interval_changes_count():
    model = EditorModel.from_ss(DEMO_SS)
    model.set_numbering(True, 5)
    every5 = {nn["index"] for nn in model.scene()["numbers"]}
    model.set_numbering(True, 10)
    every10 = {nn["index"] for nn in model.scene()["numbers"]}
    assert every10 < every5  # a wider interval -> fewer numbers


def test_numbers_absent_from_jupyter_style_payload():
    # a plain model (numbering off) never adds the optional numbers key
    model = EditorModel.from_ss(DEMO_SS)
    assert "numbers" not in model.scene()
