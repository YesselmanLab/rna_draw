"""Headless tests for per-region styling: color overrides + highlights.

Covers `EditorModel.apply_color_to_selection` / `clear_color_on_selection`,
`highlight_selection` / `clear_highlights`, the override's HIGHEST fill
precedence (it beats an active render-type scheme), the scene `highlights`
payload, save/load round-trips through the Document `extra` band, and the
MainWindow "Selection styling" panel section applying live. The pure-model
tests need no Qt; the window tests run under ``QT_QPA_PLATFORM=offscreen``.
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


def _fills(model: EditorModel) -> dict[int, str]:
    return {nt["id"]: nt["fill"] for nt in model.scene()["nucleotides"]}


# -- (a) recolor a selection; clearing restores the scheme fill -------------


def test_apply_color_override_recolors_only_selection():
    model = EditorModel.from_ss(DEMO_SS, seq=DEMO_SEQ)
    baseline = _fills(model)
    picked = [0, 1, 2, 3]
    model.select_indices(picked)
    model.apply_color_to_selection("#ff0000")
    # deselect so the transient teal tint does not mask the fill readout
    model.deselect()
    fills = _fills(model)
    for i in picked:
        assert fills[i] == "#ff0000"
    # every OTHER nucleotide keeps its original scheme fill
    for i, col in baseline.items():
        if i not in picked:
            assert fills[i] == col


def test_clear_color_on_selection_restores_scheme():
    model = EditorModel.from_ss(DEMO_SS, seq=DEMO_SEQ)
    baseline = _fills(model)
    model.select_indices([5, 6, 7])
    model.apply_color_to_selection("#00ff00")
    # re-select the same set and clear its overrides
    model.select_indices([5, 6, 7])
    model.clear_color_on_selection()
    model.deselect()
    assert _fills(model) == baseline


# -- (b) highlight a region; scene emits it; clear removes -------------------


def test_highlight_selection_emits_and_clears():
    model = EditorModel.from_ss(DEMO_SS, seq=DEMO_SEQ)
    assert "highlights" not in model.scene()  # absent when none defined
    model.select_indices([8, 9, 10, 11])
    model.highlight_selection("#ffff00")
    scene = model.scene()
    assert scene["highlights"] == [{"indices": [8, 9, 10, 11], "color": "#ffff00"}]
    model.clear_highlights()
    assert "highlights" not in model.scene()


def test_highlight_no_selection_is_noop():
    model = EditorModel.from_ss(DEMO_SS, seq=DEMO_SEQ)
    model.deselect()
    model.highlight_selection("#ffff00")
    assert model.highlights == []


# -- (c) overrides + highlights ROUND-TRIP through documents + save/load -----


def test_region_styling_round_trips_through_document():
    model = EditorModel.from_ss(DEMO_SS, seq=DEMO_SEQ)
    model.select_indices([0, 1, 2, 3])
    model.apply_color_to_selection("#ff0000")
    model.select_indices([20, 21, 22])
    model.highlight_selection("#00ffff")

    doc = model.to_document()
    assert doc.extra["color_overrides"] == {"0": "#ff0000", "1": "#ff0000",
                                            "2": "#ff0000", "3": "#ff0000"}
    assert doc.extra["highlights"] == [{"indices": [20, 21, 22], "color": "#00ffff"}]

    restored = EditorModel.from_document(doc)
    assert restored.color_overrides == {0: "#ff0000", 1: "#ff0000", 2: "#ff0000", 3: "#ff0000"}
    assert restored.highlights == [{"indices": [20, 21, 22], "color": "#00ffff"}]
    # and the restored scene reflects them
    rfills = {nt["id"]: nt["fill"] for nt in restored.scene()["nucleotides"]}
    assert rfills[0] == "#ff0000"
    assert restored.scene()["highlights"] == [{"indices": [20, 21, 22], "color": "#00ffff"}]


def test_region_styling_round_trips_through_file(tmp_path):
    model = EditorModel.from_ss(DEMO_SS, seq=DEMO_SEQ)
    model.select_indices([4, 5])
    model.apply_color_to_selection("#abcdef")
    model.select_indices([30, 31, 32])
    model.highlight_selection("#123456")

    path = tmp_path / "styled.rnadoc.json"
    model.to_document().save(path)
    restored = EditorModel.from_document(Document.load(path))
    assert restored.color_overrides == {4: "#abcdef", 5: "#abcdef"}
    assert restored.highlights == [{"indices": [30, 31, 32], "color": "#123456"}]


# -- (d) override precedence beats the render_type scheme --------------------


def test_override_beats_render_type_scheme():
    model = EditorModel.from_ss(DEMO_SS, seq=DEMO_SEQ)
    model.preset.extra.setdefault("view", {})["render_type"] = "res_type"
    model.restyle()
    scheme = _fills(model)
    # index 0 is a G -> res_type paints it #ff6666 (see test_style_panel)
    assert scheme[0] == "#ff6666"
    model.select_indices([0])
    model.apply_color_to_selection("#000000")
    model.deselect()
    assert _fills(model)[0] == "#000000"  # override wins over the scheme


# -- (e) MainWindow section present + applying via the window updates scene ---


def test_window_has_selection_styling_section():
    _, win = build_app(ss=DEMO_SS, seq=DEMO_SEQ, argv=[])
    try:
        section = win._panel.section("Selection styling")
        assert isinstance(section, CollapsibleSection)
    finally:
        win.close()


def test_window_apply_color_and_highlight_updates_scene():
    _, win = build_app(ss=DEMO_SS, seq=DEMO_SEQ, argv=[])
    try:
        win._model.select_indices([0, 1, 2, 3])
        win._sel_color_btn.set_color("#ff0000")
        win._on_color_selection()
        win._model.deselect()
        fills = {nt["id"]: nt["fill"] for nt in win._model.scene()["nucleotides"]}
        assert fills[0] == "#ff0000"

        win._model.select_indices([8, 9, 10, 11])
        win._sel_highlight_btn.set_color("#ffd300")
        win._on_highlight_selection()
        assert win._model.scene()["highlights"] == [
            {"indices": [8, 9, 10, 11], "color": "#ffd300"}
        ]
    finally:
        win.close()


def test_selection_styling_enabled_only_with_selection():
    _, win = build_app(ss=DEMO_SS, seq=DEMO_SEQ, argv=[])
    try:
        win._model.deselect()
        win._update_selection_styling()
        assert all(not w.isEnabled() for w in win._sel_apply_widgets)
        assert win._sel_readout.text() == "No selection"
        win._model.select_indices([0, 1, 2])
        win._update_selection_styling()
        assert all(w.isEnabled() for w in win._sel_apply_widgets)
        assert "3 nt selected" in win._sel_readout.text()
    finally:
        win.close()
