"""Tests for the desktop MVP: the frontend-agnostic `EditorModel` plus a
headless construction smoke of the Qt window under ``offscreen``.

The real mouse drag cannot be driven headlessly, so all editing LOGIC is
tested through `EditorModel` (pure, no Qt); the Qt layer is only smoke-tested
for construction (`QApplication` + `MainWindow` + scene build without error).
"""

from __future__ import annotations

import os

import pytest

from rna_draw.gui.model import DEMO_SS, EditorModel, Result
from rna_draw.layout.postpass import rotate_range

# The three-way junction used across these tests; helix [8,17] is the second
# stem, its junction pivot is the central multiloop center.
_SS = DEMO_SS


def test_model_scene_well_formed():
    model = EditorModel.from_ss(_SS)
    scene = model.scene()
    assert len(scene["nucleotides"]) == len(_SS)
    for i, nt in enumerate(scene["nucleotides"]):
        assert nt["id"] == i
        assert isinstance(nt["x"], float) and isinstance(nt["y"], float)
        assert nt["r"] == scene["node_r"]
    # a clean demo layout is not flagged and has no overlaps
    assert model.flagged is False
    assert scene["overlaps"] == []
    # pairs present (this structure has three hairpins + a closing stem)
    assert scene["pairs"]
    assert model.engine_name


def test_select_returns_expected_helix_and_pivot():
    model = EditorModel.from_ss(_SS)
    resolved = model.select(8)
    assert resolved is not None
    start, end, pivot = resolved
    assert (start, end) == (8, 17)
    assert model.selection == (8, 17)
    assert len(pivot) == 2
    # pivot is the parent-loop center, distinct from the clicked nt
    assert model.pivot == pivot


def test_select_unpaired_exterior_returns_none():
    # a purely exterior/unpaired structure -> no enclosing helix
    model = EditorModel.from_ss(".....")
    assert model.select(2) is None
    assert model.selection is None


def test_rotate_small_is_clean_and_moves_slice():
    model = EditorModel.from_ss(_SS)
    start, end, pivot = model.select(8)
    before = [(model._x[i], model._y[i]) for i in range(len(model._x))]
    result = model.rotate_selection(0.05)
    assert isinstance(result, Result)
    assert result.flagged is False
    assert result.overlaps == []
    # the selected slice actually moved
    moved = [i for i in range(start, end + 1) if (model._x[i], model._y[i]) != before[i]]
    assert moved, "rotation should displace the selected helix"
    # nucleotides OUTSIDE the slice are untouched (rigid, local move)
    for i in range(len(before)):
        if not (start <= i <= end):
            assert (model._x[i], model._y[i]) == before[i]


def test_rotate_matches_direct_rotate_range():
    model = EditorModel.from_ss(_SS)
    start, end, pivot = model.select(8)
    nx, ny = rotate_range(model._x, model._y, start, end, pivot[0], pivot[1], 0.05)
    model.rotate_selection(0.05)
    assert model._x == pytest.approx(nx)
    assert model._y == pytest.approx(ny)


def test_rotate_large_flags_overlaps():
    model = EditorModel.from_ss(_SS)
    model.select(8)
    result = model.rotate_selection(0.6)
    assert result.flagged is True
    assert result.overlaps  # offending indices for red tinting
    # the scene reflects the flag (never silently clean)
    scene = model.scene()
    assert scene["flagged"] is True
    assert set(scene["overlaps"]) == set(result.overlaps)


def test_translate_moves_slice_and_rechecks():
    model = EditorModel.from_ss(_SS)
    start, end, _ = model.select(8)
    before_x = list(model._x)
    result = model.translate_selection(3.0, -2.0)
    assert isinstance(result, Result)
    for i in range(start, end + 1):
        assert model._x[i] == pytest.approx(before_x[i] + 3.0)


def test_rotate_without_selection_is_noop():
    model = EditorModel.from_ss(_SS)
    result = model.rotate_selection(0.5)
    assert result.flagged is False


def test_qt_window_constructs_headless():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from rna_draw.gui.desktop.app import build_app

    app, window = build_app(argv=[])
    assert app is not None
    window.show()
    # the central view built a scene from the demo structure
    view = window.centralWidget()
    scene = view.scene()
    assert len(scene.items()) > 0
    # selecting a helix through the view stores its junction pivot (no handle)
    view.set_model(window._model)
    view._select_at(8)
    assert view._pivot is not None
    resolved = window._model.select(8)
    assert resolved is not None
    window.close()
