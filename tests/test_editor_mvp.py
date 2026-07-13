"""Headless MVP tests for the Jupyter RNA-editor anywidget (`rna_draw.gui`).

These exercise the KERNEL side end-to-end without a browser: scene shape,
click-to-select helix resolution, and the never-silent rotate/commit path
(the real `check_overlaps_native` arbiter drives `scene["flagged"]`). Pointer
events themselves cannot be driven headlessly -- that is called out in the
task's final report -- so the JS is only checked for existence + valid ESM
shape here.
"""

from __future__ import annotations

import pathlib

from rna_draw.gui import editor
from rna_draw.gui.editor_widget import RnaEditor
from rna_draw.layout.base import params_at_node_r
from rna_draw.layout.postpass import rotate_range
from rna_draw.overlap import OverlapParams
from rna_draw.overlap_native import check_overlaps_native

CLEAN_SS = "((((...((((....))))...((((....))))...))))"
PK_SS = "((((....[[[[....))))....]]]]"


def _scaled_params(w: RnaEditor) -> OverlapParams:
    return params_at_node_r(OverlapParams(node_r=10.0), w._node_r)


def test_scene_is_well_formed_dict():
    w = editor(CLEAN_SS)
    scene = w.scene
    assert isinstance(scene, dict)
    assert {"nucleotides", "pairs", "viewport", "node_r", "flagged", "overlaps", "pivot"} <= set(scene)
    assert len(scene["nucleotides"]) == len(CLEAN_SS)
    for i, nt in enumerate(scene["nucleotides"]):
        assert nt["id"] == i
        assert isinstance(nt["x"], float) and isinstance(nt["y"], float)
        assert nt["fill"].startswith("#")
    vp = scene["viewport"]
    assert {"min_x", "min_y", "w", "h"} <= set(vp)
    assert vp["w"] > 0 and vp["h"] > 0
    # every pair references valid nt indices and a known kind
    assert len(scene["pairs"]) > 0
    for p in scene["pairs"]:
        assert 0 <= p["i"] < p["j"] < len(CLEAN_SS)
        assert p["kind"] in ("nested", "crossing")


def test_pseudoknot_scene_well_formed():
    w = editor(PK_SS)
    scene = w.scene
    assert len(scene["nucleotides"]) == len(PK_SS)
    assert scene["viewport"]["w"] > 0
    # a PK cannot be drawn clean -> honestly flagged, never silently clean
    assert scene["flagged"] is True


def test_select_resolves_expected_helix_slice():
    w = editor(CLEAN_SS)
    # index 8 opens the second helix; its drawn partner is 17 -> slice [8, 17]
    w._handle_msg(w, {"type": "select", "index": 8}, [])
    assert w.selection == list(range(8, 18))
    assert w.status == "Helix 8-17 selected"
    assert w.scene["pivot"] is not None
    # the outer helix opener selects the whole structure slice
    w._handle_msg(w, {"type": "select", "index": 0}, [])
    assert w.selection == list(range(0, len(CLEAN_SS)))


def test_select_exterior_nt_selects_nothing():
    # a structure with 5' exterior tail: index 0 is inside no helix stem
    w = editor("..((((....))))..")
    w._handle_msg(w, {"type": "select", "index": 0}, [])
    assert w.selection == []
    assert "not inside a helix" in w.status


def test_rotate_matches_rotate_range_and_updates_scene():
    w = editor(CLEAN_SS)
    w._handle_msg(w, {"type": "select", "index": 8}, [])
    start, end = w._sel_start, w._sel_end
    cx, cy = w._pivot
    x0, y0 = list(w._x), list(w._y)
    scene_before = w.scene

    exp_x, exp_y = rotate_range(x0, y0, start, end, cx, cy, 0.3)
    w._handle_msg(w, {"type": "rotate", "angle": 0.3}, [])

    # (a) coords are exactly rotate_range's output
    assert w._x == exp_x
    assert w._y == exp_y
    # (b) scene dict updated (new object, coords changed inside the slice)
    assert w.scene is not scene_before
    assert w.scene["nucleotides"][8]["x"] == exp_x[8]


def test_rotate_into_overlap_flags_true():
    w = editor(CLEAN_SS)
    w._handle_msg(w, {"type": "select", "index": 8}, [])
    exp_x, exp_y = rotate_range(w._x, w._y, w._sel_start, w._sel_end, *w._pivot, 0.3)
    expected_count = check_overlaps_native(exp_x, exp_y, w._pair_map, _scaled_params(w))
    assert expected_count > 0  # sanity: this rotation really does collide

    w._handle_msg(w, {"type": "rotate", "angle": 0.3}, [])
    assert w.scene["flagged"] is True
    assert len(w.scene["overlaps"]) > 0
    assert f"{expected_count} overlap" in w.status


def test_small_rotate_stays_clean_flag_false():
    w = editor(CLEAN_SS)
    w._handle_msg(w, {"type": "select", "index": 8}, [])
    exp_x, exp_y = rotate_range(w._x, w._y, w._sel_start, w._sel_end, *w._pivot, 0.05)
    assert check_overlaps_native(exp_x, exp_y, w._pair_map, _scaled_params(w)) == 0

    w._handle_msg(w, {"type": "rotate", "angle": 0.05}, [])
    assert w.scene["flagged"] is False
    assert w.scene["overlaps"] == []
    assert "clean" in w.status


def test_editor_js_exists_and_is_valid_esm():
    import rna_draw.gui.editor_widget as ew

    js = pathlib.Path(ew.__file__).parent / "static" / "editor.js"
    assert js.is_file()
    text = js.read_text()
    assert len(text) > 0
    assert "export function render" in text
