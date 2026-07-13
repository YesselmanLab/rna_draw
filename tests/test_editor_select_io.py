"""Tests for the editor's multi-format Open, whole-helix motif selection, the
generic set-selection API, and the box/marquee select tool.

Real mouse drags cannot be driven headlessly, so the box selection is exercised
through the view's `_finish_box` commit step (a scene-space rectangle), exactly
what the rubber-band release calls. All Qt runs under ``offscreen``.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest

from rna_draw.gui.model import DEMO_SS, EditorModel
from rna_draw.io_formats import pairs_to_dot_bracket, parse_structure_file

_SCRATCH = os.environ.get(
    "RNA_DRAW_SCRATCH",
    "/private/tmp/claude-503/-Users-jyesselman2-local-code-python-developing-rna-draw/"
    "641187e9-926c-4fa2-a1cd-2ce8215f508a/scratchpad",
)
os.makedirs(_SCRATCH, exist_ok=True)

_SS = DEMO_SS


# -- (a) multi-format parsing -----------------------------------------------


def test_parse_dbn_bprna_style(tmp_path):
    path = tmp_path / "s.dbn"
    path.write_text(
        "#Name: test-hairpin\n#Length: 12\n"
        "GGGGAAAACCCC\n"
        "((((....))))\n"
    )
    ss, seq, name = parse_structure_file(str(path))
    assert ss == "((((....))))"
    assert seq == "GGGGAAAACCCC"
    assert name == "test-hairpin"


def test_parse_lone_dot_bracket(tmp_path):
    path = tmp_path / "s.dot"
    path.write_text("((((....))))\n")
    ss, seq, name = parse_structure_file(str(path))
    assert ss == "((((....))))"
    assert seq is None


def test_parse_dbn_preserves_pseudoknot(tmp_path):
    path = tmp_path / "pk.dbn"
    path.write_text("GGGCAAAUCCCAAAGGG\n(((.[[[.))).]]]..\n")
    ss, _seq, _name = parse_structure_file(str(path))
    assert "[" in ss and "]" in ss


def test_parse_fasta_with_structure(tmp_path):
    path = tmp_path / "s.fasta"
    path.write_text(">my rna\nGGGG\nAAAA\nCCCC\n((((....))))\n")
    ss, seq, name = parse_structure_file(str(path))
    assert ss == "((((....))))"
    assert seq == "GGGGAAAACCCC"
    assert name == "my rna"


def test_parse_ct(tmp_path):
    # 6-nt hairpin: 1-2 pair 6-5, 3-4 unpaired loop.
    path = tmp_path / "s.ct"
    path.write_text(
        "6  ENERGY = -1.0  hairpin\n"
        "1 G 0 2 6 1\n"
        "2 G 1 3 5 2\n"
        "3 A 2 4 0 3\n"
        "4 A 3 5 0 4\n"
        "5 C 4 6 2 5\n"
        "6 C 5 0 1 6\n"
    )
    ss, seq, name = parse_structure_file(str(path))
    assert seq == "GGAACC"
    assert ss == "((..))"
    assert name.startswith("ENERGY")


def test_parse_bpseq(tmp_path):
    path = tmp_path / "s.bpseq"
    path.write_text(
        "1 G 6\n2 G 5\n3 A 0\n4 A 0\n5 C 2\n6 C 1\n"
    )
    ss, seq, _name = parse_structure_file(str(path))
    assert seq == "GGAACC"
    assert ss == "((..))"


def test_parse_bad_file_raises(tmp_path):
    path = tmp_path / "junk.txt"
    path.write_text("this is not a structure\nnor is this line\n")
    with pytest.raises(ValueError):
        parse_structure_file(str(path))


def test_pairs_to_dot_bracket_pseudoknot():
    # crossing pairs (0,4) and (2,6) -> two bracket families.
    ss = pairs_to_dot_bracket(8, [(0, 4), (2, 6)])
    assert ss[0] == "(" and ss[4] == ")"
    assert ss[2] == "[" and ss[6] == "]"


# -- (b) whole-helix motif selection ----------------------------------------


def test_motif_on_stem_returns_full_helix():
    model = EditorModel.from_ss(_SS)
    # nt 8 sits on the 4-rung stem (7,18)..(10,15).
    sel = model.select_at(8, "motif")
    assert sel.kind == "motif"
    assert set(sel.indices) == {7, 8, 9, 10, 15, 16, 17, 18}
    # both strands, and every stacked rung's partner is present
    for i in (7, 8, 9, 10):
        assert model._pair_map[i] in sel.indices


def test_motif_on_loop_returns_loop_members():
    model = EditorModel.from_ss(_SS)
    sel = model.select_at(5, "motif")  # nt 5 is an unpaired loop nt
    assert sel.kind == "motif"
    assert 5 in sel.indices
    # a loop ring is non-contiguous (spans past the stem it encloses)
    assert sel.indices != list(range(sel.indices[0], sel.indices[-1] + 1))


# -- (c) generic set-selection API ------------------------------------------


def test_select_indices_sets_selection_and_scene():
    model = EditorModel.from_ss(_SS)
    chosen = {2, 3, 4, 20, 21}
    sel = model.select_indices(chosen)
    assert sel.kind == "range"
    assert set(sel.indices) == chosen
    assert model.sel_kind == "range"
    assert set(model.sel_indices) == chosen
    scene = model.scene()
    assert set(scene["selected"]) == chosen
    assert scene["selection_kind"] == "range"


def test_select_indices_drops_out_of_range():
    model = EditorModel.from_ss(_SS)
    sel = model.select_indices({0, 1, 999, -3})
    assert set(sel.indices) == {0, 1}


# -- (d) box / marquee select via the view ----------------------------------


@pytest.fixture
def view():
    from PySide6 import QtWidgets

    from rna_draw.gui.desktop.scene_view import RnaGraphicsView

    if QtWidgets.QApplication.instance() is None:
        QtWidgets.QApplication([])
    v = RnaGraphicsView()
    v.set_model(EditorModel.from_ss(_SS))
    yield v


def test_box_select_selects_nts_inside(view):
    from PySide6 import QtCore

    positions = view._scene.nt_positions()
    # Build a rectangle tightly around the first four nucleotides' centers.
    xs = [positions[i].x() for i in range(4)]
    ys = [positions[i].y() for i in range(4)]
    pad = view._scene.node_r()
    rect = QtCore.QRectF(
        QtCore.QPointF(min(xs) - pad, min(ys) - pad),
        QtCore.QPointF(max(xs) + pad, max(ys) + pad),
    )
    view._finish_box(rect)
    inside = set(view._model.sel_indices)
    assert {0, 1, 2, 3} <= inside
    # and every reported index really lies in the rectangle
    for i in inside:
        assert rect.contains(positions[i])
    assert view._model.sel_kind == "range"


def test_box_select_additive_unions(view):
    from PySide6 import QtCore

    positions = view._scene.nt_positions()

    def rect_around(idxs):
        xs = [positions[i].x() for i in idxs]
        ys = [positions[i].y() for i in idxs]
        pad = view._scene.node_r()
        return QtCore.QRectF(
            QtCore.QPointF(min(xs) - pad, min(ys) - pad),
            QtCore.QPointF(max(xs) + pad, max(ys) + pad),
        )

    view._finish_box(rect_around([0, 1]))
    first = set(view._model.sel_indices)
    view._finish_box(rect_around([2, 3]), additive=True)
    combined = set(view._model.sel_indices)
    assert first <= combined
    assert {2, 3} <= combined


def test_box_select_render_png(view):
    from PySide6 import QtCore, QtGui

    positions = view._scene.nt_positions()
    # A wide box spanning the first stem so several teal disks show.
    idxs = list(range(0, 11))
    xs = [positions[i].x() for i in idxs]
    ys = [positions[i].y() for i in idxs]
    pad = view._scene.node_r() * 1.5
    rect = QtCore.QRectF(
        QtCore.QPointF(min(xs) - pad, min(ys) - pad),
        QtCore.QPointF(max(xs) + pad, max(ys) + pad),
    )
    view._finish_box(rect)
    assert len(view._model.sel_indices) >= 4

    img = QtGui.QImage(700, 700, QtGui.QImage.Format.Format_ARGB32)
    img.fill(QtGui.QColor("white"))
    painter = QtGui.QPainter(img)
    view._scene.render(painter)
    painter.end()
    out = os.path.join(_SCRATCH, "box_select.png")
    assert img.save(out)
    assert os.path.getsize(out) > 0


# -- (e) MainWindow Open dialog filter + parse path --------------------------


def test_open_filter_and_dbn_load(tmp_path):
    from rna_draw.gui.desktop.app import build_app
    from rna_draw.io_formats import OPEN_FILTER

    for ext in ("dbn", "dot", "ct", "bpseq", "fasta", "fa", "txt"):
        assert f"*.{ext}" in OPEN_FILTER

    _, win = build_app(ss=DEMO_SS, seq=None)
    try:
        path = tmp_path / "load.dbn"
        path.write_text("#Name: x\nGGGGAAAACCCC\n((((....))))\n")
        ss, seq, _name = parse_structure_file(str(path))
        win.load_ss(ss, seq)
        assert win._model is not None
        assert len(win._model.scene()["nucleotides"]) == 12
        assert win._ss_edit.text() == "((((....))))"
        assert win._seq_edit.text() == "GGGGAAAACCCC"
    finally:
        win.close()
