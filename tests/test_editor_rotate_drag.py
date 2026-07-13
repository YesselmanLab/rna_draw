"""Headless tests for the desktop editor's direct drag-to-rotate gesture.

The visible rotate handle (dashed ring + orange knob) is gone. In Rotate
("move") mode the user clicks a helix -- it highlights teal -- and then DRAGS
the highlighted helix (grabbing any of its nucleotides) to rotate it about its
junction pivot. These tests drive the view's real press -> move -> release
handlers (a real mouse cannot run headlessly) and assert:

  * a press-then-drag rotates the selected helix and commits with the honest
    never-silent flag, recording exactly ONE undo entry (pushed on the first
    move, not on the press);
  * a press-then-release with NO movement selects the helix but does not
    rotate it and records no undo entry;
  * selecting a helix adds no `RotateHandle` overlay (the class is gone).

All Qt runs under ``QT_QPA_PLATFORM=offscreen``.
"""

from __future__ import annotations

import math
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest

from rna_draw.gui.model import DEMO_SS

PySide6 = pytest.importorskip("PySide6")

from PySide6 import QtCore, QtGui  # noqa: E402

from rna_draw.gui.desktop.app import build_app  # noqa: E402

_SS = DEMO_SS  # three-way junction; helix [8,17] is the second stem
_GRAB_NT = 8  # a paired nucleotide inside that helix


@pytest.fixture
def win():
    _, w = build_app(ss=_SS, seq=None)
    yield w
    w.close()


def _mouse(view, kind, scene_pt, buttons=QtCore.Qt.MouseButton.LeftButton):
    """Dispatch a synthetic mouse event whose position maps to ``scene_pt``.

    Round-trips scene -> viewport (via ``mapFromScene``) so the handler's own
    ``mapToScene`` recovers ``scene_pt`` regardless of the offscreen viewport
    size. ``kind`` is one of "press" / "move" / "release".
    """
    vpt = QtCore.QPointF(view.mapFromScene(scene_pt))
    type_ = {
        "press": QtCore.QEvent.Type.MouseButtonPress,
        "move": QtCore.QEvent.Type.MouseMove,
        "release": QtCore.QEvent.Type.MouseButtonRelease,
    }[kind]
    btn = (
        QtCore.Qt.MouseButton.LeftButton
        if kind != "move"
        else QtCore.Qt.MouseButton.NoButton
    )
    ev = QtGui.QMouseEvent(
        type_, vpt, btn, buttons, QtCore.Qt.KeyboardModifier.NoModifier
    )
    getattr(view, f"mouse{kind.capitalize()}Event")(ev)


def _scene_pt(view, idx):
    return view._scene.nt_positions()[idx]


def _rotated_about(pt, pivot_scene, angle):
    dx, dy = pt.x() - pivot_scene.x(), pt.y() - pivot_scene.y()
    ca, sa = math.cos(angle), math.sin(angle)
    return QtCore.QPointF(
        pivot_scene.x() + dx * ca - dy * sa,
        pivot_scene.y() + dx * sa + dy * ca,
    )


def test_drag_rotates_helix_and_records_one_undo(win):
    view, model = win._view, win._model
    assert win._options_panel.options().mode == "move"
    assert model.can_undo() is False

    grab = _scene_pt(view, _GRAB_NT)
    before = list(zip(model._x, model._y))

    # Press on a helix nt: it selects the helix + arms the rotate drag.
    _mouse(view, "press", grab)
    assert model.sel_kind == "helix"
    assert view._pivot is not None
    assert view._rotating is True
    assert model.can_undo() is False  # the press itself pushes nothing

    # Drag it around the junction pivot -> the helix rotates.
    px, py = view._pivot
    pivot_scene = QtCore.QPointF(px, -py)  # engine y-up -> scene y-down
    _mouse(view, "move", _rotated_about(grab, pivot_scene, math.radians(20.0)))

    assert view._rotate_pushed is True
    assert model.can_undo() is True  # exactly one entry, pushed on first move
    after = list(zip(model._x, model._y))
    moved = [i for i in range(model._sel_start, model._sel_end + 1)
             if (model._x[i], model._y[i]) != before[i]]
    assert moved, "selected helix did not move on drag"

    # A further move must NOT add a second undo entry (one gesture = one entry).
    _mouse(view, "move", _rotated_about(grab, pivot_scene, math.radians(30.0)))
    win_model_undo = len(model._undo)
    assert win_model_undo == 1

    _mouse(view, "release", grab)
    assert view._rotating is False
    # The committed flag is whatever the native checker honestly found.
    assert isinstance(model.flagged, bool)
    assert after != before


def test_click_without_movement_selects_but_does_not_rotate(win):
    view, model = win._view, win._model
    grab = _scene_pt(view, _GRAB_NT)
    before = list(zip(model._x, model._y))

    _mouse(view, "press", grab)
    _mouse(view, "release", grab)  # no move in between

    assert model.sel_kind == "helix"  # the helix is highlighted/selected
    assert view._rotating is False
    assert model.can_undo() is False  # no rotation -> no undo entry
    assert list(zip(model._x, model._y)) == before  # geometry unchanged


def test_selection_adds_no_rotate_handle(win):
    view, model = win._view, win._model
    _mouse(view, "press", _scene_pt(view, _GRAB_NT))
    _mouse(view, "release", _scene_pt(view, _GRAB_NT))

    # The handle module/class is gone; nothing handle-like is armed.
    from rna_draw.gui.desktop import scene_view as sv

    assert not hasattr(sv, "RotateHandle")
    assert not hasattr(view, "_handle")
    with pytest.raises(ImportError):
        from rna_draw.gui.desktop.handles import RotateHandle  # noqa: F401
