"""The VISIBLE rotate handle -- the discoverability fix.

The prior Jupyter attempt failed because rotation was an invisible
click-then-drag with no on-screen affordance. Here the selected helix's
junction pivot grows an obvious knob on a dashed ring: the user can SEE where
to grab and that dragging it turns the helix. `RnaGraphicsView` hit-tests the
knob and drives `EditorModel.rotate_selection` from the drag.

All coordinates are SCENE coordinates (engine x, negated engine y); the view
owns the engine<->scene mapping.
"""

from __future__ import annotations

import math

from PySide6 import QtCore, QtGui, QtWidgets

TEAL = QtGui.QColor("#1fb6a6")
KNOB_FILL = QtGui.QColor("#ff8a3d")
KNOB_EDGE = QtGui.QColor("#ffffff")


class RotateHandle(QtWidgets.QGraphicsItemGroup):
    """A dashed ring around the pivot with a grabbable knob.

    Args:
        cx: Pivot x in scene coords.
        cy: Pivot y in scene coords.
        radius: Ring radius in scene units (roughly the selected helix's
            reach, so the ring hugs the helix it rotates).
        knob_r: Knob disk radius in scene units.
    """

    def __init__(self, cx: float, cy: float, radius: float, knob_r: float = 9.0) -> None:
        super().__init__()
        self._cx = cx
        self._cy = cy
        self._radius = radius
        self._knob_r = knob_r
        self._angle = math.radians(-45.0)  # display angle of the knob on the ring
        self.setZValue(1000)

        d = 2 * radius
        self._ring = QtWidgets.QGraphicsEllipseItem(cx - radius, cy - radius, d, d)
        ring_pen = QtGui.QPen(TEAL, 1.6, QtCore.Qt.PenStyle.DashLine)
        ring_pen.setCosmetic(True)
        self._ring.setPen(ring_pen)
        self._ring.setBrush(QtCore.Qt.BrushStyle.NoBrush)
        self.addToGroup(self._ring)

        self._pivot_dot = QtWidgets.QGraphicsEllipseItem(cx - 2.5, cy - 2.5, 5.0, 5.0)
        self._pivot_dot.setBrush(QtGui.QBrush(TEAL))
        self._pivot_dot.setPen(QtGui.QPen(QtCore.Qt.PenStyle.NoPen))
        self.addToGroup(self._pivot_dot)

        self._line = QtWidgets.QGraphicsLineItem()
        line_pen = QtGui.QPen(TEAL, 1.2, QtCore.Qt.PenStyle.DashLine)
        line_pen.setCosmetic(True)
        self._line.setPen(line_pen)
        self.addToGroup(self._line)

        self._knob = QtWidgets.QGraphicsEllipseItem()
        knob_pen = QtGui.QPen(KNOB_EDGE, 2.0)
        knob_pen.setCosmetic(True)
        self._knob.setPen(knob_pen)
        self._knob.setBrush(QtGui.QBrush(KNOB_FILL))
        self.addToGroup(self._knob)

        self._place_knob()

    def _place_knob(self) -> None:
        kx = self._cx + self._radius * math.cos(self._angle)
        ky = self._cy + self._radius * math.sin(self._angle)
        r = self._knob_r
        self._knob.setRect(kx - r, ky - r, 2 * r, 2 * r)
        self._line.setLine(self._cx, self._cy, kx, ky)

    def knob_scene_pos(self) -> QtCore.QPointF:
        """Current knob center in scene coordinates."""
        return QtCore.QPointF(
            self._cx + self._radius * math.cos(self._angle),
            self._cy + self._radius * math.sin(self._angle),
        )

    def knob_radius(self) -> float:
        """Knob hit radius in scene units."""
        return self._knob_r

    def center(self) -> QtCore.QPointF:
        """Pivot center in scene coordinates."""
        return QtCore.QPointF(self._cx, self._cy)

    def set_display_angle(self, angle: float) -> None:
        """Point the knob along `angle` (radians, scene space) for drag feedback."""
        self._angle = angle
        self._place_knob()

    def display_angle(self) -> float:
        """Current knob display angle (radians, scene space)."""
        return self._angle


__all__ = ["RotateHandle"]
