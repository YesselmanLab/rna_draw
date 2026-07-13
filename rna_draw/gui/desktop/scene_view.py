"""`RnaGraphicsView` / `RnaScene`: the QGraphicsView over `EditorModel`.

Builds one `QGraphicsEllipseItem` per nucleotide (storing its index), a
backbone `QGraphicsPathItem`, pair `QGraphicsLineItem`s and routed PK-A
polylines from the model's portable scene dict. Engine space is y-up, so
scene coordinates are ``(x, -y)``; the view owns that flip.

Interaction (the discoverable joint edit):
    * left-click a nucleotide -> `EditorModel.select` -> the helix slice
      highlights teal and a VISIBLE `RotateHandle` appears at the junction;
    * drag that handle -> live `EditorModel.rotate_selection`, re-checked
      in-process every move (overlaps tint red);
    * left-click empty space -> deselect;
    * wheel -> zoom; middle/right drag -> pan.
"""

from __future__ import annotations

import math

from PySide6 import QtCore, QtGui, QtWidgets

from rna_draw.gui.model import EditorModel

from .handles import RotateHandle
from .options_panel import EditorOptions

# Palette (works acceptably on light or dark canvas).
BACKBONE = QtGui.QColor("#8a8f98")
PAIR = QtGui.QColor("#5b6068")
CROSSING = QtGui.QColor("#c0673a")
ROUTED = QtGui.QColor("#3a7bd5")
NT_EDGE = QtGui.QColor("#2b2f36")
SELECT_TEAL = QtGui.QColor("#1fb6a6")
OVERLAP_RED = QtGui.QColor("#e5484d")
LABEL = QtGui.QColor("#101418")
NUMBER = QtGui.QColor("#6b7280")

_NT_INDEX_KEY = 0

_PEN_STYLES = {
    "solid": QtCore.Qt.PenStyle.SolidLine,
    "dash": QtCore.Qt.PenStyle.DashLine,
    "dot": QtCore.Qt.PenStyle.DotLine,
}


def _color(value: str, fallback: QtGui.QColor) -> QtGui.QColor:
    """Parse a ``#rrggbb`` string to a QColor, falling back on anything invalid."""
    c = QtGui.QColor(value)
    return c if c.isValid() else fallback


class RnaScene(QtWidgets.QGraphicsScene):
    """Renders a model's scene dict; caches nt screen positions for hit tests."""

    def __init__(self) -> None:
        super().__init__()
        self._nt_positions: list[QtCore.QPointF] = []
        self._node_r: float = 10.0

    def nt_positions(self) -> list[QtCore.QPointF]:
        """Per-nucleotide center in scene coords, index-aligned."""
        return self._nt_positions

    def node_r(self) -> float:
        """Disk radius in scene units."""
        return self._node_r

    def build(
        self,
        scene: dict,
        selection: tuple[int, int] | None,
        show_overlaps: bool = True,
    ) -> None:
        """Rebuild all items from a portable scene dict.

        Args:
            scene: `EditorModel.scene()` output (engine space).
            selection: The selected helix slice ``(start, end)`` or ``None``.
            show_overlaps: When False, overlapping nucleotides are NOT tinted
                red (they keep their normal fill). The overlap set is still
                carried in the scene dict and reported honestly elsewhere --
                this flag is purely about the visual red tint. Defaults True
                so any non-desktop caller renders exactly as before.
        """
        self.clear()
        self._nt_positions = []
        nts = scene["nucleotides"]
        self._node_r = float(scene.get("node_r", 10.0))
        overlaps = set(scene.get("overlaps", [])) if show_overlaps else set()
        sel_lo, sel_hi = (selection if selection is not None else (None, None))
        # The highlighted set (any granularity): an explicit index set from the
        # scene dict (a motif ring is non-contiguous), falling back to the
        # contiguous helix range so a plain scene dict renders as before.
        selected_set = set(scene.get("selected", []))
        if not selected_set and sel_lo is not None:
            selected_set = set(range(sel_lo, sel_hi + 1))

        # Resolved visual style (colors/widths/letters). Falls back to the
        # module constants so a scene dict without a `style` block (e.g. the
        # Jupyter payload) renders exactly as before.
        style = scene.get("style") or {}
        backbone_pen = QtGui.QPen(
            _color(style.get("backbone_color", ""), BACKBONE),
            float(style.get("backbone_width", 2.0)),
        )
        pair_col = _color(style.get("pair_color", ""), PAIR)
        pair_w = float(style.get("pair_width", 1.6))
        crossing_col = _color(style.get("crossing_color", ""), CROSSING)
        crossing_w = float(style.get("crossing_width", 1.6))
        routed_col = _color(style.get("routed_color", ""), ROUTED)
        routed_w = float(style.get("routed_width", 1.6))
        routed_style = _PEN_STYLES.get(
            style.get("routed_style", "dash"), QtCore.Qt.PenStyle.DashLine
        )
        edge_col = _color(style.get("nt_edge_color", ""), NT_EDGE)
        edge_w = float(style.get("nt_edge_width", 1.0))
        default_fill = style.get("default_fill", "#9aa0a6")
        label_col = _color(style.get("letter_color", ""), LABEL)
        show_letters = bool(style.get("show_letters", True))
        letter_scale = float(style.get("letter_scale", 1.0))

        pos = [QtCore.QPointF(float(n["x"]), -float(n["y"])) for n in nts]
        self._nt_positions = pos

        # Backbone path (consecutive nts).
        if len(pos) >= 2:
            path = QtGui.QPainterPath(pos[0])
            for p in pos[1:]:
                path.lineTo(p)
            bb = self.addPath(path, backbone_pen)
            bb.setZValue(-10)

        # Base pairs.
        for pair in scene.get("pairs", []):
            i, j = pair["i"], pair["j"]
            if i >= len(pos) or j >= len(pos):
                continue
            crossing = pair.get("kind") == "crossing"
            col = crossing_col if crossing else pair_col
            wid = crossing_w if crossing else pair_w
            line = self.addLine(QtCore.QLineF(pos[i], pos[j]), QtGui.QPen(col, wid))
            line.setZValue(-5)

        # Routed PK-A polylines.
        for routed in scene.get("routed_lines", []):
            pts = routed.get("points", [])
            if len(pts) < 2:
                continue
            rpath = QtGui.QPainterPath(QtCore.QPointF(pts[0][0], -pts[0][1]))
            for px, py in pts[1:]:
                rpath.lineTo(QtCore.QPointF(px, -py))
            ritem = self.addPath(rpath, QtGui.QPen(routed_col, routed_w, routed_style))
            ritem.setZValue(-4)

        # Nucleotide disks.
        r = self._node_r
        for n in nts:
            idx = n["id"]
            center = pos[idx]
            selected = idx in selected_set
            if idx in overlaps:
                fill = OVERLAP_RED
            elif selected:
                fill = SELECT_TEAL
            else:
                fill = QtGui.QColor(n.get("fill", default_fill))
            item = QtWidgets.QGraphicsEllipseItem(center.x() - r, center.y() - r, 2 * r, 2 * r)
            edge = SELECT_TEAL if selected else edge_col
            item.setPen(QtGui.QPen(edge, 2.0 if selected else edge_w))
            item.setBrush(QtGui.QBrush(fill))
            item.setData(_NT_INDEX_KEY, idx)
            item.setZValue(0)
            self.addItem(item)
            label = n.get("label", "")
            if label and show_letters:
                text = QtWidgets.QGraphicsSimpleTextItem(str(label))
                text.setBrush(QtGui.QBrush(label_col))
                font = text.font()
                font.setPointSizeF(max(4.0, r * 0.9 * letter_scale))
                text.setFont(font)
                br = text.boundingRect()
                text.setPos(center.x() - br.width() / 2, center.y() - br.height() / 2)
                text.setZValue(1)
                self.addItem(text)

        # Residue-position numbers (VARNA-style). Optional and non-interactive:
        # small text just outside the disks; they carry no nt index so the
        # hit-test (which uses `_nt_positions`) is untouched.
        number_col = _color(style.get("number_color", ""), NUMBER)
        for num in scene.get("numbers", []):
            nx, ny = float(num["x"]), -float(num["y"])
            ntext = QtWidgets.QGraphicsSimpleTextItem(str(num.get("text", "")))
            ntext.setBrush(QtGui.QBrush(number_col))
            nfont = ntext.font()
            nfont.setPointSizeF(max(3.0, r * 0.7))
            ntext.setFont(nfont)
            nbr = ntext.boundingRect()
            ntext.setPos(nx - nbr.width() / 2, ny - nbr.height() / 2)
            ntext.setZValue(2)
            self.addItem(ntext)


class RnaGraphicsView(QtWidgets.QGraphicsView):
    """Interactive view: click-to-select, drag-the-handle-to-rotate.

    Signals:
        selectionChanged: emits ``(start, end)`` tuple or ``None``.
        overlapChanged: emits ``(flagged: bool, count: int)`` after a commit.
    """

    selectionChanged = QtCore.Signal(object)
    overlapChanged = QtCore.Signal(bool, int)

    def __init__(self, parent=None) -> None:
        self._scene = RnaScene()
        super().__init__(self._scene, parent)
        self._model: EditorModel | None = None
        self._options = EditorOptions()
        self._handle: RotateHandle | None = None
        self._pivot: tuple[float, float] | None = None  # engine-space selection pivot
        self._rotating = False
        self._last_angle = 0.0
        # Residue free-form drag (Select mode, residue granularity): the nt
        # being dragged, or None when no residue drag is in progress.
        self._drag_residue: int | None = None
        self._panning = False
        self._pan_start = QtCore.QPoint()
        self.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        self.setDragMode(QtWidgets.QGraphicsView.DragMode.NoDrag)
        self.setTransformationAnchor(QtWidgets.QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setMouseTracking(True)

    # -- model wiring -------------------------------------------------------

    def set_options(self, options: EditorOptions) -> None:
        """Adopt the editor behavior options and re-render at the new setting.

        Pure visual: overlap tinting is re-evaluated from the flags, but the
        geometry and the never-silent overlap COMPUTATION are untouched.
        """
        self._options = options
        self._rebuild()

    def _show_overlaps(self) -> bool:
        """Whether the red overlap tint is allowed to draw right now.

        Gated by the "Highlight overlaps" master switch AND, when
        "Highlight only after drag" is on, suppressed during an active drag.
        This is VISUAL ONLY -- the model always computes overlaps and the
        status bar always reports them (never-silent).
        """
        if not self._options.highlight_overlaps:
            return False
        if self._rotating and self._options.highlight_only_after_drag:
            return False
        return True

    def set_model(self, model: EditorModel) -> None:
        """Attach a model and render its initial scene."""
        self._model = model
        self._clear_handle()
        self._apply_background(model.scene().get("style", {}))
        self._rebuild()
        self.fit()
        self._emit_overlap()

    def restyle(self) -> None:
        """Re-render from the model after a visual-only style change.

        Rebuilds every item with the model's freshly restyled scene dict but
        keeps the current zoom/pan and selection (no `fit`). Also applies the
        canvas background from the scene's style block.
        """
        if self._model is None:
            return
        self._apply_background(self._model.scene().get("style", {}))
        self._rebuild()

    def _apply_background(self, style: dict) -> None:
        bg = QtGui.QColor(style.get("background", "#ffffff"))
        if bg.isValid():
            self.setBackgroundBrush(QtGui.QBrush(bg))

    def _rebuild(self) -> None:
        if self._model is None:
            return
        selection = self._model.selection
        # `scene.build()` calls `clear()`, which DELETES the handle's C++
        # object. Capture its display angle, drop the stale wrapper, rebuild,
        # then recreate a fresh handle at the (fixed) junction pivot. Re-adding
        # the deleted wrapper would raise "Internal C++ object already deleted"
        # on the very first drag move (mouseMoveEvent -> _rebuild).
        angle = self._handle.display_angle() if self._handle is not None else None
        self._handle = None
        self._scene.build(self._model.scene(), selection, self._show_overlaps())
        if selection is not None and self._pivot is not None:
            self._show_handle(selection[0], selection[1], self._pivot, angle)

    def fit(self) -> None:
        """Fit the whole structure in the viewport."""
        rect = self._scene.itemsBoundingRect()
        if rect.isNull():
            return
        self.fitInView(rect.adjusted(-20, -20, 20, 20), QtCore.Qt.AspectRatioMode.KeepAspectRatio)

    # -- selection / handle -------------------------------------------------

    def _clear_handle(self) -> None:
        if self._handle is not None:
            if self._handle.scene() is not None:
                self._scene.removeItem(self._handle)
            self._handle = None
        self._pivot = None

    def _show_handle(
        self, start: int, end: int, pivot: tuple[float, float], angle: float | None = None
    ) -> None:
        """Place the visible rotate handle at the selection's junction pivot."""
        self._clear_handle()
        self._pivot = pivot  # engine space; fixed while rotating about the junction
        cx, cy = pivot[0], -pivot[1]  # engine -> scene
        positions = self._scene.nt_positions()
        reach = 0.0
        for i in range(start, end + 1):
            if i < len(positions):
                reach = max(reach, math.hypot(positions[i].x() - cx, positions[i].y() - cy))
        radius = max(reach + 1.5 * self._scene.node_r(), 3.0 * self._scene.node_r())
        self._handle = RotateHandle(cx, cy, radius, knob_r=max(7.0, self._scene.node_r() * 0.9))
        if angle is not None:
            self._handle.set_display_angle(angle)
        self._scene.addItem(self._handle)

    def _select_at(self, idx: int) -> None:
        """Resolve a click at ``idx`` per the active mode + granularity.

        Move mode always selects a helix (rotate handle appears). Select mode
        selects at the current granularity: a helix gets the rotate handle, a
        residue becomes free-form draggable, a motif is highlight-only.
        """
        self._drag_residue = None
        if self._model is None:
            return
        if self._options.mode != "select":
            # Move mode: classic helix select + rotate handle.
            resolved = self._model.select(idx)
            if resolved is None:
                self._clear_handle()
                self._rebuild()
                self.selectionChanged.emit(None)
            else:
                start, end, pivot = resolved
                self._rebuild()
                self._show_handle(start, end, pivot)
                self.selectionChanged.emit((start, end))
            return
        # Select mode: granular selection.
        sel = self._model.select_at(idx, self._options.granularity)
        if sel is None:
            self._clear_handle()
            self._rebuild()
            self.selectionChanged.emit(None)
            return
        if sel.kind == "helix":
            self._rebuild()
            self._show_handle(sel.start, sel.end, sel.pivot)
        else:
            self._clear_handle()
            if sel.kind == "residue" and sel.indices:
                self._drag_residue = sel.indices[0]
            self._rebuild()
        self.selectionChanged.emit(
            (sel.indices[0], sel.indices[-1]) if sel.indices else None
        )

    # -- mouse interaction --------------------------------------------------

    def mousePressEvent(self, event) -> None:
        button = event.button()
        if button in (QtCore.Qt.MouseButton.MiddleButton, QtCore.Qt.MouseButton.RightButton):
            self._panning = True
            self._pan_start = event.position().toPoint()
            self.setCursor(QtCore.Qt.CursorShape.ClosedHandCursor)
            event.accept()
            return
        if button == QtCore.Qt.MouseButton.LeftButton and self._model is not None:
            spos = self.mapToScene(event.position().toPoint())
            # 1) grab the handle knob?
            if self._handle is not None:
                knob = self._handle.knob_scene_pos()
                grab = self._handle.knob_radius() * 1.8
                if math.hypot(spos.x() - knob.x(), spos.y() - knob.y()) <= grab:
                    self._rotating = True
                    self._last_angle = self._angle_from_pivot(spos)
                    self.setCursor(QtCore.Qt.CursorShape.ClosedHandCursor)
                    event.accept()
                    return
            # 2) click a nucleotide?
            idx = self._nearest_nt(spos)
            if idx is not None:
                self._select_at(idx)
                event.accept()
                return
            # 3) empty space -> deselect
            self._model.deselect()
            self._clear_handle()
            self._drag_residue = None
            self._rebuild()
            self.selectionChanged.emit(None)
            event.accept()
            return
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event) -> None:
        if self._panning:
            delta = event.position().toPoint() - self._pan_start
            self._pan_start = event.position().toPoint()
            hbar = self.horizontalScrollBar()
            vbar = self.verticalScrollBar()
            hbar.setValue(hbar.value() - delta.x())
            vbar.setValue(vbar.value() - delta.y())
            event.accept()
            return
        if (
            self._drag_residue is not None
            and self._model is not None
            and event.buttons() & QtCore.Qt.MouseButton.LeftButton
        ):
            spos = self.mapToScene(event.position().toPoint())
            # scene -> engine: engine y is up, scene y is down.
            result = self._model.move_residue(self._drag_residue, spos.x(), -spos.y())
            self._rebuild()
            self.overlapChanged.emit(result.flagged, len(result.overlaps))
            event.accept()
            return
        if self._rotating and self._model is not None and self._handle is not None:
            spos = self.mapToScene(event.position().toPoint())
            cur = self._angle_from_pivot(spos)
            # scene y is flipped vs engine y, so an engine-CCW rotation is a
            # scene-CW one: negate the scene-space delta to keep the visual
            # drag direction matching the helix's motion.
            delta_scene = cur - self._last_angle
            self._last_angle = cur
            result = self._model.rotate_selection(-delta_scene)
            self._handle.set_display_angle(cur)
            self._rebuild()
            self.overlapChanged.emit(result.flagged, len(result.overlaps))
            event.accept()
            return
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event) -> None:
        if self._panning:
            self._panning = False
            self.setCursor(QtCore.Qt.CursorShape.ArrowCursor)
            event.accept()
            return
        if self._rotating:
            self._rotating = False
            self.setCursor(QtCore.Qt.CursorShape.ArrowCursor)
            # The drag has settled: rebuild so "highlight only after drag" can
            # now reveal any red tint that was suppressed mid-drag.
            self._rebuild()
            event.accept()
            return
        if self._drag_residue is not None and event.button() == QtCore.Qt.MouseButton.LeftButton:
            # Residue drag settled; keep it selected for a subsequent drag.
            self._rebuild()
            event.accept()
            return
        super().mouseReleaseEvent(event)

    def wheelEvent(self, event) -> None:
        factor = 1.15 if event.angleDelta().y() > 0 else 1.0 / 1.15
        self.scale(factor, factor)
        event.accept()

    # -- helpers ------------------------------------------------------------

    def _angle_from_pivot(self, spos: QtCore.QPointF) -> float:
        center = self._handle.center()
        return math.atan2(spos.y() - center.y(), spos.x() - center.x())

    def _nearest_nt(self, spos: QtCore.QPointF) -> int | None:
        positions = self._scene.nt_positions()
        if not positions:
            return None
        threshold = self._scene.node_r() * 1.6
        best_i, best_d = None, threshold
        for i, p in enumerate(positions):
            d = math.hypot(spos.x() - p.x(), spos.y() - p.y())
            if d <= best_d:
                best_d, best_i = d, i
        return best_i

    def _emit_overlap(self) -> None:
        if self._model is not None:
            self.overlapChanged.emit(self._model.flagged, len(self._model.scene()["overlaps"]))


__all__ = ["RnaGraphicsView", "RnaScene"]
