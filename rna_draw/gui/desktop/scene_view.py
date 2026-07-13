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
from .theme import SELECT_TEAL as SELECT_TEAL_HEX
from .theme import canvas_colors

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
        backbone_col = _color(style.get("backbone_color", ""), BACKBONE)
        backbone_w = float(style.get("backbone_width", 2.0))
        backbone_pen = QtGui.QPen(backbone_col, backbone_w)
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
        # Live display scale on the drawn disk radius (NOT the layout node_r the
        # hit-test / handle / halos key off) and whether to draw disks at all.
        show_spheres = bool(style.get("show_spheres", True))
        sphere_scale = float(style.get("sphere_scale", 1.0))

        pos = [QtCore.QPointF(float(n["x"]), -float(n["y"])) for n in nts]
        self._nt_positions = pos
        # Overall render style ("spheres" = classic disks; "letters" = the
        # RFview-style colored-letters mode). A plain string so more styles can
        # be added later without touching the hit-test / interaction paths
        # (which key off the cached nt CENTERS below, not the drawn items).
        render_style = str(style.get("render_style", "spheres"))

        # Backbone + base pairs for the SPHERES style (the letters style draws
        # its own gapped backbone + typed pair symbols in `_build_letters`).
        if render_style != "letters":
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

        # Persistent highlight halos: a larger, semi-transparent disk behind
        # each highlighted nucleotide (zvalue BELOW the opaque disks) so the
        # region visibly "glows" as a ring around the letters without hiding
        # them. Non-interactive (carries no nt index).
        r = self._node_r
        halo_r = r * 1.9
        for hl in scene.get("highlights", []):
            col = _color(hl.get("color", "#ffff00"), QtGui.QColor("#ffff00"))
            col.setAlpha(90)
            brush = QtGui.QBrush(col)
            for idx in hl.get("indices", []):
                if idx >= len(pos):
                    continue
                c = pos[idx]
                halo = QtWidgets.QGraphicsEllipseItem(
                    c.x() - halo_r, c.y() - halo_r, 2 * halo_r, 2 * halo_r
                )
                halo.setBrush(brush)
                halo.setPen(QtGui.QPen(QtCore.Qt.PenStyle.NoPen))
                halo.setZValue(-2)
                self.addItem(halo)

        if render_style == "letters":
            # RFview-style: colored letters, a gapped backbone, and typed pair
            # symbols. No disks. Draws its own backbone + pairs (skipped above).
            self._build_letters(
                nts,
                pos,
                scene,
                selected_set,
                overlaps,
                default_fill,
                letter_scale,
                backbone_col,
                backbone_w,
                pair_col,
                pair_w,
                crossing_col,
                crossing_w,
            )
        else:
            # Nucleotide disks + letters. The DRAWN disk radius is the layout
            # radius scaled by the live sphere display scale; when "Show spheres"
            # is off the disk items are skipped entirely (letters, numbers,
            # backbone and pairs still render). Click-select keys off cached nt
            # centers, not disk items, so interaction is unaffected by the drawn
            # size or visibility.
            draw_r = r * sphere_scale
            for n in nts:
                idx = n["id"]
                center = pos[idx]
                selected = idx in selected_set
                if show_spheres:
                    if idx in overlaps:
                        fill = OVERLAP_RED
                    elif selected:
                        fill = SELECT_TEAL
                    else:
                        fill = QtGui.QColor(n.get("fill", default_fill))
                    item = QtWidgets.QGraphicsEllipseItem(
                        center.x() - draw_r, center.y() - draw_r, 2 * draw_r, 2 * draw_r
                    )
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

    # -- "Letters" render style (RFview-like) -------------------------------

    def _build_letters(
        self,
        nts: list[dict],
        pos: list[QtCore.QPointF],
        scene: dict,
        selected_set: set[int],
        overlaps: set[int],
        default_fill: str,
        letter_scale: float,
        backbone_col: QtGui.QColor,
        backbone_w: float,
        pair_col: QtGui.QColor,
        pair_w: float,
        crossing_col: QtGui.QColor,
        crossing_w: float,
    ) -> None:
        """Render the scene as RFview-style COLORED LETTERS (no disks).

        Each nucleotide is its letter, colored by its resolved ``fill`` (or teal
        when selected / red when overlapping so selection + overlap feedback
        works without disks); a blank-sequence position becomes a small dot. The
        backbone is a thin line drawn only BETWEEN the letters (shortened at both
        ends so it stops clear of each glyph), and each base pair is a TYPE
        symbol keyed off its ``bond`` (=, -, wobble, LW placeholder). Letters sit
        on top (zvalue 1); backbone/pairs are drawn below.
        """
        r = self._node_r
        # 1) Letters first so each glyph's radius can be measured for the gaps.
        #    `glyph_r[i]` is a half-extent of nt i's drawn glyph (or dot).
        glyph_r = [r * 0.5] * len(pos)
        font_pt = max(4.0, r * 0.95 * letter_scale)
        for n in nts:
            idx = n["id"]
            center = pos[idx]
            if idx in overlaps:
                col = OVERLAP_RED
            elif idx in selected_set:
                col = SELECT_TEAL
            else:
                col = QtGui.QColor(n.get("fill", default_fill))
            label = n.get("label", "")
            if label and str(label).strip():
                text = QtWidgets.QGraphicsSimpleTextItem(str(label))
                text.setBrush(QtGui.QBrush(col))
                font = text.font()
                font.setPointSizeF(font_pt)
                font.setBold(True)
                text.setFont(font)
                br = text.boundingRect()
                text.setPos(center.x() - br.width() / 2, center.y() - br.height() / 2)
                text.setZValue(1)
                self.addItem(text)
                glyph_r[idx] = 0.5 * math.hypot(br.width(), br.height())
            else:
                # No letter (blank sequence): a small dot marks the position.
                dot_r = r * 0.28
                dot = QtWidgets.QGraphicsEllipseItem(
                    center.x() - dot_r, center.y() - dot_r, 2 * dot_r, 2 * dot_r
                )
                dot.setBrush(QtGui.QBrush(col))
                dot.setPen(QtGui.QPen(QtCore.Qt.PenStyle.NoPen))
                dot.setZValue(1)
                self.addItem(dot)
                glyph_r[idx] = dot_r * 1.4

        # Extra breathing room so a connector never touches the glyph box.
        margin = max(1.5, r * 0.28)

        # 2) Gapped backbone: a thin line between CONSECUTIVE letters, stopping
        #    short of each glyph at both ends ("line only between letters").
        bb_pen = QtGui.QPen(backbone_col, backbone_w)
        for a in range(len(pos) - 1):
            self._gapped_line(
                pos[a], pos[a + 1], glyph_r[a] + margin, glyph_r[a + 1] + margin,
                bb_pen, -10.0,
            )

        # 3) Typed pair symbols between paired letters.
        for pair in scene.get("pairs", []):
            i, j = pair["i"], pair["j"]
            if i >= len(pos) or j >= len(pos):
                continue
            crossing = pair.get("kind") == "crossing"
            col = crossing_col if crossing else pair_col
            wid = crossing_w if crossing else pair_w
            self._pair_symbol(
                pos[i], pos[j], glyph_r[i] + margin, glyph_r[j] + margin,
                pair.get("bond", "other"), col, wid,
            )

    def _gapped_line(
        self,
        p_a: QtCore.QPointF,
        p_b: QtCore.QPointF,
        gap_a: float,
        gap_b: float,
        pen: QtGui.QPen,
        z: float,
    ) -> tuple[QtCore.QPointF, QtCore.QPointF] | None:
        """Draw a line from ``p_a`` to ``p_b`` shortened by a gap at each end.

        Returns the shortened endpoints, or ``None`` when the endpoints are too
        close for anything to remain after the gaps (nothing drawn).
        """
        dx, dy = p_b.x() - p_a.x(), p_b.y() - p_a.y()
        length = math.hypot(dx, dy)
        if length <= gap_a + gap_b:
            return None
        ux, uy = dx / length, dy / length
        s = QtCore.QPointF(p_a.x() + ux * gap_a, p_a.y() + uy * gap_a)
        e = QtCore.QPointF(p_b.x() - ux * gap_b, p_b.y() - uy * gap_b)
        line = self.addLine(QtCore.QLineF(s, e), pen)
        line.setZValue(z)
        return s, e

    def _pair_symbol(
        self,
        p_i: QtCore.QPointF,
        p_j: QtCore.QPointF,
        gap_i: float,
        gap_j: float,
        bond: str,
        col: QtGui.QColor,
        wid: float,
    ) -> None:
        """Draw the pair rung as a TYPE symbol keyed off ``bond``.

        ``gc`` -> a double line (``=``, two parallel offset lines); ``au`` -> a
        single line (``-``); ``gu`` -> a single line + an OPEN circle at the
        midpoint (wobble); anything else -> a single line + a filled circle
        marker (a Leontis-Westhof cis-WC-style placeholder). The rung is
        shortened at both ends so it stops clear of the letters. A small helper
        so more bond styles are trivial to add. Drawn below the letters (z=-5).
        """
        dx, dy = p_j.x() - p_i.x(), p_j.y() - p_i.y()
        length = math.hypot(dx, dy)
        if length <= gap_i + gap_j:
            return
        ux, uy = dx / length, dy / length
        s = QtCore.QPointF(p_i.x() + ux * gap_i, p_i.y() + uy * gap_i)
        e = QtCore.QPointF(p_j.x() - ux * gap_j, p_j.y() - uy * gap_j)
        px, py = -uy, ux  # unit perpendicular to the rung
        r = self._node_r
        pen = QtGui.QPen(col, wid)
        z = -5.0
        if bond == "gc":
            # Double line: two short parallel lines offset perpendicular.
            off = max(1.2, r * 0.16)
            for sgn in (-1.0, 1.0):
                ox, oy = px * off * sgn, py * off * sgn
                line = self.addLine(
                    QtCore.QLineF(s.x() + ox, s.y() + oy, e.x() + ox, e.y() + oy), pen
                )
                line.setZValue(z)
            return
        # Single line for au / gu / other.
        line = self.addLine(QtCore.QLineF(s, e), pen)
        line.setZValue(z)
        if bond in ("gu", "other"):
            mx, my = (s.x() + e.x()) / 2.0, (s.y() + e.y()) / 2.0
            cr = max(1.6, r * 0.22)
            circ = QtWidgets.QGraphicsEllipseItem(mx - cr, my - cr, 2 * cr, 2 * cr)
            circ.setPen(pen)
            # gu wobble = open circle; other (LW placeholder) = filled marker.
            circ.setBrush(
                QtGui.QBrush(QtCore.Qt.BrushStyle.NoBrush)
                if bond == "gu"
                else QtGui.QBrush(QtGui.QColor(col))
            )
            circ.setZValue(z + 0.5)
            self.addItem(circ)


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
        # Undo coalescing: a whole rotate/residue drag must record exactly ONE
        # undo entry. These flags mark whether `push_undo` has already fired for
        # the in-progress gesture, so only the FIRST mutating move pushes (the
        # snapshot is thus the pre-drag pose) and per-frame moves do not.
        self._rotate_pushed = False
        self._residue_pushed = False
        # Residue free-form drag (Select mode, residue granularity): the nt
        # being dragged, or None when no residue drag is in progress.
        self._drag_residue: int | None = None
        self._panning = False
        self._pan_start = QtCore.QPoint()
        # Box / marquee selection (Select mode, drag on empty space): a live,
        # clearly-visible teal rectangle drawn straight into the scene (a
        # QGraphicsRectItem styles reliably across platforms where a
        # QRubberBand does not), its scene-space origin, and whether the drag
        # adds to (Shift) rather than replaces the current selection.
        self._sel_rect: QtWidgets.QGraphicsRectItem | None = None
        self._box_origin_scene = QtCore.QPointF()
        self._box_additive = False
        # Themed canvas base + dot-grid colors (pushed in by `set_theme`).
        self._canvas_color = QtGui.QColor(canvas_colors("dark")[0])
        self._grid_color = QtGui.QColor(canvas_colors("dark")[1])
        # A user-picked canvas background (from the style panel) overrides the
        # theme canvas color; None means "use the theme".
        self._user_bg: QtGui.QColor | None = None
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

    def set_theme(self, name: str) -> None:
        """Adopt the active theme's canvas base + dot-grid colors and repaint."""
        canvas, grid = canvas_colors(name)
        self._canvas_color = QtGui.QColor(canvas)
        self._grid_color = QtGui.QColor(grid)
        if self.viewport() is not None:
            self.viewport().update()

    def _apply_background(self, style: dict) -> None:
        """Record any explicit user canvas color; the theme paints the default.

        A default white (`#ffffff`) background means "no override" -- the
        canvas then follows the active theme (near-black in dark). Anything
        else is a deliberate user choice and is honored as the canvas base.
        """
        bg = QtGui.QColor(style.get("background", "#ffffff"))
        self._user_bg = bg if (bg.isValid() and bg.name().lower() != "#ffffff") else None
        if self.viewport() is not None:
            self.viewport().update()

    def drawBackground(self, painter: QtGui.QPainter, rect: QtCore.QRectF) -> None:
        """Paint the themed canvas + a subtle dot-grid behind the scene.

        Base fill is the user's canvas color if set, else the theme canvas
        color; small dots on an ~18-unit lattice give the mockup's texture.
        Purely decorative -- nucleotide/pair items are drawn on top unchanged.
        """
        base = self._user_bg or self._canvas_color
        painter.fillRect(rect, base)
        spacing = 18.0
        r = 0.9
        painter.setPen(QtCore.Qt.PenStyle.NoPen)
        painter.setBrush(QtGui.QBrush(self._grid_color))
        left = spacing * math.floor(rect.left() / spacing)
        top = spacing * math.floor(rect.top() / spacing)
        y = top
        while y <= rect.bottom():
            x = left
            while x <= rect.right():
                painter.drawEllipse(QtCore.QPointF(x, y), r, r)
                x += spacing
            y += spacing

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
                self._residue_pushed = False
            self._rebuild()
        self.selectionChanged.emit(
            (sel.indices[0], sel.indices[-1]) if sel.indices else None
        )

    def _box_begin(self, view_pt: QtCore.QPoint, additive: bool = False) -> None:
        """Start a visible marquee at `view_pt` (viewport coords).

        Draws a teal-tinted, teal-bordered `QGraphicsRectItem` into the scene
        so the selection box is clearly visible while dragging (the user asked
        to "visually see the selection box"). Factored out so a headless test
        can drive begin/update/finish without a real mouse.
        """
        self._box_origin_scene = self.mapToScene(view_pt)
        self._box_additive = bool(additive)
        teal = QtGui.QColor(SELECT_TEAL_HEX)
        pen = QtGui.QPen(teal, 1.5, QtCore.Qt.PenStyle.DashLine)
        pen.setCosmetic(True)  # constant 1.5px border at any zoom
        fill = QtGui.QColor(teal)
        fill.setAlpha(48)  # translucent teal wash
        rect_item = QtWidgets.QGraphicsRectItem(
            QtCore.QRectF(self._box_origin_scene, self._box_origin_scene)
        )
        rect_item.setPen(pen)
        rect_item.setBrush(QtGui.QBrush(fill))
        rect_item.setZValue(900)  # above nts, below the rotate handle (1000)
        self._scene.addItem(rect_item)
        self._sel_rect = rect_item

    def _box_update(self, view_pt: QtCore.QPoint) -> None:
        """Grow the marquee to the current pointer position (viewport coords)."""
        if self._sel_rect is None:
            return
        cur = self.mapToScene(view_pt)
        self._sel_rect.setRect(QtCore.QRectF(self._box_origin_scene, cur).normalized())

    def _box_finish(self) -> None:
        """Commit the marquee: select enclosed nts, then remove the box item."""
        if self._sel_rect is None:
            return
        scene_rect = self._sel_rect.rect().normalized()
        additive = self._box_additive
        self._scene.removeItem(self._sel_rect)
        self._sel_rect = None
        self._finish_box(scene_rect, additive)

    def _finish_box(self, scene_rect: QtCore.QRectF, additive: bool = False) -> None:
        """Select every nucleotide whose center lies inside ``scene_rect``.

        The rubber-band tool's commit step, factored out so it can be driven
        directly (a real mouse drag cannot run headlessly). Hit-tests the
        cached nt centers against the SCENE-space rectangle and sets the
        model's selection to that index set (adding to the current selection
        when ``additive``). Highlight-only: no rotate handle is armed.

        Args:
            scene_rect: The selection rectangle in scene coordinates.
            additive: When True, union with the current selection.
        """
        if self._model is None:
            return
        positions = self._scene.nt_positions()
        hit = {i for i, p in enumerate(positions) if scene_rect.contains(p)}
        if additive:
            hit |= set(self._model.sel_indices)
        self._clear_handle()
        self._drag_residue = None
        self._model.select_indices(hit)
        self._rebuild()
        self.selectionChanged.emit((min(hit), max(hit)) if hit else None)

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
                    self._rotate_pushed = False
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
            # 3) empty space: Select mode -> start a box/marquee; else deselect
            if self._options.mode == "select":
                additive = bool(
                    event.modifiers() & QtCore.Qt.KeyboardModifier.ShiftModifier
                )
                self._box_begin(event.position().toPoint(), additive)
                event.accept()
                return
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
        if self._sel_rect is not None and event.buttons() & QtCore.Qt.MouseButton.LeftButton:
            self._box_update(event.position().toPoint())
            event.accept()
            return
        if (
            self._drag_residue is not None
            and self._model is not None
            and event.buttons() & QtCore.Qt.MouseButton.LeftButton
        ):
            spos = self.mapToScene(event.position().toPoint())
            # Coalesce the whole residue drag into ONE undo entry: push the
            # pre-drag pose on the first mutating move only.
            if not self._residue_pushed:
                self._model.push_undo()
                self._residue_pushed = True
            # scene -> engine: engine y is up, scene y is down.
            result = self._model.move_residue(self._drag_residue, spos.x(), -spos.y())
            self._rebuild()
            self.overlapChanged.emit(result.flagged, len(result.overlaps))
            event.accept()
            return
        if self._rotating and self._model is not None and self._handle is not None:
            spos = self.mapToScene(event.position().toPoint())
            # Coalesce the whole rotate drag into ONE undo entry: push the
            # pre-rotate pose on the first mutating move only.
            if not self._rotate_pushed:
                self._model.push_undo()
                self._rotate_pushed = True
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
        if self._sel_rect is not None:
            self._box_finish()
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
            # Residue drag settled; keep it selected for a subsequent drag, but
            # arm a fresh undo entry for that next drag.
            self._residue_pushed = False
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

    def selection_view_rect(self) -> QtCore.QRect | None:
        """Bounding rect of the current selection's disks in VIEWPORT pixels.

        Unions the cached scene-space centers of every selected nucleotide
        (padded by the disk radius), maps that scene rectangle into the
        viewport's pixel coordinates, and returns it -- the anchor a floating
        panel positions against. Returns ``None`` when nothing is selected or
        no positions are cached yet.
        """
        if self._model is None:
            return None
        indices = self._model.sel_indices
        if not indices:
            return None
        positions = self._scene.nt_positions()
        pts = [positions[i] for i in indices if 0 <= i < len(positions)]
        if not pts:
            return None
        r = self._scene.node_r()
        xs = [p.x() for p in pts]
        ys = [p.y() for p in pts]
        scene_rect = QtCore.QRectF(
            min(xs) - r, min(ys) - r, (max(xs) - min(xs)) + 2 * r, (max(ys) - min(ys)) + 2 * r
        )
        return self.mapFromScene(scene_rect).boundingRect()


__all__ = ["RnaGraphicsView", "RnaScene"]
