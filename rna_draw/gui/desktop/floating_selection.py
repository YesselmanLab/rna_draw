"""`FloatingSelectionPanel`: a contextual popover for styling a selection.

A small, frameless, rounded card that floats OVER the canvas (it is a child of
the `RnaGraphicsView` viewport) right next to whatever the user just selected.
It offers the same per-region styling the left "Selection styling" section does
-- recolor the selection, add a highlight halo, clear -- but inline, so the
user never has to travel to the side panel for a quick tweak.

The panel is intentionally dumb: it renders a readout for the current selection
and emits a signal per action carrying the chosen color. `MainWindow` wires
those signals to the SAME model calls + `push_undo` the side panel uses, so the
floating and side-panel edits are identical and undoable. It never touches
geometry (purely visual styling), so the never-silent contract is untouched.
"""

from __future__ import annotations

from PySide6 import QtCore, QtWidgets

from .style_panel import ColorButton
from .theme import THEMES

# Gap (px) kept between the selection's bounding box and the panel, and the
# minimum inset kept from the viewport edges when clamping.
_GAP = 12
_EDGE = 6

# Default swatch colors, mirrored from the left "Selection styling" section so
# the two entry points open on the same colors.
_DEFAULT_COLOR = "#e5484d"
_DEFAULT_HIGHLIGHT = "#ffd300"


class FloatingSelectionPanel(QtWidgets.QFrame):
    """A floating, themed popover with inline per-selection styling actions.

    Signals:
        colorRequested: emits the chosen ``#rrggbb`` to recolor the selection.
        highlightRequested: emits the chosen ``#rrggbb`` to highlight it.
        clearRequested: emits when the user clears the selection's styling.
        closed: emits when the user dismisses the panel with the "x".
    """

    colorRequested = QtCore.Signal(str)
    highlightRequested = QtCore.Signal(str)
    clearRequested = QtCore.Signal()
    closed = QtCore.Signal()

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setObjectName("floatingSelection")
        self.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        # Never pull keyboard focus away from the canvas; show without stealing
        # activation from the main window.
        self.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_ShowWithoutActivating, True)
        self._build_ui()
        self.apply_theme("dark")
        self.hide()

    # -- construction -------------------------------------------------------

    def _build_ui(self) -> None:
        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(10, 8, 10, 10)
        outer.setSpacing(7)

        # Header: selection readout + a close "x".
        header = QtWidgets.QHBoxLayout()
        header.setContentsMargins(0, 0, 0, 0)
        header.setSpacing(8)
        self._readout = QtWidgets.QLabel("No selection")
        self._readout.setObjectName("floatReadout")
        header.addWidget(self._readout, 1)
        close_btn = QtWidgets.QToolButton()
        close_btn.setObjectName("floatClose")
        close_btn.setText("×")  # times sign
        close_btn.setToolTip("Dismiss this panel")
        close_btn.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        close_btn.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)
        close_btn.clicked.connect(self.closed.emit)
        header.addWidget(close_btn, 0)
        outer.addLayout(header)

        # Color row: swatch + action.
        color_row = QtWidgets.QHBoxLayout()
        color_row.setContentsMargins(0, 0, 0, 0)
        color_row.setSpacing(6)
        self._color_btn = ColorButton(_DEFAULT_COLOR)
        self._color_btn.setToolTip("Color to paint the selected nucleotides")
        self._color_btn.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)
        color_act = QtWidgets.QPushButton("Color")
        color_act.setToolTip("Recolor the selection (overrides the scheme)")
        color_act.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)
        color_act.clicked.connect(lambda: self.colorRequested.emit(self._color_btn.color()))
        color_row.addWidget(self._color_btn, 0)
        color_row.addWidget(color_act, 1)
        outer.addLayout(color_row)

        # Highlight row: swatch + action.
        hl_row = QtWidgets.QHBoxLayout()
        hl_row.setContentsMargins(0, 0, 0, 0)
        hl_row.setSpacing(6)
        self._highlight_btn = ColorButton(_DEFAULT_HIGHLIGHT)
        self._highlight_btn.setToolTip("Halo color for the highlighted region")
        self._highlight_btn.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)
        hl_act = QtWidgets.QPushButton("Highlight")
        hl_act.setToolTip("Add a translucent halo over the selection")
        hl_act.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)
        hl_act.clicked.connect(
            lambda: self.highlightRequested.emit(self._highlight_btn.color())
        )
        hl_row.addWidget(self._highlight_btn, 0)
        hl_row.addWidget(hl_act, 1)
        outer.addLayout(hl_row)

        # Clear action (color + highlight on the selection).
        clear_act = QtWidgets.QPushButton("Clear")
        clear_act.setToolTip("Remove color + highlight styling")
        clear_act.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)
        clear_act.clicked.connect(self.clearRequested.emit)
        outer.addWidget(clear_act)

    # -- selection readout --------------------------------------------------

    def set_selection(self, kind: str | None, indices) -> None:
        """Update the header readout from the current selection kind + indices."""
        self._readout.setText(self._readout_text(kind, list(indices)))

    @staticmethod
    def _readout_text(kind: str | None, indices: list[int]) -> str:
        """Compact human-readable label for a selection (1-based positions)."""
        n = len(indices)
        if n == 0:
            return "No selection"
        if kind == "helix":
            lo, hi = min(indices) + 1, max(indices) + 1
            return f"Helix {lo}–{hi} · {n} nt"
        if kind == "residue" and n == 1:
            return f"Residue {indices[0] + 1}"
        label = {"motif": "Motif", "range": "Selection"}.get(kind or "", "Selection")
        return f"{label} · {n} nt selected"

    # -- positioning --------------------------------------------------------

    def place_near(self, anchor: QtCore.QRect, viewport: QtCore.QRect) -> None:
        """Position the panel just outside `anchor`, clamped inside `viewport`.

        Prefers to sit ABOVE the selection's bounding box (`anchor`, in the
        viewport's pixel coordinates); if there is no room above, drops it
        below. Then clamps both axes so the whole panel stays within
        `viewport`. Pure geometry -- safe to call headlessly.
        """
        self.adjustSize()
        w, h = self.width(), self.height()

        x = anchor.center().x() - w // 2
        y = anchor.top() - h - _GAP
        if y < viewport.top() + _EDGE:
            below = anchor.bottom() + _GAP
            # Use below only if it fits better; otherwise clamp above anyway.
            if below + h <= viewport.bottom() - _EDGE or below >= y:
                y = below

        x = max(viewport.left() + _EDGE, min(x, viewport.right() - w - _EDGE))
        y = max(viewport.top() + _EDGE, min(y, viewport.bottom() - h - _EDGE))
        self.move(int(x), int(y))

    # -- theming ------------------------------------------------------------

    def apply_theme(self, name: str) -> None:
        """Restyle the popover with the active theme tokens (dark/light)."""
        t = THEMES.get(name, THEMES["dark"])
        self.setStyleSheet(_panel_qss(t))


def _panel_qss(t: dict) -> str:
    """QSS for the popover from a theme token dict (see `theme.THEMES`)."""
    return f"""
QFrame#floatingSelection {{
    background: {t['panel']};
    border: 1px solid {t['line']};
    border-radius: 10px;
}}
QFrame#floatingSelection QLabel {{ color: {t['ink']}; background: transparent; }}
QLabel#floatReadout {{ color: {t['ink']}; font-weight: 600; }}
QFrame#floatingSelection QPushButton {{
    padding: 4px 12px;
    border-radius: 6px;
    background: {t['field']};
    color: {t['ink']};
    border: 1px solid {t['line']};
}}
QFrame#floatingSelection QPushButton:hover {{
    background: {t['hover']};
    border: 1px solid {t['accent']};
}}
QToolButton#floatClose {{
    border: none;
    background: transparent;
    color: {t['muted']};
    font-size: 16px;
    font-weight: 600;
    padding: 0 3px;
}}
QToolButton#floatClose:hover {{ color: {t['ink']}; }}
"""


__all__ = ["FloatingSelectionPanel"]
