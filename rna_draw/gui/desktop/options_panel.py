"""`OptionsPanel`: the right-dock editor BEHAVIOR settings (not visual style).

Distinct from `StylePanel` (which edits how the drawing LOOKS): this panel
edits how the editor BEHAVES while you interact with it. The toggles are pure
app state the view reads on every rebuild -- they never touch geometry and
never suppress the never-silent overlap COMPUTATION, only whether an overlap
is drawn in red.

The panel carries one live `EditorOptions` dataclass and emits
`optionsChanged` whenever a toggle flips; `MainWindow` wires that to
`RnaGraphicsView.set_options` + a refresh.
"""

from __future__ import annotations

from dataclasses import dataclass

from PySide6 import QtCore, QtWidgets


@dataclass
class EditorOptions:
    """App-behavior flags the view reads (no geometry, no visual style).

    Args:
        highlight_overlaps: Master switch for the red overlap tinting.
            Default OFF -- overlaps are computed and reported honestly in the
            status bar but not drawn red unless the user opts in.
        highlight_only_after_drag: When highlighting is enabled, still
            suppress the red tint DURING an active drag (show it only once the
            drag settles). Default ON.
        snap_back_on_overlap: Reserved -- revert a move that overlaps.
            Default OFF (placeholder; not yet wired into the move logic).
        mode: Interaction mode governing how clicks/drags are read.
            ``"move"`` (default) = arrange: click selects a helix, its rotate
            handle appears, drag rotates. ``"select"`` = inspect: click
            selects at ``granularity``; a residue can be dragged; no rotate
            handle unless a helix is selected. (``"edit"`` reserved.)
        granularity: What a click selects in ``"select"`` mode:
            ``"residue"`` | ``"helix"`` (default) | ``"motif"``. Ignored in
            ``"move"`` mode (which always selects a helix).
    """

    highlight_overlaps: bool = False
    highlight_only_after_drag: bool = True
    snap_back_on_overlap: bool = False
    mode: str = "move"
    granularity: str = "helix"


class OptionsPanel(QtWidgets.QWidget):
    """Grouped behavior toggles over one live `EditorOptions`.

    Signals:
        optionsChanged: a toggle flipped; the view should re-read options and
            refresh (a purely visual rebuild -- no relayout, no re-check).
    """

    optionsChanged = QtCore.Signal()

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._options = EditorOptions()
        self._build_ui()

    # -- public API ---------------------------------------------------------

    def options(self) -> EditorOptions:
        """The live `EditorOptions` (mutated in place by the checkboxes)."""
        return self._options

    # -- construction -------------------------------------------------------

    def _build_ui(self) -> None:
        # No outer group-box: this panel is embedded in the left dock's
        # collapsible "Editing" section, whose header already titles it.
        form = QtWidgets.QVBoxLayout(self)
        form.setContentsMargins(6, 4, 6, 6)
        form.setSpacing(6)

        self._highlight_cb = QtWidgets.QCheckBox("Highlight overlaps")
        self._highlight_cb.setToolTip(
            "Tint overlapping nucleotides red. The overlap is ALWAYS checked and "
            "reported in the status bar regardless of this switch."
        )
        self._highlight_cb.setChecked(self._options.highlight_overlaps)
        self._highlight_cb.toggled.connect(self._on_highlight)
        form.addWidget(self._highlight_cb)

        self._after_drag_cb = QtWidgets.QCheckBox("Highlight only after drag (not during)")
        self._after_drag_cb.setToolTip(
            "When highlighting is on, keep disks their normal fill while dragging; "
            "show the red tint only once the drag settles."
        )
        self._after_drag_cb.setChecked(self._options.highlight_only_after_drag)
        self._after_drag_cb.toggled.connect(self._on_after_drag)
        form.addWidget(self._after_drag_cb)

        self._snap_cb = QtWidgets.QCheckBox("Snap back on overlap")
        self._snap_cb.setToolTip("Reserved: revert a move that creates an overlap (not yet wired).")
        self._snap_cb.setChecked(self._options.snap_back_on_overlap)
        self._snap_cb.toggled.connect(self._on_snap_back)
        form.addWidget(self._snap_cb)

    # -- change handlers ----------------------------------------------------

    def _on_highlight(self, on: bool) -> None:
        self._options.highlight_overlaps = bool(on)
        self.optionsChanged.emit()

    def _on_after_drag(self, on: bool) -> None:
        self._options.highlight_only_after_drag = bool(on)
        self.optionsChanged.emit()

    def _on_snap_back(self, on: bool) -> None:
        self._options.snap_back_on_overlap = bool(on)
        self.optionsChanged.emit()


__all__ = ["OptionsPanel", "EditorOptions"]
