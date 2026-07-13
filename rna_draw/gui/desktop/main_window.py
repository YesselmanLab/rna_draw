"""`MainWindow`: toolbar + graphics view + one folding LEFT control panel.

The window is a thin shell: it owns an `EditorModel`, hands it to the
`RnaGraphicsView`, and reflects the never-silent verdict in the status bar.
All controls live in a single LEFT dock of collapsible sections (Structure &
Sequence, Nucleotides, Connectors, Colors / scheme, Canvas, Editing) built on
`StylePanel` + injected sections -- there is no separate right Options dock.
The persistent hint keeps the core interaction discoverable at all times.
"""

from __future__ import annotations

from PySide6 import QtCore, QtGui, QtWidgets

from rna_draw.gui.model import DEMO_SEQ, DEMO_SS, EditorModel

from .collapsible import CollapsibleSection
from .options_panel import OptionsPanel
from .scene_view import RnaGraphicsView
from .style_panel import StylePanel

HINT = "Click a helix, then drag its orange handle to rotate it about its junction."


class MainWindow(QtWidgets.QMainWindow):
    """The desktop RNA editor's main window."""

    def __init__(self, ss: str = DEMO_SS, seq: str | None = None) -> None:
        super().__init__()
        self.setWindowTitle("rna_draw -- RNA structure editor")
        self.resize(1000, 760)

        # The demo opens WITH its sequence so letters show on launch; a
        # user-supplied structure with no sequence stays blank.
        if seq is None and ss == DEMO_SS:
            seq = DEMO_SEQ

        self._model: EditorModel | None = None
        # Residue-numbering state, re-applied to every freshly laid-out model.
        self._show_numbers = True
        self._number_interval = 10

        self._view = RnaGraphicsView()
        self.setCentralWidget(self._view)
        self._view.selectionChanged.connect(self._on_selection)
        self._view.overlapChanged.connect(self._on_overlap)

        self._build_toolbar()
        self._build_left_panel()
        self._build_mode_toolbar()
        self._build_statusbar()

        # hand the view its initial behavior options (highlighting OFF by default)
        self._view.set_options(self._options_panel.options())

        self.load_ss(ss, seq)

    # -- chrome -------------------------------------------------------------

    def _build_toolbar(self) -> None:
        bar = QtWidgets.QToolBar("Main")
        bar.setMovable(False)
        bar.setIconSize(QtCore.QSize(16, 16))
        self.addToolBar(bar)

        open_act = QtGui.QAction("Open", self)
        open_act.setToolTip("Open a text file whose first line is a dot-bracket structure")
        open_act.triggered.connect(self._open_file)
        bar.addAction(open_act)

        fit_act = QtGui.QAction("Fit", self)
        fit_act.setToolTip("Fit the whole structure in view")
        fit_act.triggered.connect(self._view.fit)
        bar.addAction(fit_act)

        bar.addSeparator()
        self._toggle_panel_act = QtGui.QAction("Controls", self)
        self._toggle_panel_act.setToolTip("Show/hide the left controls panel")
        self._toggle_panel_act.setCheckable(True)
        self._toggle_panel_act.setChecked(True)
        self._toggle_panel_act.toggled.connect(self._toggle_panel_dock)
        bar.addAction(self._toggle_panel_act)

    def _build_mode_toolbar(self) -> None:
        """The interaction-mode bar (Move / Select / Edit) + granularity combo.

        The active mode governs how clicks/drags are read in the view (see
        `RnaGraphicsView._select_at`). Move (default) preserves the current
        click-a-helix-then-rotate behavior; Select selects at the chosen
        granularity for inspection (and lets a single residue be dragged). Edit
        is a disabled placeholder so the bar is visibly extensible.
        """
        bar = QtWidgets.QToolBar("Mode")
        bar.setObjectName("modeToolbar")
        bar.setMovable(False)
        self.addToolBar(bar)

        self._mode_group = QtGui.QActionGroup(self)
        self._mode_group.setExclusive(True)

        self._move_act = QtGui.QAction("Move", self)
        self._move_act.setCheckable(True)
        self._move_act.setChecked(True)
        self._move_act.setToolTip("Arrange: click a helix, drag its handle to rotate it")
        self._move_act.triggered.connect(lambda: self._set_mode("move"))

        self._select_act = QtGui.QAction("Select", self)
        self._select_act.setCheckable(True)
        self._select_act.setToolTip(
            "Inspect: click selects at the chosen granularity; drag a residue to move it"
        )
        self._select_act.triggered.connect(lambda: self._set_mode("select"))

        self._edit_act = QtGui.QAction("Edit", self)
        self._edit_act.setCheckable(True)
        self._edit_act.setEnabled(False)  # placeholder: future structure editing
        self._edit_act.setToolTip("Structure editing (coming soon)")

        for act in (self._move_act, self._select_act, self._edit_act):
            self._mode_group.addAction(act)
            bar.addAction(act)

        bar.addSeparator()
        bar.addWidget(QtWidgets.QLabel("Select: "))
        self._gran_combo = QtWidgets.QComboBox()
        self._gran_combo.addItems(["Residue", "Helix", "Motif"])
        self._gran_combo.setCurrentText("Helix")  # default preserves current behavior
        self._gran_combo.setToolTip("What a click selects in Select mode")
        self._gran_combo.currentTextChanged.connect(self._on_granularity)
        bar.addWidget(self._gran_combo)

    def _set_mode(self, mode: str) -> None:
        """Adopt an interaction mode and push it to the view (visual only)."""
        opts = self._options_panel.options()
        opts.mode = mode
        self._view.set_options(opts)

    def _on_granularity(self, text: str) -> None:
        """Adopt a selection granularity and push it to the view (visual only)."""
        opts = self._options_panel.options()
        opts.granularity = text.strip().lower()
        self._view.set_options(opts)

    def _build_left_panel(self) -> None:
        """Build the single LEFT dock: StylePanel with injected top/bottom sections."""
        self._panel = StylePanel()
        self._build_structure_section()  # -> top of the panel
        self._build_editing_section()  # -> bottom of the panel

        self._dock = QtWidgets.QDockWidget("Controls", self)
        self._dock.setObjectName("controlsDock")
        self._dock.setAllowedAreas(
            QtCore.Qt.DockWidgetArea.LeftDockWidgetArea
            | QtCore.Qt.DockWidgetArea.RightDockWidgetArea
        )
        self._dock.setWidget(self._panel)
        # A comfortable min width so combos ("res_type"/"paired") never clip.
        self._dock.setMinimumWidth(300)
        self.addDockWidget(QtCore.Qt.DockWidgetArea.LeftDockWidgetArea, self._dock)

        self._panel.visualChanged.connect(self._on_style_visual)
        self._panel.relayoutRequested.connect(self._on_style_relayout)
        self._options_panel.optionsChanged.connect(self._on_options_changed)
        self._dock.visibilityChanged.connect(
            lambda vis: self._toggle_panel_act.setChecked(vis)
        )

        self._view_menu = self.menuBar().addMenu("&View")
        self._view_menu.addAction(self._dock.toggleViewAction())

    def _build_structure_section(self) -> None:
        """VARNA-style Structure + Sequence on separate lines + a Draw button."""
        section = CollapsibleSection("Structure & Sequence", expanded=True)
        form = QtWidgets.QFormLayout()
        form.setLabelAlignment(QtCore.Qt.AlignmentFlag.AlignRight)
        form.setFieldGrowthPolicy(
            QtWidgets.QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow
        )
        form.setContentsMargins(6, 4, 6, 8)
        form.setSpacing(6)

        self._ss_edit = QtWidgets.QLineEdit()
        self._ss_edit.setPlaceholderText("Dot-bracket, e.g. ((((....))))")
        self._ss_edit.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Fixed
        )
        self._ss_edit.returnPressed.connect(self._draw_from_edit)
        form.addRow("Structure:", self._ss_edit)

        self._seq_edit = QtWidgets.QLineEdit()
        self._seq_edit.setPlaceholderText("Optional sequence, e.g. GGGG....CCCC")
        self._seq_edit.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Fixed
        )
        self._seq_edit.returnPressed.connect(self._draw_from_edit)
        form.addRow("Sequence:", self._seq_edit)

        draw_btn = QtWidgets.QPushButton("Draw")
        draw_btn.setToolTip("Lay out and draw the structure with its sequence letters")
        draw_btn.clicked.connect(self._draw_from_edit)
        form.addRow("", draw_btn)

        section.set_content_layout(form)
        self._panel.add_section_top(section)

    def _build_editing_section(self) -> None:
        """Editing behavior toggles (folded in from the old right dock) + numbering."""
        section = CollapsibleSection("Editing", expanded=False)
        col = QtWidgets.QVBoxLayout()
        col.setContentsMargins(6, 4, 6, 8)
        col.setSpacing(6)

        self._options_panel = OptionsPanel()
        col.addWidget(self._options_panel)

        # Residue-numbering controls (drive the model's scene numbers).
        num_row = QtWidgets.QWidget()
        num_lay = QtWidgets.QHBoxLayout(num_row)
        num_lay.setContentsMargins(0, 0, 0, 0)
        num_lay.setSpacing(6)
        self._numbers_cb = QtWidgets.QCheckBox("Show residue numbers")
        self._numbers_cb.setChecked(self._show_numbers)
        self._numbers_cb.toggled.connect(self._on_numbering_changed)
        num_lay.addWidget(self._numbers_cb)
        num_lay.addStretch(1)
        num_lay.addWidget(QtWidgets.QLabel("every"))
        self._number_interval_sp = QtWidgets.QSpinBox()
        self._number_interval_sp.setRange(1, 100)
        self._number_interval_sp.setValue(self._number_interval)
        self._number_interval_sp.valueChanged.connect(self._on_numbering_changed)
        num_lay.addWidget(self._number_interval_sp)
        num_lay.addWidget(QtWidgets.QLabel("nt"))
        col.addWidget(num_row)

        section.set_content_layout(col)
        self._panel.add_section_bottom(section)

    def _toggle_panel_dock(self, on: bool) -> None:
        self._dock.setVisible(on)

    def _on_options_changed(self) -> None:
        """A behavior toggle flipped: re-read options + refresh (visual only)."""
        self._view.set_options(self._options_panel.options())

    def _on_numbering_changed(self, *_args) -> None:
        """Residue-numbering toggle/interval changed: re-apply to the model (visual)."""
        self._show_numbers = self._numbers_cb.isChecked()
        self._number_interval = int(self._number_interval_sp.value())
        if self._model is not None:
            self._model.set_numbering(self._show_numbers, self._number_interval)
            self._view.restyle()

    def _build_statusbar(self) -> None:
        sb = self.statusBar()
        self._engine_label = QtWidgets.QLabel("")
        self._overlap_label = QtWidgets.QLabel("")
        self._hint_label = QtWidgets.QLabel(HINT)
        self._hint_label.setStyleSheet("color: #1fb6a6;")
        sb.addWidget(self._engine_label)
        sb.addWidget(self._overlap_label)
        sb.addPermanentWidget(self._hint_label)

    # -- actions ------------------------------------------------------------

    def load_ss(self, ss: str, seq: str | None = None) -> None:
        """Lay out `ss` and render it, or show the error in the status bar."""
        try:
            model = EditorModel.from_ss(ss, seq=seq)
        except Exception as exc:  # keep the app alive on a bad structure
            self._overlap_label.setText(f"Could not lay out structure: {exc}")
            return
        # carry the current panel styling + numbering onto the fresh model
        model.set_style(self._panel.preset())
        model.set_numbering(self._show_numbers, self._number_interval)
        self._model = model
        self._ss_edit.setText(ss)
        self._seq_edit.setText(seq or "")
        self._view.set_model(model)
        self._engine_label.setText(f"{len(ss)} nt  |  engine: {model.engine_name}")
        self._on_overlap(model.flagged, len(model.scene()["overlaps"]))

    def _draw_from_edit(self) -> None:
        """Apply BOTH the Structure and Sequence fields via a fresh layout."""
        ss = self._ss_edit.text().strip()
        if not ss:
            return
        seq = self._seq_edit.text().strip() or None
        self.load_ss(ss, seq)

    def _apply_sequence(self) -> None:
        """Relabel the current drawing from the Sequence field (no re-layout, no move)."""
        if self._model is None:
            return
        note = self._model.set_sequence(self._seq_edit.text())
        self._view.restyle()
        if note:
            self.statusBar().showMessage(note, 6000)

    def _open_file(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Open dot-bracket structure", "", "Text files (*.txt *.dat *.ss);;All files (*)"
        )
        if not path:
            return
        try:
            with open(path, encoding="utf-8") as fh:
                lines = [ln.strip() for ln in fh if ln.strip()]
        except OSError as exc:
            self._overlap_label.setText(f"Could not read file: {exc}")
            return
        if not lines:
            self._overlap_label.setText("File contained no structure line")
            return
        ss = next((ln for ln in lines if set(ln) <= set("().[]{}<>")), lines[0])
        seq = None
        if len(lines) >= 2 and len(lines[0]) == len(lines[1]):
            # convention: seq on the first line, structure on the second
            if set(lines[1]) <= set("().[]{}<>"):
                seq, ss = lines[0], lines[1]
        self.load_ss(ss, seq)

    # -- style editor -------------------------------------------------------

    def _on_style_visual(self) -> None:
        """A live-visual style change: restyle the canvas without re-layout."""
        if self._model is None:
            return
        self._model.set_style(self._panel.preset())
        self._view.restyle()

    def _on_style_relayout(self) -> None:
        """Apply staged geometry (node_r/spacing): re-run layout + re-check."""
        if self._model is None:
            return
        result = self._model.relayout(self._panel.preset())
        self._view.set_model(self._model)
        self._engine_label.setText(
            f"{len(self._model.scene()['nucleotides'])} nt  |  engine: {self._model.engine_name}"
        )
        self._on_overlap(result.flagged, len(result.overlaps))

    # -- signals ------------------------------------------------------------

    def _on_selection(self, selection) -> None:
        if selection is None:
            self._overlap_label.setText("No helix selected")
        else:
            self._on_overlap(self._model.flagged, len(self._model.scene()["overlaps"]))

    def _on_overlap(self, flagged: bool, count: int) -> None:
        if flagged:
            self._overlap_label.setText(f"⚠ {count} overlap(s) flagged")
            self._overlap_label.setStyleSheet("color: #e5484d; font-weight: 600;")
        else:
            self._overlap_label.setText("✓ clean")
            self._overlap_label.setStyleSheet("color: #2ea043; font-weight: 600;")


__all__ = ["MainWindow"]
