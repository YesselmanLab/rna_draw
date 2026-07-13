"""`MainWindow`: toolbar + graphics view + one folding LEFT control panel.

The window is a thin shell: it owns an `EditorModel`, hands it to the
`RnaGraphicsView`, and reflects the never-silent verdict in the status bar.
All controls live in a single LEFT dock of collapsible sections (Structure &
Sequence, Nucleotides, Connectors, Colors / scheme, Canvas, Editing) built on
`StylePanel` + injected sections -- there is no separate right Options dock.
The persistent hint keeps the core interaction discoverable at all times.
"""

from __future__ import annotations

import json
from datetime import datetime

from PySide6 import QtCore, QtGui, QtWidgets

from rna_draw.document import Document
from rna_draw.gui.model import DEMO_SEQ, DEMO_SS, EditorModel
from rna_draw.io_formats import OPEN_FILTER, parse_structure_file

from .collapsible import CollapsibleSection
from .options_panel import OptionsPanel
from .scene_view import RnaGraphicsView
from .style_panel import ColorButton, StylePanel

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
        # Path of the currently open `.rnadoc.json` (Save target); None until a
        # document is opened or saved-as.
        self._current_path: str | None = None
        # In-memory version-history snapshots (each a full `Document`), also
        # persisted into a saved file's `extra["versions"]`.
        self._versions: list[Document] = []

        self._view = RnaGraphicsView()
        self.setCentralWidget(self._view)
        self._view.selectionChanged.connect(self._on_selection)
        self._view.overlapChanged.connect(self._on_overlap)

        self._build_toolbar()
        self._build_file_menu()
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
        open_act.setToolTip(
            "Open a .rnadoc.json document or an RNA structure file "
            "(.dbn, .dot, .ct, .bpseq, .fasta, .txt)"
        )
        open_act.triggered.connect(self._open_file)
        bar.addAction(open_act)

        save_act = QtGui.QAction("Save", self)
        save_act.setToolTip("Save the current drawing as a .rnadoc.json document")
        save_act.triggered.connect(self._save_document)
        bar.addAction(save_act)

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

    def _build_file_menu(self) -> None:
        """The File menu: Open / Save / Save As over the `.rnadoc` document layer."""
        menu = self.menuBar().addMenu("&File")

        open_act = QtGui.QAction("&Open...", self)
        open_act.setShortcut(QtGui.QKeySequence.StandardKey.Open)
        open_act.triggered.connect(self._open_file)
        menu.addAction(open_act)

        menu.addSeparator()

        save_act = QtGui.QAction("&Save", self)
        save_act.setShortcut(QtGui.QKeySequence.StandardKey.Save)
        save_act.triggered.connect(self._save_document)
        menu.addAction(save_act)

        save_as_act = QtGui.QAction("Save &As...", self)
        save_as_act.setShortcut(QtGui.QKeySequence.StandardKey.SaveAs)
        save_as_act.triggered.connect(self._save_document_as)
        menu.addAction(save_as_act)

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
        self._build_selection_styling_section()  # -> bottom of the panel
        self._build_versions_section()  # -> bottom of the panel
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

    def _build_selection_styling_section(self) -> None:
        """Per-region styling: recolor the selection or add a highlight halo.

        Enabled only when a selection exists (any granularity). "Color
        selection" applies a per-nt color override (wins over the scheme);
        "Highlight selection" adds a persistent translucent halo. Both are
        purely visual -- no geometry moves, so the never-silent contract is
        untouched. Wired live: apply -> model update -> canvas rebuild.
        """
        section = CollapsibleSection("Selection styling", expanded=False)
        col = QtWidgets.QVBoxLayout()
        col.setContentsMargins(6, 4, 6, 8)
        col.setSpacing(6)

        self._sel_readout = QtWidgets.QLabel("No selection")
        self._sel_readout.setStyleSheet("color:#7a828d; font-style:italic;")
        col.addWidget(self._sel_readout)

        # Recolor row: swatch + apply button.
        color_row = QtWidgets.QWidget()
        color_lay = QtWidgets.QHBoxLayout(color_row)
        color_lay.setContentsMargins(0, 0, 0, 0)
        color_lay.setSpacing(6)
        self._sel_color_btn = ColorButton("#e5484d")
        self._sel_color_btn.setToolTip("Color to paint the selected nucleotides")
        color_apply = QtWidgets.QPushButton("Color selection")
        color_apply.setToolTip("Recolor the selected nucleotides (overrides the scheme)")
        color_apply.clicked.connect(self._on_color_selection)
        color_lay.addWidget(self._sel_color_btn)
        color_lay.addWidget(color_apply, 1)
        col.addWidget(color_row)

        # Highlight row: swatch + apply button.
        hl_row = QtWidgets.QWidget()
        hl_lay = QtWidgets.QHBoxLayout(hl_row)
        hl_lay.setContentsMargins(0, 0, 0, 0)
        hl_lay.setSpacing(6)
        self._sel_highlight_btn = ColorButton("#ffd300")
        self._sel_highlight_btn.setToolTip("Halo color for the highlighted region")
        hl_apply = QtWidgets.QPushButton("Highlight selection")
        hl_apply.setToolTip("Add a persistent translucent halo over the selected region")
        hl_apply.clicked.connect(self._on_highlight_selection)
        hl_lay.addWidget(self._sel_highlight_btn)
        hl_lay.addWidget(hl_apply, 1)
        col.addWidget(hl_row)

        # Clear buttons.
        clear_row = QtWidgets.QWidget()
        clear_lay = QtWidgets.QHBoxLayout(clear_row)
        clear_lay.setContentsMargins(0, 0, 0, 0)
        clear_lay.setSpacing(6)
        clear_color = QtWidgets.QPushButton("Clear color")
        clear_color.setToolTip("Revert the selected nucleotides to the scheme color")
        clear_color.clicked.connect(self._on_clear_color)
        clear_hl = QtWidgets.QPushButton("Clear highlights")
        clear_hl.setToolTip("Remove every highlight halo")
        clear_hl.clicked.connect(self._on_clear_highlights)
        clear_lay.addWidget(clear_color)
        clear_lay.addWidget(clear_hl)
        col.addWidget(clear_row)

        # Widgets that only make sense with a live selection.
        self._sel_apply_widgets = [color_apply, hl_apply, clear_color]

        section.set_content_layout(col)
        self._panel.add_section_bottom(section)
        self._update_selection_styling()

    def _update_selection_styling(self) -> None:
        """Refresh the Selection-styling readout + enablement from the model."""
        has_sel = bool(self._model is not None and self._model.sel_indices)
        for w in getattr(self, "_sel_apply_widgets", []):
            w.setEnabled(has_sel)
        if has_sel:
            n = len(self._model.sel_indices)
            kind = self._model.sel_kind or "selection"
            self._sel_readout.setText(f"{n} nt selected ({kind})")
        else:
            self._sel_readout.setText("No selection")

    def _on_color_selection(self) -> None:
        """Apply a per-nt color override to the current selection (visual)."""
        if self._model is None:
            return
        self._model.apply_color_to_selection(self._sel_color_btn.color())
        self._view.restyle()

    def _on_highlight_selection(self) -> None:
        """Add a persistent highlight halo over the current selection (visual)."""
        if self._model is None:
            return
        self._model.highlight_selection(self._sel_highlight_btn.color())
        self._view.restyle()

    def _on_clear_color(self) -> None:
        """Drop color overrides on the current selection (revert to scheme)."""
        if self._model is None:
            return
        self._model.clear_color_on_selection()
        self._view.restyle()

    def _on_clear_highlights(self) -> None:
        """Remove every highlight halo (visual)."""
        if self._model is None:
            return
        self._model.clear_highlights()
        self._view.restyle()

    def _build_versions_section(self) -> None:
        """The 'Versions' section: snapshot the current state + restore snapshots.

        Each snapshot is a full `Document` (coords + style), held in memory and
        round-tripped into a saved file's `extra["versions"]`. Restoring one
        re-runs the never-silent gate (via `EditorModel.from_document`).
        """
        section = CollapsibleSection("Versions", expanded=False)
        col = QtWidgets.QVBoxLayout()
        col.setContentsMargins(6, 4, 6, 8)
        col.setSpacing(6)

        snap_btn = QtWidgets.QPushButton("Snapshot current")
        snap_btn.setToolTip("Capture the current drawing as a version you can jump back to")
        snap_btn.clicked.connect(self._take_snapshot)
        col.addWidget(snap_btn)

        self._versions_list = QtWidgets.QListWidget()
        self._versions_list.setToolTip("Double-click a version to restore it")
        self._versions_list.itemDoubleClicked.connect(self._restore_version_item)
        col.addWidget(self._versions_list)

        restore_btn = QtWidgets.QPushButton("Restore selected")
        restore_btn.clicked.connect(self._restore_selected_version)
        col.addWidget(restore_btn)

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
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open", "", OPEN_FILTER)
        if not path:
            return
        if path.endswith(".rnadoc.json") or path.endswith(".rnadoc"):
            self._open_document(path)
            return
        try:
            ss, seq, _name = parse_structure_file(path)
        except (OSError, ValueError) as exc:
            self._overlap_label.setText(f"Could not read structure: {exc}")
            return
        self.load_ss(ss, seq)

    # -- documents (.rnadoc.json) -------------------------------------------

    def _open_document(self, path: str) -> None:
        """Load a `.rnadoc.json`, restore its layout+style+version history."""
        try:
            doc = Document.load(path)
        except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
            self._overlap_label.setText(f"Could not load document: {exc}")
            return
        self._load_document(doc)
        self._current_path = path
        self._restore_versions_from_doc(doc)

    def _load_document(self, doc: Document) -> None:
        """Adopt `doc` into a fresh model + view (never-silent gate re-run).

        The model's `from_document` re-validates the stored coordinates live;
        an overlapping layout is surfaced flagged (status message) rather than
        drawn silently clean.
        """
        try:
            model = EditorModel.from_document(doc)
        except Exception as exc:  # keep the app alive on a malformed document
            self._overlap_label.setText(f"Could not restore document: {exc}")
            return
        model.set_numbering(self._show_numbers, self._number_interval)
        self._model = model
        self._ss_edit.setText(model._ss or "")
        self._seq_edit.setText(model._seq or "")
        self._view.set_model(model)
        # Adopt the document's style into the panel (so the controls reflect it
        # and stay bound to the same preset the model draws with).
        self._panel.load_preset(model.preset)
        self._engine_label.setText(
            f"{len(model.scene()['nucleotides'])} nt  |  engine: {model.engine_name}"
        )
        self._on_overlap(model.flagged, len(model.scene()["overlaps"]))
        if model.load_note:
            self.statusBar().showMessage(model.load_note, 8000)

    def _save_document(self) -> None:
        """Save to the current path, or prompt for one on first save."""
        if self._current_path:
            self._write_document(self._current_path)
        else:
            self._save_document_as()

    def _save_document_as(self) -> None:
        """Prompt for a path and save the drawing (+ version history) there."""
        if self._model is None:
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save RNA document", "drawing.rnadoc.json", "RNA documents (*.rnadoc.json)"
        )
        if not path:
            return
        if not (path.endswith(".rnadoc.json") or path.endswith(".rnadoc")):
            path += ".rnadoc.json"
        self._current_path = path
        self._write_document(path)

    def _write_document(self, path: str) -> None:
        """Serialize the live model (+ snapshot history) to `path`."""
        if self._model is None:
            return
        doc = self._model.to_document()
        doc.extra["versions"] = [snap.to_dict() for snap in self._versions]
        try:
            doc.save(path)
        except OSError as exc:
            self._overlap_label.setText(f"Could not save: {exc}")
            return
        self.statusBar().showMessage(f"Saved {path}", 4000)

    # -- version history ----------------------------------------------------

    def _take_snapshot(self, _checked: bool = False, name: str | None = None) -> None:
        """Capture the current editor state as a named version snapshot."""
        if self._model is None:
            return
        if not name:
            name = f"v{len(self._versions) + 1} - {datetime.now():%H:%M:%S}"
        self._versions.append(self._model.to_document(name=name))
        self._versions_list.addItem(name)

    def _restore_version_item(self, item: QtWidgets.QListWidgetItem) -> None:
        self._restore_version(self._versions_list.row(item))

    def _restore_selected_version(self) -> None:
        self._restore_version(self._versions_list.currentRow())

    def _restore_version(self, row: int) -> None:
        """Load the snapshot at `row` back into the model (never-silent gate)."""
        if 0 <= row < len(self._versions):
            self._load_document(self._versions[row])

    def _restore_versions_from_doc(self, doc: Document) -> None:
        """Rebuild the snapshot list from a loaded file's `extra["versions"]`."""
        self._versions = []
        self._versions_list.clear()
        for entry in doc.extra.get("versions", []):
            try:
                snap = Document.from_dict(entry)
            except (KeyError, TypeError, ValueError):
                continue
            name = snap.extra.get("name") or f"v{len(self._versions) + 1}"
            self._versions.append(snap)
            self._versions_list.addItem(str(name))

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
        self._update_selection_styling()
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
