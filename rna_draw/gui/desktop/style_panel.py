"""`StylePanel`: the left-dock style editor bound to a `StylePreset`.

A single scrollable widget of grouped controls (nucleotides, connectors,
palette/scheme, canvas) that mutate one live `StylePreset` in place. Visual
changes (colors, widths, letters, background, palette, render-type) emit
`visualChanged` for an immediate, layout-free canvas restyle; geometry knobs
(node_r, spacing) only stage the preset and are applied on the explicit
`relayoutRequested` -- never silently rescaling stored coordinates (B3).

The panel imports Qt but no layout/checker code: it edits a preset and emits;
`MainWindow` wires the emissions to `EditorModel.set_style` / `relayout`.
"""

from __future__ import annotations

from PySide6 import QtCore, QtGui, QtWidgets

from rna_draw.gui.model import DEFAULT_VIEW, view_style
from rna_draw.style import StylePreset, default_preset

from .collapsible import CollapsibleSection

_PALETTE_LABELS = {
    "r": "red",
    "g": "green",
    "b": "blue",
    "y": "yellow",
    "c": "cyan",
    "m": "magenta",
    "w": "white",
    "e": "grey",
    "o": "orange",
}
_RENDER_TYPES = ["none", "res_type", "paired"]
_DATA_PALETTES = ["viridis", "plasma", "inferno", "magma", "cividis", "coolwarm", "Reds"]
_LINE_STYLES = ["solid", "dash", "dot"]


class ColorButton(QtWidgets.QPushButton):
    """A swatch button that opens `QColorDialog` and reports the chosen hex.

    Signals:
        colorChanged: emits the new ``#rrggbb`` string when the user picks one.
    """

    colorChanged = QtCore.Signal(str)

    def __init__(self, color: str = "#ffffff", parent=None) -> None:
        super().__init__(parent)
        self._color = "#ffffff"
        self.setFixedSize(QtCore.QSize(46, 22))
        self.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        self.clicked.connect(self._pick)
        self.set_color(color)

    def color(self) -> str:
        """The current ``#rrggbb`` swatch color."""
        return self._color

    def set_color(self, color: str) -> None:
        """Set the swatch color WITHOUT emitting (for programmatic sync)."""
        c = QtGui.QColor(color)
        if not c.isValid():
            return
        self._color = c.name()
        self.setStyleSheet(
            f"background:{self._color}; border:1px solid #9098a3; border-radius:4px;"
        )

    def _pick(self) -> None:
        chosen = QtWidgets.QColorDialog.getColor(QtGui.QColor(self._color), self, "Pick a color")
        if chosen.isValid():
            self.set_color(chosen.name())
            self.colorChanged.emit(self._color)


class StylePanel(QtWidgets.QWidget):
    """Grouped style controls over one live `StylePreset`.

    Signals:
        visualChanged: a live-visual property changed; restyle the canvas
            without re-laying-out (safe -- no geometry moved).
        relayoutRequested: the user asked to apply staged GEOMETRY changes
            (node_r / spacing); the canvas must re-run layout + re-check.
    """

    visualChanged = QtCore.Signal()
    relayoutRequested = QtCore.Signal()

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._preset: StylePreset = self._fresh_default()
        self._loading = False
        self._swatches: dict[str, ColorButton] = {}
        self._sections: dict[str, CollapsibleSection] = {}
        self._build_ui()
        self._sync_widgets()

    # -- public API ---------------------------------------------------------

    def preset(self) -> StylePreset:
        """The live preset (mutated in place by the controls)."""
        return self._preset

    def load_preset(self, preset: StylePreset) -> None:
        """Adopt `preset`, refresh every control, and signal a visual restyle."""
        self._preset = preset
        self._preset.extra.setdefault("view", dict(DEFAULT_VIEW))
        # ensure all view keys present
        merged = dict(DEFAULT_VIEW)
        merged.update(self._preset.extra["view"])
        self._preset.extra["view"] = merged
        self._rebuild_palette_group()
        self._sync_widgets()
        self.visualChanged.emit()

    def load_preset_file(self, path: str) -> None:
        """Load a `.rnastyle.json` from `path` (reuses `StylePreset.load`)."""
        self.load_preset(StylePreset.load(path))

    def save_preset_file(self, path: str) -> None:
        """Save the live preset to `path` as `.rnastyle.json`."""
        self._preset.save(path)

    def reset(self) -> None:
        """Restore the built-in default preset and signal a restyle."""
        self.load_preset(self._fresh_default())

    @staticmethod
    def _fresh_default() -> StylePreset:
        """Default preset with the panel's view extras + letters on by default."""
        preset = default_preset()
        preset.extra["view"] = dict(DEFAULT_VIEW)
        preset.layout_defaults.render_in_letters = True
        return preset

    # -- convenience setters (also used by tests) ---------------------------

    def _view(self) -> dict:
        return self._preset.extra.setdefault("view", dict(DEFAULT_VIEW))

    def set_render_type(self, name: str) -> None:
        """Set the coloring scheme (none/res_type/paired) and restyle."""
        self._view()["render_type"] = name
        self._sync_widgets()
        self.visualChanged.emit()

    def set_palette_color(self, key: str, hex_color: str) -> None:
        """Recolor palette entry `key` and restyle."""
        rgb = QtGui.QColor(hex_color)
        self._preset.palette[key] = [rgb.redF(), rgb.greenF(), rgb.blueF()]
        if key in self._swatches:
            self._swatches[key].set_color(hex_color)
        self.visualChanged.emit()

    def set_show_letters(self, on: bool) -> None:
        """Toggle nucleotide letters and restyle."""
        self._preset.layout_defaults.render_in_letters = bool(on)
        self._sync_widgets()
        self.visualChanged.emit()

    def set_default_fill(self, hex_color: str) -> None:
        """Set the default nucleotide fill (hex) and restyle."""
        self._preset.default_color = hex_color
        self._sync_widgets()
        self.visualChanged.emit()

    def set_node_r(self, value: float) -> None:
        """Stage a new disk radius (geometry); apply via `relayoutRequested`."""
        self._preset.layout_defaults.node_r = float(value)
        self._sync_widgets()

    # -- construction -------------------------------------------------------

    def _build_ui(self) -> None:
        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        body = QtWidgets.QWidget()
        self._form = QtWidgets.QVBoxLayout(body)
        self._form.setContentsMargins(8, 8, 8, 8)
        self._form.setSpacing(4)
        scroll.setWidget(body)
        outer.addWidget(scroll, 1)

        # Structure & Sequence (top) and Editing (bottom) sections are injected
        # by MainWindow via `add_section_top` / `add_section_bottom`; the style
        # groups sit between them. Nucleotides opens by default, rest folded.
        self._build_nucleotides_group()
        self._build_connectors_group()
        self._build_scheme_group()
        self._build_canvas_group()
        self._form.addStretch(1)
        self._build_buttons()
        outer.addWidget(self._button_row)

    def section(self, title: str) -> CollapsibleSection | None:
        """The `CollapsibleSection` for `title` (for tests / default folding)."""
        return self._sections.get(title)

    def add_section_top(self, section: CollapsibleSection) -> None:
        """Insert an externally-built section above the style groups."""
        self._sections[section._title] = section
        self._form.insertWidget(0, section)

    def add_section_bottom(self, section: CollapsibleSection) -> None:
        """Insert an externally-built section below the style groups (above the stretch)."""
        self._sections[section._title] = section
        # The trailing item in `_form` is the stretch; insert just before it.
        self._form.insertWidget(self._form.count() - 1, section)

    def _group(self, title: str, expanded: bool = False) -> QtWidgets.QFormLayout:
        section = CollapsibleSection(title, expanded=expanded)
        form = QtWidgets.QFormLayout()
        form.setLabelAlignment(QtCore.Qt.AlignmentFlag.AlignRight)
        form.setFieldGrowthPolicy(
            QtWidgets.QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow
        )
        form.setContentsMargins(6, 4, 6, 8)
        form.setSpacing(6)
        section.set_content_layout(form)
        self._sections[title] = section
        self._form.addWidget(section)
        return form

    def _color_row(self, on_change) -> ColorButton:
        btn = ColorButton()
        btn.colorChanged.connect(on_change)
        return btn

    def _vset(self, key: str):
        """A change-handler that writes `key` into the view extras."""
        return lambda value: self._set_view(key, value)

    def _dspin(self, lo, hi, step, on_change) -> QtWidgets.QDoubleSpinBox:
        sp = QtWidgets.QDoubleSpinBox()
        sp.setRange(lo, hi)
        sp.setSingleStep(step)
        sp.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Fixed
        )
        sp.valueChanged.connect(on_change)
        return sp

    def _combo(self, items) -> QtWidgets.QComboBox:
        """A combobox that expands to fill its field so long item text never clips."""
        cb = QtWidgets.QComboBox()
        cb.addItems(items)
        cb.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Fixed
        )
        cb.setMinimumContentsLength(10)
        return cb

    def _build_nucleotides_group(self) -> None:
        f = self._group("Nucleotides", expanded=True)
        self._default_fill_btn = self._color_row(self._on_default_fill)
        f.addRow("Default fill", self._default_fill_btn)
        self._edge_color_btn = self._color_row(lambda c: self._set_view("nt_edge_color", c))
        f.addRow("Outline color", self._edge_color_btn)
        self._edge_width_sp = self._dspin(0.0, 8.0, 0.5, self._vset("nt_edge_width"))
        f.addRow("Outline width", self._edge_width_sp)
        self._letters_cb = QtWidgets.QCheckBox("Show letters")
        self._letters_cb.toggled.connect(self._on_letters)
        f.addRow("", self._letters_cb)
        self._text_size_sp = self._dspin(4.0, 150.0, 2.0, self._on_text_size)
        f.addRow("Letter size", self._text_size_sp)
        self._letter_color_btn = self._color_row(lambda c: self._set_view("letter_color", c))
        f.addRow("Letter color", self._letter_color_btn)

    def _build_connectors_group(self) -> None:
        f = self._group("Connectors")
        self._pair_color_btn = self._color_row(lambda c: self._set_connector("nested_pair", c))
        f.addRow("Base-pair color", self._pair_color_btn)
        self._pair_width_sp = self._dspin(0.2, 8.0, 0.2, self._vset("pair_width"))
        f.addRow("Base-pair width", self._pair_width_sp)
        self._pk_conn_btn = self._color_row(lambda c: self._set_connector("pk_connector", c))
        f.addRow("PK connector color", self._pk_conn_btn)
        self._pk_line_btn = self._color_row(lambda c: self._set_connector("pk_line", c))
        f.addRow("PK routed color", self._pk_line_btn)
        self._routed_width_sp = self._dspin(0.2, 8.0, 0.2, self._vset("routed_width"))
        f.addRow("PK routed width", self._routed_width_sp)
        self._routed_style_cb = self._combo(_LINE_STYLES)
        self._routed_style_cb.currentTextChanged.connect(self._vset("routed_style"))
        f.addRow("PK routed style", self._routed_style_cb)
        self._backbone_color_btn = self._color_row(lambda c: self._set_view("backbone_color", c))
        f.addRow("Backbone color", self._backbone_color_btn)
        self._backbone_width_sp = self._dspin(0.2, 8.0, 0.2, self._vset("backbone_width"))
        f.addRow("Backbone width", self._backbone_width_sp)

    def _build_scheme_group(self) -> None:
        f = self._group("Colors / scheme")
        self._render_cb = self._combo(_RENDER_TYPES)
        self._render_cb.currentTextChanged.connect(self._on_render_type)
        f.addRow("Render type", self._render_cb)
        self._data_pal_cb = self._combo(_DATA_PALETTES)
        self._data_pal_cb.currentTextChanged.connect(lambda t: self._set_view("data_palette", t))
        f.addRow("Data palette", self._data_pal_cb)
        # editable single-letter palette swatches, rebuilt on preset load.
        self._palette_box = QtWidgets.QGroupBox("Palette")
        self._palette_grid = QtWidgets.QGridLayout(self._palette_box)
        self._palette_grid.setSpacing(4)
        f.addRow(self._palette_box)
        self._rebuild_palette_group()

    def _rebuild_palette_group(self) -> None:
        while self._palette_grid.count():
            item = self._palette_grid.takeAt(0)
            w = item.widget()
            if w is not None:
                w.deleteLater()
        self._swatches = {}
        keys = [k for k in self._preset.palette if len(k) == 1]
        for i, key in enumerate(keys):
            label = QtWidgets.QLabel(f"{key} ({_PALETTE_LABELS.get(key, key)})")
            btn = ColorButton()
            btn.colorChanged.connect(lambda c, k=key: self.set_palette_color(k, c))
            self._swatches[key] = btn
            row, col = divmod(i, 2)
            self._palette_grid.addWidget(label, row, col * 2)
            self._palette_grid.addWidget(btn, row, col * 2 + 1)

    def _build_canvas_group(self) -> None:
        f = self._group("Canvas")
        self._bg_btn = self._color_row(lambda c: self._set_view("background", c))
        f.addRow("Background", self._bg_btn)
        note = QtWidgets.QLabel("Geometry below needs Relayout:")
        note.setStyleSheet("color:#7a828d; font-style:italic;")
        f.addRow(note)
        self._node_r_sp = self._dspin(2.0, 40.0, 1.0, self._on_node_r)
        f.addRow("Disk size (node_r)", self._node_r_sp)
        self._primary_sp = self._dspin(4.0, 80.0, 1.0, self._on_primary_space)
        f.addRow("Primary space", self._primary_sp)
        self._pair_sp = self._dspin(4.0, 80.0, 1.0, self._on_pair_space)
        f.addRow("Pair space", self._pair_sp)
        self._relayout_btn = QtWidgets.QPushButton("Relayout")
        self._relayout_btn.setToolTip("Re-run layout at the new geometry (re-checks overlaps)")
        self._relayout_btn.clicked.connect(self.relayoutRequested.emit)
        f.addRow("", self._relayout_btn)

    def _build_buttons(self) -> None:
        self._button_row = QtWidgets.QWidget()
        row = QtWidgets.QHBoxLayout(self._button_row)
        row.setContentsMargins(10, 6, 10, 10)
        reset_btn = QtWidgets.QPushButton("Reset")
        reset_btn.clicked.connect(self.reset)
        load_btn = QtWidgets.QPushButton("Load...")
        load_btn.clicked.connect(self._load_dialog)
        save_btn = QtWidgets.QPushButton("Save...")
        save_btn.clicked.connect(self._save_dialog)
        row.addWidget(reset_btn)
        row.addWidget(load_btn)
        row.addWidget(save_btn)

    # -- widget <-> preset sync ---------------------------------------------

    def _sync_widgets(self) -> None:
        """Push every preset value into the widgets without re-emitting."""
        self._loading = True
        try:
            resolved = view_style(self._preset)
            ld = self._preset.layout_defaults
            view = {**DEFAULT_VIEW, **self._preset.extra.get("view", {})}
            self._default_fill_btn.set_color(resolved["default_fill"])
            self._edge_color_btn.set_color(view["nt_edge_color"])
            self._edge_width_sp.setValue(float(view["nt_edge_width"]))
            self._letters_cb.setChecked(bool(ld.render_in_letters))
            self._text_size_sp.setValue(float(ld.text_size))
            self._letter_color_btn.set_color(view["letter_color"])
            self._pair_color_btn.set_color(resolved["pair_color"])
            self._pair_width_sp.setValue(float(view["pair_width"]))
            self._pk_conn_btn.set_color(resolved["crossing_color"])
            self._pk_line_btn.set_color(resolved["routed_color"])
            self._routed_width_sp.setValue(float(view["routed_width"]))
            self._routed_style_cb.setCurrentText(view["routed_style"])
            self._backbone_color_btn.set_color(view["backbone_color"])
            self._backbone_width_sp.setValue(float(view["backbone_width"]))
            self._render_cb.setCurrentText(view["render_type"])
            self._data_pal_cb.setCurrentText(view["data_palette"])
            self._bg_btn.set_color(view["background"])
            self._node_r_sp.setValue(float(ld.node_r))
            self._primary_sp.setValue(float(ld.primary_space))
            self._pair_sp.setValue(float(ld.pair_space))
            for key, btn in self._swatches.items():
                rgb = self._preset.palette.get(key)
                if rgb is not None:
                    btn.set_color(QtGui.QColor.fromRgbF(*rgb[:3]).name())
        finally:
            self._loading = False

    # -- change handlers (emit only on genuine user edits) ------------------

    def _set_view(self, key: str, value) -> None:
        if self._loading:
            return
        self._view()[key] = value
        self.visualChanged.emit()

    def _set_connector(self, key: str, hex_color: str) -> None:
        if self._loading:
            return
        self._preset.connector_colors[key] = hex_color
        self.visualChanged.emit()

    def _on_default_fill(self, hex_color: str) -> None:
        if self._loading:
            return
        self._preset.default_color = hex_color
        self.visualChanged.emit()

    def _on_letters(self, on: bool) -> None:
        if self._loading:
            return
        self._preset.layout_defaults.render_in_letters = bool(on)
        self.visualChanged.emit()

    def _on_text_size(self, value: float) -> None:
        if self._loading:
            return
        self._preset.layout_defaults.text_size = float(value)
        self.visualChanged.emit()

    def _on_render_type(self, name: str) -> None:
        if self._loading:
            return
        self._view()["render_type"] = name
        self.visualChanged.emit()

    def _on_node_r(self, value: float) -> None:
        if self._loading:
            return
        self._preset.layout_defaults.node_r = float(value)

    def _on_primary_space(self, value: float) -> None:
        if self._loading:
            return
        self._preset.layout_defaults.primary_space = float(value)

    def _on_pair_space(self, value: float) -> None:
        if self._loading:
            return
        self._preset.layout_defaults.pair_space = float(value)

    # -- file dialogs -------------------------------------------------------

    def _load_dialog(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Load style preset", "", "RNA style (*.rnastyle.json *.json);;All files (*)"
        )
        if path:
            try:
                self.load_preset_file(path)
            except (OSError, ValueError):
                pass

    def _save_dialog(self) -> None:
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save style preset", "style.rnastyle.json", "RNA style (*.rnastyle.json *.json)"
        )
        if path:
            try:
                self.save_preset_file(path)
            except OSError:
                pass


__all__ = ["StylePanel", "ColorButton"]
