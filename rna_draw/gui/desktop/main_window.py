"""`MainWindow`: toolbar + graphics view + status bar around `EditorModel`.

The window is a thin shell: it owns an `EditorModel`, hands it to the
`RnaGraphicsView`, and reflects the never-silent verdict in the status bar.
The persistent hint ("Click a helix, then drag its handle to rotate") keeps
the core interaction discoverable at all times.
"""

from __future__ import annotations

from PySide6 import QtCore, QtGui, QtWidgets

from rna_draw.gui.model import DEMO_SS, EditorModel

from .scene_view import RnaGraphicsView

HINT = "Click a helix, then drag its orange handle to rotate it about its junction."


class MainWindow(QtWidgets.QMainWindow):
    """The desktop RNA editor's main window."""

    def __init__(self, ss: str = DEMO_SS, seq: str | None = None) -> None:
        super().__init__()
        self.setWindowTitle("rna_draw -- RNA structure editor")
        self.resize(1000, 760)

        self._view = RnaGraphicsView()
        self.setCentralWidget(self._view)
        self._view.selectionChanged.connect(self._on_selection)
        self._view.overlapChanged.connect(self._on_overlap)

        self._build_toolbar()
        self._build_statusbar()

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
        bar.addWidget(QtWidgets.QLabel(" Structure: "))
        self._ss_edit = QtWidgets.QLineEdit()
        self._ss_edit.setPlaceholderText("Enter dot-bracket, e.g. ((((....))))")
        self._ss_edit.setMinimumWidth(360)
        self._ss_edit.returnPressed.connect(self._draw_from_edit)
        bar.addWidget(self._ss_edit)

        draw_act = QtGui.QAction("Draw", self)
        draw_act.setToolTip("Lay out and draw the entered dot-bracket structure")
        draw_act.triggered.connect(self._draw_from_edit)
        bar.addAction(draw_act)

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
        self._model = model
        self._ss_edit.setText(ss)
        self._view.set_model(model)
        self._engine_label.setText(f"{len(ss)} nt  |  engine: {model.engine_name}")
        self._on_overlap(model.flagged, len(model.scene()["overlaps"]))

    def _draw_from_edit(self) -> None:
        ss = self._ss_edit.text().strip()
        if ss:
            self.load_ss(ss)

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
