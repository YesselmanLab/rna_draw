"""`CollapsibleSection`: a foldable titled container for the left control panel.

A header `QToolButton` (a triangle arrow -- pointing down when expanded,
right when collapsed -- next to a bold title) toggles the visibility of a
single content widget. Several sections stack in one scrolling column and any
number may be open at once, replacing the old two-dock split with one panel.

The section is intentionally dumb: it owns no app state, just shows/hides its
content. Callers put a layout of controls inside via `set_content_layout`
(or drop widgets in with `add_widget`) and read/set the fold state with
`is_expanded` / `set_expanded`.
"""

from __future__ import annotations

from PySide6 import QtCore, QtWidgets


class CollapsibleSection(QtWidgets.QWidget):
    """A titled, foldable section wrapping one content widget.

    Signals:
        toggled: emits the new expanded state (bool) when the user folds or
            unfolds the section.
    """

    toggled = QtCore.Signal(bool)

    def __init__(self, title: str, expanded: bool = True, parent=None) -> None:
        super().__init__(parent)
        self._title = title

        self._button = QtWidgets.QToolButton()
        self._button.setText(title)
        self._button.setCheckable(True)
        self._button.setChecked(expanded)
        self._button.setToolButtonStyle(QtCore.Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self._button.setArrowType(
            QtCore.Qt.ArrowType.DownArrow if expanded else QtCore.Qt.ArrowType.RightArrow
        )
        self._button.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Fixed
        )
        self._button.setStyleSheet(
            "QToolButton { border: none; font-weight: 600; padding: 4px 2px; text-align: left; }"
        )
        self._button.toggled.connect(self._on_toggled)

        self._content = QtWidgets.QWidget()
        self._content.setVisible(expanded)

        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)
        outer.addWidget(self._button)
        outer.addWidget(self._content)

    # -- content ------------------------------------------------------------

    def set_content_layout(self, layout: QtWidgets.QLayout) -> None:
        """Install `layout` as the section's content (replacing any prior one)."""
        old = self._content.layout()
        if old is not None:
            # Reparent the old layout onto a throwaway widget so Qt drops it.
            QtWidgets.QWidget().setLayout(old)
        self._content.setLayout(layout)

    def add_widget(self, widget: QtWidgets.QWidget) -> None:
        """Append `widget` to the content, creating a vbox layout on demand."""
        if self._content.layout() is None:
            self._content.setLayout(QtWidgets.QVBoxLayout())
        self._content.layout().addWidget(widget)

    def content_widget(self) -> QtWidgets.QWidget:
        """The widget whose visibility the header toggles."""
        return self._content

    # -- fold state ---------------------------------------------------------

    def is_expanded(self) -> bool:
        """Whether the section is currently open."""
        return self._button.isChecked()

    def set_expanded(self, on: bool) -> None:
        """Open (`True`) or fold (`False`) the section programmatically."""
        self._button.setChecked(bool(on))

    def _on_toggled(self, on: bool) -> None:
        self._content.setVisible(on)
        self._button.setArrowType(
            QtCore.Qt.ArrowType.DownArrow if on else QtCore.Qt.ArrowType.RightArrow
        )
        self.toggled.emit(on)


__all__ = ["CollapsibleSection"]
