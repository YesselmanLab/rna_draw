"""Interactive RNA-structure editor frontends for rna_draw.

The Jupyter anywidget (`RnaEditor`) is the first frontend of the planned
cross-platform editor. Everything authoritative -- layout and the
never-silent overlap check -- runs kernel-side on the real
`rna_draw.layout`/`rna_draw.overlap` engines; the browser view is a thin,
dependency-free SVG renderer of the shared `rna_draw.scene.Scene`.
"""

from __future__ import annotations

from .editor_widget import RnaEditor, editor

__all__ = ["RnaEditor", "editor"]
