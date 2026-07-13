"""`EditorModel`: the frontend-agnostic RNA-helix editing controller.

This holds the authoritative geometry (`_x`, `_y`, the structure tree) and
performs the only two rigid edits the editors expose -- rotate and translate
a selected helix about its junction pivot. It imports NO GUI toolkit, so the
editing logic is unit-testable headlessly (see `tests/test_desktop_mvp.py`);
the Qt view in `rna_draw.gui.desktop` and the Jupyter widget in
`rna_draw.gui.editor_widget` are thin shells over this same logic.

Never-silent contract (identical to the Jupyter widget): every committed move
is re-checked in-process by the REAL native overlap kernel
(`rna_draw.overlap_native.check_overlaps_native`, the exact twin of the frozen
arbiter) at the layout's own SCALED `node_r` -- the same radius the pipeline
gated at, NOT a fresh default. Overlaps are flagged (never presented as
clean); offending nucleotide indices come from the definitional checker's
witnesses so a view can tint them red. Coordinates are ENGINE space (y-up); a
view flips y for screen.
"""

from __future__ import annotations

from dataclasses import dataclass

from rna_draw.gui.hinge import redistribute_loop
from rna_draw.layout.base import params_at_node_r
from rna_draw.layout.pipeline import layout_guaranteed, resolve_engine
from rna_draw.layout.postpass import (
    rotate_range,
    translate_range,
    witness_nucleotides,
)
from rna_draw.layout.structure_tree import Loop, build_structure_tree, loop_center
from rna_draw.overlap import OverlapParams, check_overlaps
from rna_draw.overlap_native import check_overlaps_native
from rna_draw.render_rna import get_pairmap_from_secstruct
from rna_draw.scene import build_scene

# The target render radius; the pipeline scales down from here when gating.
_TARGET_NODE_R = 10.0

# The demo structure shown on launch so the window is populated immediately: a
# three-way junction (three hairpins off a central multiloop).
DEMO_SS = "((((...((((....))))...((((....))))...))))"


@dataclass
class Result:
    """Outcome of a committed rigid move.

    Args:
        flagged: True iff the native checker found >=1 overlap in the
            committed geometry (never presented as clean when True).
        overlaps: Offending nucleotide indices (for red tinting), empty
            when `flagged` is False.
    """

    flagged: bool
    overlaps: list[int]


class EditorModel:
    """Authoritative click-to-select, rotate/translate RNA-helix editor.

    Pure logic, no Qt import. Holds engine-space coordinates plus the
    topological structure tree and mutates them via exactly-rigid moves,
    re-checking every commit with the native kernel.
    """

    def __init__(self) -> None:
        self._x: list[float] = []
        self._y: list[float] = []
        self._node_r: float = _TARGET_NODE_R
        self._pair_map: list[int] = []
        self._seq: str | None = None
        self._colors = None
        self._crossing_pairs: list = []
        self._crossing_lines: list = []
        self._tree = None
        self._engine_name: str = ""
        # Current selection state.
        self._sel_start: int | None = None
        self._sel_end: int | None = None
        self._pivot: tuple[float, float] | None = None
        # Parent loop of the current selection, redistributed on rotation so
        # the loop stays a clean circle (VARNA hinge; see rna_draw.gui.hinge).
        self._sel_loop: Loop | None = None
        # Cached checker verdict for the current geometry.
        self._flagged: bool = False
        self._overlaps: list[int] = []

    # -- construction -------------------------------------------------------

    @classmethod
    def from_ss(cls, ss: str, seq: str | None = None, engine: str = "auto") -> "EditorModel":
        """Build a model by running the guaranteed layout on a dot-bracket.

        Args:
            ss: Dot-bracket secondary structure.
            seq: Optional sequence for per-nt labels; must match ``len(ss)``.
            engine: Layout-engine selection string for `resolve_engine`.

        Returns:
            A ready-to-edit `EditorModel`.
        """
        result = layout_guaranteed(
            ss, engine=resolve_engine(engine), params=OverlapParams(node_r=_TARGET_NODE_R)
        )
        pair_map = result.pair_map
        if pair_map is None:
            pair_map = get_pairmap_from_secstruct(ss)

        colors = None
        if seq is not None:
            try:
                from rna_draw.colorer import Colorer

                colors = Colorer().get_rgb_colors(seq, ss, default_color=None)
            except Exception:
                colors = None

        model = cls()
        model._x = [float(v) for v in result.x]
        model._y = [float(v) for v in result.y]
        model._node_r = float(result.node_r)
        model._pair_map = [int(v) for v in pair_map]
        model._seq = seq
        model._colors = colors
        model._crossing_pairs = list(result.crossing_pairs or [])
        model._crossing_lines = list(result.crossing_lines or [])
        model._tree = build_structure_tree(model._pair_map)
        model._engine_name = result.engine_name
        # Honest initial verdict from the real checker.
        model._flagged = bool(result.flagged) or len(model._overlap_indices()) > 0
        model._overlaps = model._overlap_indices() if model._flagged else []
        return model

    # -- read-only accessors ------------------------------------------------

    @property
    def engine_name(self) -> str:
        """Name of the engine that produced the current geometry."""
        return self._engine_name

    @property
    def flagged(self) -> bool:
        """Whether the current committed geometry has overlaps."""
        return self._flagged

    @property
    def selection(self) -> tuple[int, int] | None:
        """The selected helix slice ``(start, end)`` or ``None``."""
        if self._sel_start is None or self._sel_end is None:
            return None
        return self._sel_start, self._sel_end

    @property
    def pivot(self) -> tuple[float, float] | None:
        """The current rotation pivot in engine space, or ``None``."""
        return self._pivot

    # -- checker plumbing ---------------------------------------------------

    def _scaled_params(self) -> OverlapParams:
        """Checker params at the layout's actual ``node_r`` (the gate radius)."""
        return params_at_node_r(OverlapParams(node_r=_TARGET_NODE_R), self._node_r)

    def _overlap_indices(self) -> list[int]:
        """Offending nucleotide indices from the definitional checker's witnesses."""
        if len(self._x) < 2:
            return []
        report = check_overlaps(self._x, self._y, self._pair_map, self._scaled_params())
        idx: set[int] = set()
        for w in report.witnesses:
            ka, kb = witness_nucleotides(w, self._pair_map)
            idx.add(ka)
            idx.add(kb)
        return sorted(idx)

    def _centroid(self) -> tuple[float, float]:
        """Fallback pivot: the whole-structure centroid."""
        n = len(self._x)
        if n == 0:
            return 0.0, 0.0
        return sum(self._x) / n, sum(self._y) / n

    # -- selection ----------------------------------------------------------

    def select(self, index: int) -> tuple[int, int, tuple[float, float]] | None:
        """Select the innermost helix containing ``index`` and set its pivot.

        Args:
            index: Clicked nucleotide index.

        Returns:
            ``(start, end, (cx, cy))`` -- the contiguous helix slice and its
            parent-loop-center pivot -- or ``None`` if ``index`` is an
            exterior/unpaired nucleotide inside no helix (selection cleared).
        """
        if self._tree is None or not (0 <= index < len(self._pair_map)):
            self._clear_selection()
            return None
        stack = self._tree.enclosing_pairs.get(index, [])
        if not stack:
            self._clear_selection()
            return None
        start, end = stack[-1]
        parent_pair = stack[-2] if len(stack) >= 2 else None
        parent_loop = self._tree.loop_by_closing_pair.get(parent_pair)
        if parent_loop is not None and parent_loop.members:
            try:
                pivot = loop_center(parent_loop.members, self._x, self._y)
            except ValueError:
                pivot = self._centroid()
        else:
            pivot = self._centroid()
        self._sel_start, self._sel_end, self._pivot = start, end, pivot
        self._sel_loop = parent_loop
        return start, end, pivot

    def deselect(self) -> None:
        """Clear the current selection and pivot."""
        self._clear_selection()

    def _clear_selection(self) -> None:
        self._sel_start = self._sel_end = None
        self._pivot = None
        self._sel_loop = None

    # -- rigid moves (never-silent) -----------------------------------------

    def rotate_selection(self, angle: float, redistribute: bool = True) -> Result:
        """Rotate the selected helix about its pivot, redistribute, re-check, commit.

        The helix's rigid slice is rotated about the parent-loop center
        exactly as before; then, unless ``redistribute`` is False, the
        parent loop's unpaired members are re-spaced evenly on the loop
        circle so the loop stays a clean, compact ring instead of distorting
        (the VARNA hinge; see `rna_draw.gui.hinge.redistribute_loop`). Only
        the dragged helix moves -- sibling helices and the loop's closing
        pair are untouched. When the parent loop has no unpaired members
        (e.g. a stacked-pair continuation), redistribution is a no-op and
        this is a pure rigid rotation.

        Args:
            angle: Rotation angle in radians (counter-clockwise positive).
            redistribute: When True (default), re-space the parent loop's
                unpaired members after rotating. Pass False for a pure
                rigid-slice rotation (the pre-hinge behaviour).

        Returns:
            A `Result` reflecting the native checker's verdict on the
            committed geometry.
        """
        if self._sel_start is None or self._pivot is None:
            return Result(self._flagged, list(self._overlaps))
        cx, cy = self._pivot
        nx, ny = rotate_range(self._x, self._y, self._sel_start, self._sel_end, cx, cy, angle)
        if redistribute and self._sel_loop is not None:
            nx, ny = redistribute_loop(self._sel_loop, nx, ny, (cx, cy))
        return self._commit(nx, ny)

    def translate_selection(self, dx: float, dy: float) -> Result:
        """Rigidly translate the selected helix, re-check, commit.

        Args:
            dx: Translation along engine x.
            dy: Translation along engine y.

        Returns:
            A `Result` reflecting the native checker's verdict.
        """
        if self._sel_start is None:
            return Result(self._flagged, list(self._overlaps))
        nx, ny = translate_range(self._x, self._y, self._sel_start, self._sel_end, dx, dy)
        return self._commit(nx, ny)

    def _commit(self, nx: list[float], ny: list[float]) -> Result:
        """Commit new coordinates after the never-silent native re-check.

        The authoritative verdict always comes from
        `check_overlaps_native` at the scaled gate radius -- a committed
        layout is never reported clean without that check.
        """
        count = check_overlaps_native(nx, ny, self._pair_map, self._scaled_params())
        self._x, self._y = nx, ny
        self._flagged = count > 0
        self._overlaps = self._overlap_indices() if self._flagged else []
        return Result(self._flagged, list(self._overlaps))

    # -- scene --------------------------------------------------------------

    def scene(self) -> dict:
        """Serialize the current geometry + verdict to a portable scene dict."""
        scene = build_scene(
            self._x,
            self._y,
            self._pair_map,
            self._node_r,
            seq=self._seq,
            colors=self._colors,
            crossing_pairs=self._crossing_pairs,
            routed_lines=self._crossing_lines,
            flagged=self._flagged,
            overlaps=self._overlaps,
            pivot=self._pivot,
        )
        return scene.to_dict()


__all__ = ["EditorModel", "Result", "DEMO_SS"]
