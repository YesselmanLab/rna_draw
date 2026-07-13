"""`RnaEditor`: the Jupyter anywidget for interactive RNA-helix editing.

The kernel is the single source of truth. The browser view (``static/
editor.js``) renders a `rna_draw.scene.Scene`, lets the user click a helix and
drag to rotate it, and sends only *intents* (``select`` / ``rotate``) back.
Every committed move is re-checked by the REAL overlap kernel
(`rna_draw.overlap_native.check_overlaps_native`, the exact twin of the
frozen arbiter) at the layout's own scaled ``node_r`` -- overlaps are flagged
(red + status), never silently drawn. The optimistic client-side drag is
cosmetic only; the authoritative scene pushed back after a move always
reflects the checker.
"""

from __future__ import annotations

import math
import pathlib

import anywidget
import traitlets

from rna_draw.layout.base import params_at_node_r
from rna_draw.layout.postpass import rotate_range, witness_nucleotides
from rna_draw.layout.pipeline import layout_guaranteed
from rna_draw.layout.structure_tree import build_structure_tree, loop_center
from rna_draw.overlap import OverlapParams, check_overlaps
from rna_draw.overlap_native import check_overlaps_native
from rna_draw.render_rna import get_pairmap_from_secstruct
from rna_draw.scene import build_scene

_STATIC = pathlib.Path(__file__).parent / "static"
# The target render radius; the pipeline scales down from here when gating.
_TARGET_NODE_R = 10.0


class RnaEditor(anywidget.AnyWidget):
    """Interactive click-to-select, drag-to-rotate RNA helix editor.

    Synced (kernel<->browser): `scene`, `selection`, `status`. Everything
    authoritative (coordinates, structure tree, overlap check) lives in the
    non-synced kernel attributes below and is never trusted from the client.
    """

    _esm = _STATIC / "editor.js"

    scene = traitlets.Dict().tag(sync=True)
    selection = traitlets.List().tag(sync=True)
    status = traitlets.Unicode("").tag(sync=True)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # Authoritative kernel-side model (NOT synced).
        self._x: list[float] = []
        self._y: list[float] = []
        self._node_r: float = _TARGET_NODE_R
        self._pair_map: list[int] = []
        self._seq: str | None = None
        self._colors = None
        self._crossing_pairs = []
        self._crossing_lines = []
        self._tree = None
        # Current selection state.
        self._sel_start: int | None = None
        self._sel_end: int | None = None
        self._pivot: tuple[float, float] | None = None
        self.on_msg(self._handle_msg)

    # -- kernel model setup -------------------------------------------------

    def _load(self, x, y, node_r, pair_map, seq, colors, crossing_pairs, crossing_lines):
        """Populate the authoritative kernel model and build the structure tree."""
        self._x = [float(v) for v in x]
        self._y = [float(v) for v in y]
        self._node_r = float(node_r)
        self._pair_map = [int(v) for v in pair_map]
        self._seq = seq
        self._colors = colors
        self._crossing_pairs = list(crossing_pairs or [])
        self._crossing_lines = list(crossing_lines or [])
        self._tree = build_structure_tree(self._pair_map)

    def _scaled_params(self) -> OverlapParams:
        """Checker params at the layout's actual ``node_r`` (the gate radius)."""
        return params_at_node_r(OverlapParams(node_r=_TARGET_NODE_R), self._node_r)

    def _centroid(self) -> tuple[float, float]:
        """Fallback pivot: the whole-structure centroid."""
        n = len(self._x)
        if n == 0:
            return 0.0, 0.0
        return sum(self._x) / n, sum(self._y) / n

    def _resolve_branch(self, k: int) -> tuple[int, int, tuple[float, float]] | None:
        """Innermost helix containing nucleotide ``k`` and its junction pivot.

        Args:
            k: Clicked nucleotide index.

        Returns:
            ``(start, end, (cx, cy))`` -- the contiguous helix slice ``[start,
            end]`` and its parent-loop-center pivot -- or ``None`` if ``k``
            is an exterior/unpaired nucleotide inside no helix.
        """
        if self._tree is None or not (0 <= k < len(self._pair_map)):
            return None
        stack = self._tree.enclosing_pairs.get(k, [])
        if not stack:
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
        return start, end, pivot

    # -- scene (re)build ----------------------------------------------------

    def _build_scene(self, flagged: bool, overlaps) -> dict:
        """Rebuild the scene dict from the current authoritative coordinates."""
        scene = build_scene(
            self._x,
            self._y,
            self._pair_map,
            self._node_r,
            seq=self._seq,
            colors=self._colors,
            crossing_pairs=self._crossing_pairs,
            routed_lines=self._crossing_lines,
            flagged=flagged,
            overlaps=overlaps,
            pivot=self._pivot,
        )
        return scene.to_dict()

    def _overlap_indices(self) -> list[int]:
        """Offending nucleotide indices from the definitional checker's witnesses."""
        if len(self._x) < 2:
            return []
        report = check_overlaps(self._x, self._y, self._pair_map, self._scaled_params())
        idx = set()
        for w in report.witnesses:
            ka, kb = witness_nucleotides(w, self._pair_map)
            idx.add(ka)
            idx.add(kb)
        return sorted(idx)

    # -- message handling ---------------------------------------------------

    def _handle_msg(self, widget, content, buffers):
        """Dispatch a custom message from the browser view.

        Args:
            widget: The widget instance (ipywidgets passes ``self``).
            content: The message payload dict (``type`` plus fields).
            buffers: Binary buffers (unused).
        """
        msg_type = content.get("type")
        if msg_type == "select":
            self._on_select(int(content["index"]))
        elif msg_type == "rotate":
            self._on_rotate(float(content["angle"]))

    def _on_select(self, k: int):
        """Select the innermost helix containing ``k`` and set its pivot."""
        resolved = self._resolve_branch(k)
        if resolved is None:
            self._sel_start = self._sel_end = None
            self._pivot = None
            self.selection = []
            self.status = f"Nucleotide {k} is not inside a helix"
            return
        start, end, pivot = resolved
        self._sel_start, self._sel_end, self._pivot = start, end, pivot
        self.selection = list(range(start, end + 1))
        # Refresh the scene so the view gets the current pivot for angle math.
        self.scene = self._build_scene(self.scene.get("flagged", False), self.scene.get("overlaps", []))
        self.status = f"Helix {start}-{end} selected"

    def _on_rotate(self, angle: float):
        """Authoritatively rotate the selected helix, re-check, and commit."""
        if self._sel_start is None or self._pivot is None:
            self.status = "Select a helix before rotating"
            return
        start, end = self._sel_start, self._sel_end
        cx, cy = self._pivot
        nx, ny = rotate_range(self._x, self._y, start, end, cx, cy, angle)
        # THE never-silent check: the REAL native kernel twin at the gate radius.
        count = check_overlaps_native(nx, ny, self._pair_map, self._scaled_params())
        # Commit the checked geometry (VARNA-style: flag, do not snap back).
        self._x, self._y = nx, ny
        flagged = count > 0
        overlaps = self._overlap_indices() if flagged else []
        self.scene = self._build_scene(flagged, overlaps)
        deg = math.degrees(angle)
        if flagged:
            self.status = f"Rotated {deg:+.1f} deg -- {count} overlap(s) flagged"
        else:
            self.status = f"Rotated {deg:+.1f} deg -- clean"


def editor(ss: str, seq: str | None = None, engine: str = "auto") -> RnaEditor:
    """Build an interactive `RnaEditor` for a dot-bracket structure.

    Args:
        ss: Dot-bracket secondary structure.
        seq: Optional sequence (per-nucleotide labels); must match ``len(ss)``.
        engine: Layout engine name for `resolve_engine` (default ``"auto"``).

    Returns:
        A ready-to-display `RnaEditor` widget.
    """
    from rna_draw.layout.pipeline import resolve_engine

    result = layout_guaranteed(
        ss, engine=resolve_engine(engine), params=OverlapParams(node_r=_TARGET_NODE_R)
    )
    pair_map = result.pair_map if result.pair_map is not None else get_pairmap_from_secstruct(ss)

    colors = None
    if seq is not None:
        try:
            from rna_draw.colorer import Colorer

            colors = Colorer().get_rgb_colors(seq, ss, default_color=None)
        except Exception:
            colors = None

    w = RnaEditor()
    w._load(
        result.x,
        result.y,
        result.node_r,
        pair_map,
        seq,
        colors,
        result.crossing_pairs,
        result.crossing_lines,
    )
    w.scene = w._build_scene(result.flagged, w._overlap_indices() if result.flagged else [])
    w.status = f"{len(ss)} nt loaded ({result.engine_name}) -- click a helix to select"
    return w


__all__ = ["RnaEditor", "editor"]
