"""Portable, serializable render-intent model shared by all rna_draw frontends.

A `Scene` is the single geometry-and-style payload every future editor
frontend (the Jupyter anywidget today, a desktop/web app later) renders and
round-trips. It is deliberately dumb: pure data, JSON-serializable via
`Scene.to_dict()`, no engine or checker logic. Coordinates are stored in
ENGINE space (y-up, the same space `rna_draw.layout` and
`rna_draw.overlap` operate in); a view is responsible for flipping y to
screen space so the authoritative kernel-side overlap check always sees the
exact coordinates the engine produced.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

DEFAULT_FILL = "#9aa0a6"


def _rgb_to_hex(rgb: Sequence[float]) -> str:
    """Convert an ``(r, g, b)`` triple in ``0..1`` to a ``#rrggbb`` string.

    Args:
        rgb: Three floats in ``[0, 1]`` (extra channels such as alpha are
            ignored).

    Returns:
        A lowercase ``#rrggbb`` hex color string.
    """
    r, g, b = (max(0, min(255, round(float(c) * 255))) for c in rgb[:3])
    return f"#{r:02x}{g:02x}{b:02x}"


@dataclass
class Scene:
    """The full render intent for one RNA secondary-structure drawing.

    Args:
        nucleotides: One dict per nucleotide, index order:
            ``{"id", "x", "y", "r", "fill", "label"}``. Coordinates are in
            engine space (y-up).
        pairs: Drawn base pairs, ``{"i", "j", "kind"}`` where ``kind`` is
            ``"nested"`` or ``"crossing"``.
        routed_lines: PK-A routed connectors, ``{"i", "j", "points"}`` with
            ``points`` a list of ``[x, y]`` in engine space.
        viewport: ``{"min_x", "min_y", "w", "h"}`` bounding box (engine
            space, padded).
        node_r: Disk radius the layout was gated/rendered at.
        flagged: True if the current geometry has overlaps (or is a
            fallback tier) -- never presented as clean.
        overlaps: Nucleotide indices to tint red (the offending nts).
        pivot: ``[cx, cy]`` current rotation pivot in engine space, or
            ``None`` when nothing is selected.
    """

    nucleotides: list[dict] = field(default_factory=list)
    pairs: list[dict] = field(default_factory=list)
    routed_lines: list[dict] = field(default_factory=list)
    viewport: dict = field(default_factory=dict)
    node_r: float = 10.0
    flagged: bool = False
    overlaps: list[int] = field(default_factory=list)
    pivot: list[float] | None = None

    def to_dict(self) -> dict:
        """Serialize to a plain JSON-friendly dict (for traitlet sync)."""
        return {
            "nucleotides": self.nucleotides,
            "pairs": self.pairs,
            "routed_lines": self.routed_lines,
            "viewport": self.viewport,
            "node_r": self.node_r,
            "flagged": self.flagged,
            "overlaps": list(self.overlaps),
            "pivot": list(self.pivot) if self.pivot is not None else None,
        }


def _viewport(x: Sequence[float], y: Sequence[float], node_r: float) -> dict:
    """Padded bounding box of coordinates in engine space.

    Args:
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.
        node_r: Disk radius (drives the padding, ~3x).

    Returns:
        ``{"min_x", "min_y", "w", "h"}``; a unit box for empty input.
    """
    if not x:
        return {"min_x": 0.0, "min_y": 0.0, "w": 1.0, "h": 1.0}
    pad = 3.0 * node_r
    min_x, max_x = min(x) - pad, max(x) + pad
    min_y, max_y = min(y) - pad, max(y) + pad
    return {
        "min_x": min_x,
        "min_y": min_y,
        "w": max(max_x - min_x, 1.0),
        "h": max(max_y - min_y, 1.0),
    }


def build_scene(
    x: Sequence[float],
    y: Sequence[float],
    pair_map: Sequence[int],
    node_r: float,
    seq: str | None = None,
    colors: Sequence[Sequence[float]] | None = None,
    crossing_pairs: Sequence[tuple[int, int]] | None = None,
    routed_lines: Sequence[object] | None = None,
    flagged: bool = False,
    overlaps: Sequence[int] | None = None,
    pivot: Sequence[float] | None = None,
) -> Scene:
    """Assemble a `Scene` from raw coordinates + drawn topology.

    Shared by `scene_from_layout` (initial build) and the editor widget
    (rebuild after each committed rigid move), so both produce byte-identical
    scene shapes.

    Args:
        x: Nucleotide x-coordinates (engine space).
        y: Nucleotide y-coordinates (engine space).
        pair_map: Entry ``i`` is the partner of nucleotide ``i`` or ``-1``.
        node_r: Disk radius the layout was gated at.
        seq: Optional sequence for per-nt labels; defaults to blank labels.
        colors: Optional per-nt ``(r, g, b)`` in ``0..1``; else default gray.
        crossing_pairs: PK-B crossing pairs to mark ``kind="crossing"``.
        routed_lines: PK-A `RoutedLine`-like objects (``.i``, ``.j``,
            ``.points``).
        flagged: Whether the current geometry has overlaps.
        overlaps: Nucleotide indices to tint red.
        pivot: ``(cx, cy)`` current rotation pivot, or ``None``.

    Returns:
        A populated `Scene`.
    """
    crossing_set = {tuple(sorted(p)) for p in (crossing_pairs or [])}

    nucleotides = []
    for i in range(len(x)):
        fill = _rgb_to_hex(colors[i]) if colors is not None else DEFAULT_FILL
        label = seq[i] if seq is not None and i < len(seq) else ""
        nucleotides.append(
            {"id": i, "x": float(x[i]), "y": float(y[i]), "r": float(node_r), "fill": fill, "label": label}
        )

    pairs = []
    seen = set()
    for i, j in enumerate(pair_map):
        if j < 0 or i > j:
            continue
        key = (i, j)
        if key in seen:
            continue
        seen.add(key)
        kind = "crossing" if key in crossing_set else "nested"
        pairs.append({"i": i, "j": j, "kind": kind})
    # crossing pairs may live outside the nested pair_map view
    for i, j in crossing_set:
        if (i, j) not in seen:
            pairs.append({"i": i, "j": j, "kind": "crossing"})

    routed = []
    for line in routed_lines or []:
        routed.append(
            {"i": line.i, "j": line.j, "points": [[float(px), float(py)] for px, py in line.points]}
        )

    return Scene(
        nucleotides=nucleotides,
        pairs=pairs,
        routed_lines=routed,
        viewport=_viewport(x, y, node_r),
        node_r=float(node_r),
        flagged=bool(flagged),
        overlaps=[int(k) for k in (overlaps or [])],
        pivot=[float(pivot[0]), float(pivot[1])] if pivot is not None else None,
    )


def scene_from_layout(result, seq: str | None = None, colors=None) -> Scene:
    """Build a `Scene` from a `rna_draw.layout` `LayoutResult`.

    Args:
        result: A `LayoutResult` (``x``, ``y``, ``node_r``, ``flagged``,
            ``pair_map``, ``crossing_pairs``, ``crossing_lines``).
        seq: Optional sequence for per-nt labels.
        colors: Optional per-nt ``(r, g, b)`` in ``0..1``.

    Returns:
        A populated `Scene` in engine space.
    """
    from rna_draw.render_rna import get_pairmap_from_secstruct  # local: avoid import cost

    pair_map = result.pair_map
    if pair_map is None:
        # Reconstruct from the nucleotide count is not possible; caller-side
        # scene_from_layout is only used for the initial build where the ss is
        # known, so the editor factory passes pair_map explicitly via build_scene.
        # Fall back to an all-unpaired map if somehow absent.
        pair_map = [-1] * len(result.x)
    return build_scene(
        result.x,
        result.y,
        pair_map,
        result.node_r,
        seq=seq,
        colors=colors,
        crossing_pairs=result.crossing_pairs,
        routed_lines=result.crossing_lines,
        flagged=result.flagged,
    )


__all__ = ["Scene", "build_scene", "scene_from_layout"]
