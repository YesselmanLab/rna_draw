"""The layered "never draw a silent overlap" gate for loaded/persisted layouts.

Two checks, both reusing the FROZEN checker's own building blocks
unchanged (`overlap.py`/`geometry.py`/`pseudoknot/validate.py` are never
edited):

1. `overlap.check_overlaps` -- disks, backbone capsules, pair capsules.
2. `check_routed_lines` -- every stored PK-A `RoutedLine`, tested with the
   EXISTING `pseudoknot.validate.polyline_is_clean` (RC2): each segment
   against every base primitive, every other stored line, and its own
   non-adjacent siblings. A hand-rolled line-vs-disk-only check would miss
   line-vs-pair-capsule, line-vs-backbone, line-vs-line, and self-crossing
   overlaps -- exactly the blind spot a tampered/edited-coords document
   must not have.

`never_silent_gate` (RC1) re-validates at the pipeline's own SCALED
half-widths, `params_at_node_r(OverlapParams(node_r=target_node_r),
node_r)` -- NOT a fresh default `OverlapParams(node_r=node_r)`. The
pipeline's adaptive search (`pipeline._largest_clean_node_r`,
`pseudoknot/engine.py`'s `gate_params = params_at_node_r(params,
nested_result.node_r)`) can shrink a layout's disk radius below the target
(e.g. target 10.0 -> stored 8.0, half-widths scaling from 7.5 to 6.0); a
gate re-checking at the default `OverlapParams(node_r=8.0)` (half-widths
still 7.5) would see LESS clearance than the pipeline actually validated
and drew, and would falsely raise on a perfectly clean, shrunk layout.
"""

from __future__ import annotations

from rna_draw.layout.base import RoutedLine, params_at_node_r
from rna_draw.layout.pseudoknot.validate import polyline_capsules, polyline_is_clean
from rna_draw.overlap import (
    OverlapParams,
    OverlapReport,
    build_primitives,
    check_overlaps,
)


class LayoutOverlapError(Exception):
    """Raised when a loaded document's stored coordinates overlap.

    Distinct from the batch pipeline's `LayoutResult.flagged` contract: a
    loaded document asserts its coordinates are drawable as-is, so an
    overlap here is an error, never a silent draw.
    """


def check_routed_lines(
    x: list[float],
    y: list[float],
    pair_map: list[int],
    crossing_lines: list[RoutedLine],
    params: OverlapParams,
) -> list[RoutedLine]:
    """Which stored PK-A routed lines fail the full layered validator (RC2).

    Reuses `pseudoknot.validate.polyline_is_clean` -- the SAME validator
    `pseudoknot/engine.py`'s own `_clean_routed_lines` backstop uses -- so a
    stored line is tested against every disk/backbone/pair capsule the
    layout draws, every other stored line, and its own non-adjacent
    segments (self-crossing), not merely line-vs-disk.

    Args:
        x: Nucleotide x-coordinates (UNSHIFTED; translation-invariant).
        y: Nucleotide y-coordinates (UNSHIFTED).
        pair_map: The FULL drawn pair map (nested + PK-B), or an all-`-1`
            map if this layout has none.
        crossing_lines: Every stored PK-A `RoutedLine`.
        params: Geometry to validate at -- callers MUST already have
            scaled this to the layout's own `node_r` via `params_at_node_r`
            (see `never_silent_gate`'s docstring, RC1); this function does
            not rescale anything itself.

    Returns:
        The subset of `crossing_lines` that overlap something (empty if
        every line is clean).
    """
    if not crossing_lines:
        return []
    base_primitives = build_primitives(x, y, pair_map, params)
    committed: list = []
    dirty: list[RoutedLine] = []
    for uid, line in enumerate(crossing_lines, start=1):
        segments = polyline_capsules(line.points, line.i, line.j, params.pair_half_width, uid)
        if polyline_is_clean(segments, base_primitives, committed, pair_map, params.tol):
            committed.extend(segments)
        else:
            dirty.append(line)
    return dirty


def never_silent_gate(
    x: list[float],
    y: list[float],
    pair_map: list[int],
    crossing_lines: list[RoutedLine],
    node_r: float,
    target_node_r: float,
) -> tuple[bool, OverlapReport, list[RoutedLine]]:
    """The load/draw arbiter: frozen `check_overlaps` AND the RC2 routed-line test.

    RC1: gates at `params_at_node_r(OverlapParams(node_r=target_node_r),
    node_r)` -- the pipeline's own scaled half-widths for THIS layout's
    stored `node_r` -- reproducing exactly what `layout_guaranteed`'s
    adaptive search validated the layout against (and what the renderer
    draws at), never the library's unscaled default.

    Args:
        x: Nucleotide x-coordinates (UNSHIFTED). `check_overlaps` is
            translation-invariant (every `geometry.py` predicate uses only
            coordinate differences), so this yields the same verdict the
            pipeline got on shifted, rendered coordinates.
        y: Nucleotide y-coordinates (UNSHIFTED).
        pair_map: The FULL drawn pair map.
        crossing_lines: Every stored PK-A `RoutedLine`.
        node_r: This layout's OWN stored disk radius (`LayoutResult.node_r`
            / `LayoutBand.node_r`) -- may be smaller than `target_node_r`
            if the pipeline's adaptive search shrank it.
        target_node_r: The search ceiling this layout was originally
            produced against (the `DrawParameters.NODE_R` / preset
            `layout_defaults.node_r` in effect at layout time) -- the
            reference `OverlapParams(node_r=target_node_r)`'s default
            half-width-to-radius ratio is what `params_at_node_r` scales
            from.

    Returns:
        `(clean, report, dirty_lines)`: `clean` is True iff `report.passed`
        and `dirty_lines` is empty.
    """
    params = params_at_node_r(OverlapParams(node_r=target_node_r), node_r)
    report = check_overlaps(x, y, pair_map, params)
    dirty = check_routed_lines(x, y, pair_map, crossing_lines, params)
    return report.passed and not dirty, report, dirty


__all__ = ["LayoutOverlapError", "check_routed_lines", "never_silent_gate"]
