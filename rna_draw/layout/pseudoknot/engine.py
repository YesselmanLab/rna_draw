"""Pseudoknot layout orchestrator (M3): parse -> max-nested extraction ->
nested layout (reusing the existing checker-gated pipeline) -> crossing
placement -> assembled `LayoutResult`.

KEY GEOMETRIC DECISION (locked, see the plan): crossing nucleotides are
NEVER re-placed. Every nucleotide is placed exactly once by the nested
layout (Phase 2); a crossing pair is purely an ADDITIONAL connector -- a
straight in-plane PK-B or a routed PK-A line -- between two already-placed
nucleotides (Phase 3, `.placement`).
"""

from __future__ import annotations

from rna_draw.geometry import Capsule
from rna_draw.layout.base import (
    EngineError,
    LayoutResult,
    RoutedLine,
    is_pseudoknot_free,
    params_at_node_r,
)
from rna_draw.overlap import OverlapParams, build_primitives, check_overlaps
from rna_draw.render_rna import get_pairmap_from_secstruct

from .extraction import max_nested_subset
from .parsing import Stem, group_stems, nested_secstruct, parse_all_pairs, stem_pairs
from .placement import PlacementResult, place_crossings
from .proximity import bias_crossing_proximity
from .validate import polyline_capsules, polyline_is_clean


def layout_pseudoknot(secstruct: str, params: OverlapParams) -> LayoutResult:
    """Lay out a (possibly pseudoknotted) structure end-to-end.

    Args:
        secstruct: Dot-bracket secondary structure, any of `.()[]{}<>`.
        params: Target geometry the nested subset is laid out against.

    Returns:
        A `LayoutResult` (`engine_name="pseudoknot"`). `flagged=False`
        iff every crossing stem became a clean in-plane PK-B connector and
        the augmented `check_overlaps` passes; `flagged=True` if any
        crossing escalated to a routed PK-A line or was left unplaced --
        never a silent overlap either way.

    Raises:
        EngineError: If the max-nested extraction's own pseudoknot-free
            round-trip guard fails (an internal invariant violation, see
            `_assert_round_trip`). The caller
            (`pipeline._try_pseudoknot`) falls through to the circle
            `SafeFallbackEngine` on this.
    """
    n = len(secstruct)
    stems = group_stems(parse_all_pairs(secstruct))
    retained, crossing = max_nested_subset(stems)
    nested = nested_secstruct(n, stem_pairs(retained))
    _assert_round_trip(n, retained, nested)

    nested_result = _layout_nested(nested, params)
    gate_params = params_at_node_r(params, nested_result.node_r)
    base_pair_map = get_pairmap_from_secstruct(nested)
    bx, by = _bias_nested_layout(nested_result, base_pair_map, crossing, gate_params)
    placement = place_crossings(bx, by, base_pair_map, crossing, gate_params)

    return _assemble(bx, by, base_pair_map, placement, gate_params)


def _bias_nested_layout(
    nested_result: LayoutResult,
    base_pair_map: list[int],
    crossing: list[Stem],
    gate_params: OverlapParams,
) -> tuple[list[float], list[float]]:
    """Phase 2b entry gate: only bias an already checker-clean nested layout.

    A fallback/dirty nested layout (`not nested_result.report.passed or
    nested_result.flagged`) is left untouched -- `bias_crossing_proximity`
    itself also guards on `check_overlaps`, so this is belt-and-suspenders,
    matching the plan's "never bias a layout that isn't already clean"
    entry condition exactly.

    Args:
        nested_result: `_layout_nested`'s own result.
        base_pair_map: The nested subset's own pair_map.
        crossing: Every crossing stem.
        gate_params: The resolved geometry to gate against.

    Returns:
        `(x, y)`, possibly proximity-biased.
    """
    if not nested_result.report.passed or nested_result.flagged:
        return nested_result.x, nested_result.y
    return bias_crossing_proximity(
        nested_result.x, nested_result.y, base_pair_map, crossing, gate_params
    )


def _layout_nested(nested: str, params: OverlapParams) -> LayoutResult:
    """Lay out the pseudoknot-free nested subset with the existing pipeline.

    Local import: breaks the `pipeline` <-> `pseudoknot` module cycle
    (`pipeline._try_pseudoknot` imports `layout_pseudoknot` at module
    scope; this call only needs `layout_guaranteed` once `pipeline` has
    finished importing, which is true by the time any function here
    actually runs). `nested` is guaranteed pseudoknot-free by
    `_assert_round_trip`, so this can never recurse back into the
    pseudoknot tier (`layout_guaranteed` only calls `_try_pseudoknot` when
    `not is_pseudoknot_free(...)`).

    Args:
        nested: The `()`-only reconstructed dot-bracket string.
        params: Target geometry.

    Returns:
        `layout_guaranteed`'s result for `nested` -- checker-clean or the
        circle fallback, never a silent overlap. Reused as-is either way:
        even a circle-fallback nested layout still gets real crossing
        connectors drawn on top of it in Phase 3, a strict improvement
        over today's pair-less circle for a nested subset too large for
        the constructive tier.
    """
    from rna_draw.layout.pipeline import layout_guaranteed

    return layout_guaranteed(nested, params=params)


def _assert_round_trip(n: int, retained: list[Stem], nested: str) -> None:
    """Guard the Phase 1 -> Phase 2 handoff (non-negotiable, see the plan).

    Args:
        n: Total structure length.
        retained: The stems `max_nested_subset` retained.
        nested: `nested_secstruct(n, stem_pairs(retained))`.

    Raises:
        EngineError: If `nested` is not pseudoknot-free, or its own
            round-trip pair_map does not exactly equal `retained`'s pairs
            -- an extraction bug (MWIS's independent-set property should
            make this unreachable), but feeding a mis-extracted,
            non-nested subset to the tree engine would be a silent
            correctness break.
    """
    if not is_pseudoknot_free(nested):
        raise EngineError(
            f"pseudoknot extraction produced a non-nested subset (internal "
            f"invariant violation): {nested!r}"
        )
    expected = [-1] * n
    for i, j in stem_pairs(retained):
        expected[i], expected[j] = j, i
    if get_pairmap_from_secstruct(nested) != expected:
        raise EngineError("pseudoknot nested_secstruct round-trip mismatch (internal invariant)")


def _assemble(
    x: list[float],
    y: list[float],
    base_pair_map: list[int],
    placement: PlacementResult,
    gate_params: OverlapParams,
) -> LayoutResult:
    """Fold crossing placement into the final pseudoknot `LayoutResult`.

    Args:
        x: The nucleotide x-coordinates the layout was actually drawn/
            checked at -- the nested layout's own coordinates, or Phase
            2b's proximity-biased coordinates when biasing ran.
        y: The corresponding y-coordinates.
        base_pair_map: The nested subset's own pair_map.
        placement: Every crossing stem's PK-B/PK-A/unplaced outcome
            (already computed against `x`/`y`).
        gate_params: The SAME geometry the nested layout resolved at
            (MUST-FIX #3: gate radius == render radius).

    Returns:
        `LayoutResult(engine_name="pseudoknot", ...)` with `pair_map` =
        the FULL drawn pairs (nested + PK-B), `crossing_pairs` = the PK-B
        pairs alone (for the renderer's distinct connector color), and
        `crossing_lines` = every PK-A `RoutedLine` that survived the final
        mutual-validation backstop (see `_clean_routed_lines`).
    """
    final_pair_map = list(base_pair_map)
    for i, j in placement.pk_b_pairs:
        final_pair_map[i], final_pair_map[j] = j, i
    report = check_overlaps(x, y, final_pair_map, gate_params)
    clean_lines, dropped_any = _clean_routed_lines(
        x, y, final_pair_map, placement.pk_a_lines, gate_params
    )
    flagged = bool(clean_lines) or bool(placement.unplaced) or dropped_any or not report.passed

    return LayoutResult(
        x=x,
        y=y,
        engine_name="pseudoknot",
        report=report,
        flagged=flagged,
        node_r=gate_params.node_r,
        pair_map=final_pair_map,
        crossing_pairs=placement.pk_b_pairs,
        crossing_lines=clean_lines,
    )


def _clean_routed_lines(
    x: list[float],
    y: list[float],
    final_pair_map: list[int],
    lines: list[RoutedLine],
    params: OverlapParams,
) -> tuple[list[RoutedLine], bool]:
    """Belt-and-suspenders backstop (BUG 1 fix): drop any routed line the
    incremental placement logic somehow let through with an overlap.

    Re-validates every PK-A line's own capsules against `final_pair_map`'s
    base primitives (nested pairs AND PK-B rungs) AND every other routed
    line, using the same frozen predicates `placement`/`validate` use.
    Incremental validation (`placement._try_in_plane`/`_route_stem`)
    should make this always pass; if it ever doesn't, that is an
    incremental-logic bug -- the offending line is simply never drawn
    rather than silently overlapping.

    Args:
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.
        final_pair_map: The nested + PK-B pair_map (`build_primitives`'s
            input covers disks, backbone, and every pair capsule).
        lines: Every PK-A `RoutedLine` the placement phase committed.
        params: Geometry to validate against (the render radius).

    Returns:
        `(clean_lines, dropped_any)`: the subset of `lines` that are
        mutually clean against everything, and whether any line had to be
        dropped.
    """
    base_primitives = build_primitives(x, y, final_pair_map, params)
    committed: list[Capsule] = []
    clean: list[RoutedLine] = []
    dropped_any = False
    for uid, line in enumerate(lines, start=1):
        segments = polyline_capsules(line.points, line.i, line.j, params.pair_half_width, uid)
        if polyline_is_clean(segments, base_primitives, committed, final_pair_map, params.tol):
            clean.append(line)
            committed.extend(segments)
        else:
            dropped_any = True
    return clean, dropped_any


__all__ = ["layout_pseudoknot"]
