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

from rna_draw.layout.base import EngineError, LayoutResult, is_pseudoknot_free, params_at_node_r
from rna_draw.overlap import OverlapParams, check_overlaps
from rna_draw.render_rna import get_pairmap_from_secstruct

from .extraction import max_nested_subset
from .parsing import Stem, group_stems, nested_secstruct, parse_all_pairs, stem_pairs
from .placement import PlacementResult, place_crossings


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
    placement = place_crossings(
        nested_result.x, nested_result.y, base_pair_map, crossing, gate_params
    )

    return _assemble(nested_result, base_pair_map, placement, gate_params)


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
    nested_result: LayoutResult,
    base_pair_map: list[int],
    placement: PlacementResult,
    gate_params: OverlapParams,
) -> LayoutResult:
    """Fold crossing placement into the final pseudoknot `LayoutResult`.

    Args:
        nested_result: The nested subset's own (checker-clean-or-flagged)
            layout.
        base_pair_map: The nested subset's own pair_map.
        placement: Every crossing stem's PK-B/PK-A/unplaced outcome.
        gate_params: The SAME geometry the nested layout resolved at
            (MUST-FIX #3: gate radius == render radius).

    Returns:
        `LayoutResult(engine_name="pseudoknot", ...)` with `pair_map` =
        the FULL drawn pairs (nested + PK-B), `crossing_pairs` = the PK-B
        pairs alone (for the renderer's distinct connector color), and
        `crossing_lines` = every PK-A `RoutedLine`.
    """
    final_pair_map = list(base_pair_map)
    for i, j in placement.pk_b_pairs:
        final_pair_map[i], final_pair_map[j] = j, i
    report = check_overlaps(nested_result.x, nested_result.y, final_pair_map, gate_params)
    flagged = bool(placement.pk_a_lines) or bool(placement.unplaced) or not report.passed

    return LayoutResult(
        x=nested_result.x,
        y=nested_result.y,
        engine_name="pseudoknot",
        report=report,
        flagged=flagged,
        node_r=gate_params.node_r,
        pair_map=final_pair_map,
        crossing_pairs=placement.pk_b_pairs,
        crossing_lines=placement.pk_a_lines,
    )


__all__ = ["layout_pseudoknot"]
