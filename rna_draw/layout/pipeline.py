"""`layout_guaranteed`: the checker-gated layout pipeline + engine selection.

Ties every engine to the M2 checker (`check_overlaps`) so `rna_draw`'s
real output is honestly either checker-clean or explicitly flagged --
never a silent overlap. See `LayoutResult`'s docstring for the exact
contract.
"""

from __future__ import annotations

from rna_draw.overlap import OverlapParams, OverlapReport, check_overlaps
from rna_draw.render_rna import get_pairmap_from_secstruct

from .base import (
    EngineError,
    LayoutEngine,
    LayoutResult,
    empty_report,
    is_pseudoknot_free,
)
from .fallback import SafeFallbackEngine
from .legacy import LegacyEngine
from .puzzler import PuzzlerEngine
from .vienna import (
    EXPECTED_PUZZLER_OPTIONS_SIZEOF,
    EXPECTED_VIENNA_ABI_VERSION,
    ViennaNaviewEngine,
    ViennaPuzzlerEngine,
    ViennaTurtleEngine,
)

# Re-exported from `.vienna` (where the ABI guard now lives, so it runs on
# every engine construction, not just the `resolve_engine` path). Kept
# importable here for back-compat with callers/tests using
# `pipeline.EXPECTED_*`.
__all__ = ["EXPECTED_PUZZLER_OPTIONS_SIZEOF", "EXPECTED_VIENNA_ABI_VERSION"]

MIN_NODE_R_FRACTION = 0.8  # never shrink readable disks below 80% of target
NODE_R_STEP = 0.25  # search granularity, in layout units

_VIENNA_ENGINE_FACTORIES: dict[str, type[LayoutEngine]] = {
    "vienna_puzzler": ViennaPuzzlerEngine,
    "naview": ViennaNaviewEngine,
    "turtle": ViennaTurtleEngine,
}


def default_engine() -> LayoutEngine:
    """Pick the default engine: puzzler when available, else legacy.

    Returns:
        A `PuzzlerEngine` if `RNAplot` is on `PATH`, else a `LegacyEngine`
        (so installs without ViennaRNA keep working).
    """
    if PuzzlerEngine.is_available():
        return PuzzlerEngine()
    return LegacyEngine()


def resolve_engine(name: str) -> LayoutEngine | None:
    """Resolve a CLI/API engine-selection string to an engine instance.

    Args:
        name: One of `"auto"`, `"legacy"`, `"puzzler"`, `"vienna_puzzler"`,
            `"naview"`, `"turtle"`.

    Returns:
        `None` for `"auto"` (the pipeline uses `default_engine()`), or a
        fresh engine instance otherwise.

    Raises:
        ValueError: If `name` is none of the above.
        EngineUnavailableError: If `name` selects an in-process ViennaRNA
            engine and the compiled binding has drifted from the ABI this
            was built against (the guard runs in the engine's `__init__`,
            see `rna_draw.layout.vienna._assert_vienna_abi`).
    """
    if name == "auto":
        return None
    if name == "legacy":
        return LegacyEngine()
    if name == "puzzler":
        return PuzzlerEngine()
    if name in _VIENNA_ENGINE_FACTORIES:
        return _VIENNA_ENGINE_FACTORIES[name]()
    raise ValueError(f"unknown layout engine: {name!r}")


def layout_guaranteed(
    secstruct: str,
    engine: LayoutEngine | None = None,
    params: OverlapParams | None = None,
) -> LayoutResult:
    """Lay out `secstruct`, guaranteeing a checker-clean or flagged result.

    Args:
        secstruct: Dot-bracket secondary structure.
        engine: Engine to try first; `None` uses `default_engine()`.
        params: Target geometry to gate against; defaults to
            `OverlapParams()`.

    Returns:
        A `LayoutResult` that is either checker-clean (`flagged is False`)
        or the safe fallback (`flagged is True`) -- never a silent
        overlap.
    """
    params = params or OverlapParams()
    if len(secstruct) == 0:
        return LayoutResult([], [], "empty", empty_report(), flagged=False, node_r=params.node_r)

    engine = engine or default_engine()
    primary = _try_primary(engine, secstruct, params)
    if primary is not None:
        return primary
    return _fallback_result(secstruct, params)


def _try_primary(
    engine: LayoutEngine, secstruct: str, params: OverlapParams
) -> LayoutResult | None:
    """Attempt a checker-clean layout from a real engine.

    Args:
        engine: The engine to try.
        secstruct: Dot-bracket secondary structure.
        params: Target geometry to search for a clean radius within.

    Returns:
        A `flagged=False` `LayoutResult`, or `None` if the input is a
        pseudoknot, the engine raised `EngineError`, or no radius in
        `[floor, target]` is checker-clean (the pipeline should fall back).
    """
    if not is_pseudoknot_free(secstruct):
        return None
    try:
        x, y = engine.layout(secstruct)
    except EngineError:
        return None

    pair_map = get_pairmap_from_secstruct(secstruct)
    found = _largest_clean_node_r(x, y, pair_map, params)
    if found is None:
        return None
    used, report = found
    return LayoutResult(x, y, engine.name, report, flagged=False, node_r=used.node_r)


def _largest_clean_node_r(
    x: list[float], y: list[float], pair_map: list[int], params: OverlapParams
) -> tuple[OverlapParams, OverlapReport] | None:
    """Find the largest disk radius in `[floor, target]` that passes.

    Rather than gate at a single fixed `node_r` (which would dump
    otherwise-good layouts to the fallback over tiny loop-local
    near-touches), search downward from `params.node_r` and render at
    exactly the radius that clears the checker -- gate radius == render
    radius, so a clean result is never silently overlapping.

    Half-widths scale with `node_r` (their ratio to `target` is preserved)
    so the whole footprint model shrinks coherently.

    Args:
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.
        params: Target geometry; `params.node_r` is the search ceiling.

    Returns:
        `(params_at_that_radius, report)` for the largest clean radius, or
        `None` if even the floor radius fails.
    """
    target = params.node_r
    floor = target * MIN_NODE_R_FRACTION
    bb_ratio = params.backbone_half_width / target if target else 0.75
    pr_ratio = params.pair_half_width / target if target else 0.75

    radius = target
    while radius >= floor - 1e-9:
        candidate = OverlapParams(
            node_r=radius,
            backbone_half_width=bb_ratio * radius,
            pair_half_width=pr_ratio * radius,
            tol=params.tol,
        )
        report = check_overlaps(x, y, pair_map, candidate)
        if report.passed:
            return candidate, report
        radius -= NODE_R_STEP
    return None


def _fallback_result(secstruct: str, params: OverlapParams) -> LayoutResult:
    """Lay out `secstruct` with the always-clean `SafeFallbackEngine`.

    Args:
        secstruct: Dot-bracket secondary structure.
        params: Target geometry the fallback is clean by construction at.

    Returns:
        A `flagged=True` `LayoutResult`. For a pseudoknot, `pair_map`
        omits the crossing rungs, so `report` is advisory -- `flagged`
        remains the honest signal that this is not a full layout.
    """
    x, y = SafeFallbackEngine(params=params).layout(secstruct)
    pair_map = get_pairmap_from_secstruct(secstruct)
    report = check_overlaps(x, y, pair_map, params)
    return LayoutResult(x, y, "fallback", report, flagged=True, node_r=params.node_r)


__all__ = [
    "layout_guaranteed",
    "default_engine",
    "resolve_engine",
    "EXPECTED_VIENNA_ABI_VERSION",
    "EXPECTED_PUZZLER_OPTIONS_SIZEOF",
]
