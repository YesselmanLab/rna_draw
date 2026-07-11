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
    iter_adaptive_params,
)
from .constructive import ConstructiveEngine
from .fallback import SafeFallbackEngine
from .legacy import LegacyEngine
from .production import production_engine
from .puzzler import PuzzlerEngine
from .vienna import (
    EXPECTED_PUZZLER_OPTIONS_SIZEOF,
    EXPECTED_VIENNA_ABI_VERSION,
    ViennaNaviewEngine,
    ViennaPuzzlerEngine,
    ViennaTurtleEngine,
)

_VIENNA_ENGINE_FACTORIES: dict[str, type[LayoutEngine]] = {
    "vienna_puzzler": ViennaPuzzlerEngine,
    "naview": ViennaNaviewEngine,
    "turtle": ViennaTurtleEngine,
}


def _production_available() -> bool:
    """Whether the in-process production engine can be constructed.

    Checks that `rna_draw._vienna_layout` imports and its ABI matches what
    `rna_draw.layout.vienna` was built against -- the same comparison
    `vienna._assert_vienna_abi` makes, but returning a bool instead of
    raising, so `default_engine()` can silently fall through to the
    subprocess puzzler or legacy engine on drift instead. Kept as a
    standalone, monkeypatchable module-level function so tests can force
    either branch of `default_engine()`'s precedence.

    Returns:
        True iff the extension imports and neither `abi_version()` nor
        `sizeof_puzzler_options()` has drifted from `EXPECTED_*`.
    """
    try:
        from rna_draw import _vienna_layout
    except ImportError:
        return False
    return (
        _vienna_layout.abi_version() == EXPECTED_VIENNA_ABI_VERSION
        and _vienna_layout.sizeof_puzzler_options() == EXPECTED_PUZZLER_OPTIONS_SIZEOF
    )


def default_engine() -> LayoutEngine:
    """Pick the default engine: production, else subprocess puzzler, else legacy.

    Returns:
        `production_engine()` (in-process clearance escalation + a
        wall-clock-bounded local overlap post-pass, see
        `rna_draw.layout.production`) if `_production_available()`; else a
        `PuzzlerEngine` if `RNAplot` is on `PATH`; else a `LegacyEngine`
        (so installs without a working ViennaRNA setup keep working).
    """
    if _production_available():
        return production_engine()
    if PuzzlerEngine.is_available():
        return PuzzlerEngine()
    return LegacyEngine()


def resolve_engine(name: str) -> LayoutEngine | None:
    """Resolve a CLI/API engine-selection string to an engine instance.

    Args:
        name: One of `"auto"`, `"legacy"`, `"puzzler"`, `"production"`,
            `"constructive"`, `"vienna_puzzler"`, `"naview"`, `"turtle"`.

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
    if name == "production":
        return production_engine()
    if name == "constructive":
        return ConstructiveEngine()
    if name in _VIENNA_ENGINE_FACTORIES:
        return _VIENNA_ENGINE_FACTORIES[name]()
    raise ValueError(f"unknown layout engine: {name!r}")


def layout_guaranteed(
    secstruct: str,
    engine: LayoutEngine | None = None,
    params: OverlapParams | None = None,
) -> LayoutResult:
    """Lay out `secstruct`, guaranteeing a checker-clean or flagged result.

    Three-tier chain, each tier checker-gated (never a silent overlap):
    1. `engine` (the compact production primary by default) -- if
       checker-clean, `flagged=False`.
    2. `ConstructiveEngine` -- a compact, conventional-looking layout that
       is clean by construction and checker-verified again here; tried
       only because tier 1 was not clean, so it is reported `flagged=True`
       even though `report.passed` is also `True`.
    3. The circle `SafeFallbackEngine` -- the guaranteed-terminating last
       resort, `flagged=True`.

    Args:
        secstruct: Dot-bracket secondary structure.
        engine: Engine to try first; `None` uses `default_engine()`.
        params: Target geometry to gate against; defaults to
            `OverlapParams()`.

    Returns:
        A `LayoutResult` that is either the checker-clean primary
        (`flagged is False`) or a checker-verified fallback tier
        (`flagged is True`, `engine_name` is `"constructive"` or
        `"fallback"`) -- never a silent overlap.
    """
    params = params or OverlapParams()
    if len(secstruct) == 0:
        return LayoutResult([], [], "empty", empty_report(), flagged=False, node_r=params.node_r)

    engine = engine or default_engine()
    primary = _try_primary(engine, secstruct, params)
    if primary is not None:
        return primary
    constructive = _try_constructive_fallback(secstruct, params)
    if constructive is not None:
        return constructive
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
    for candidate in iter_adaptive_params(params):
        report = check_overlaps(x, y, pair_map, candidate)
        if report.passed:
            return candidate, report
    return None


def _try_constructive_fallback(secstruct: str, params: OverlapParams) -> LayoutResult | None:
    """Attempt the compact `ConstructiveEngine` as the middle fallback tier.

    Checker-gated exactly like `_try_primary` (same pseudoknot guard and
    `_largest_clean_node_r` search) -- the only difference is that the
    result, even though checker-clean, is honestly reported as a fallback
    (`flagged=True`): it is only ever tried because the compact production
    primary was NOT checker-clean, so it is more compact and conventional
    than the circle `SafeFallbackEngine`, but it is still not the primary.

    Args:
        secstruct: Dot-bracket secondary structure.
        params: Target geometry to search for a clean radius within.

    Returns:
        A `flagged=True` `LayoutResult` with `engine_name="constructive"`,
        or `None` if the constructive engine also can't produce a
        checker-clean layout (a pseudoknot, `EngineError` -- e.g. over its
        node-count guard -- or no clean radius in range) -- the pipeline
        should fall through to the circle `SafeFallbackEngine`.
    """
    result = _try_primary(ConstructiveEngine(), secstruct, params)
    if result is None:
        return None
    return LayoutResult(
        result.x, result.y, result.engine_name, result.report, flagged=True, node_r=result.node_r
    )


def _fallback_result(secstruct: str, params: OverlapParams) -> LayoutResult:
    """Lay out `secstruct` with the always-clean `SafeFallbackEngine` (last resort).

    Only reached once both the primary engine AND `ConstructiveEngine` have
    failed to produce a checker-clean layout (`_try_constructive_fallback`
    returned `None`) -- the final, guaranteed-terminating backstop.

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


# `EXPECTED_*` are re-exported from `.vienna` (where the ABI guard now
# lives, so it runs on every engine construction, not just the
# `resolve_engine` path). Kept importable here for back-compat with
# callers/tests using `pipeline.EXPECTED_*`.
__all__ = [
    "layout_guaranteed",
    "default_engine",
    "resolve_engine",
    "EXPECTED_VIENNA_ABI_VERSION",
    "EXPECTED_PUZZLER_OPTIONS_SIZEOF",
]
