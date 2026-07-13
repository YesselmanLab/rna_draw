"""The measured-best M5 layout composition, wired for real production use.

Per-structure clearance escalation (`EscalatingClearanceEngine`) wrapped by
a monotone, time-bounded local overlap post-pass (`PostPassEngine`,
`rna_draw.layout.postpass.remove_overlaps`). On the 450-structure hard set
(`benchmarks/hard_gate.py`) this composition takes stock puzzler's ~30%
clean to ~87% clean, with 0/435 structures made worse.

This module is the SINGLE source of truth for that composition -- it used
to live only in `benchmarks/engines.py`; that module now re-imports these
names rather than keeping a parallel copy (see its module docstring).

NATIVE PRIMARY (Milestone A step 10): `EscalatingClearanceEngine`'s inner
puzzler call (`_layout_at_clearance`) now calls the owned native
`rna_layout` port (`rna_draw._layout_core.plot_coords_puzzler_full`) --
verified to reproduce the vendored `_vienna_layout.plot_coords_puzzler_opts`
call it replaces to a max diff of 3.6e-12 across the production clearance
ladder (`tests/test_native_parity.py`, Milestone A step 9). `_vienna_layout`
is only compiled under the `RNA_DRAW_BUILD_ORACLE` CMake option (Milestone
A step 11) and is used here ONLY as a fallback if `_layout_core` is somehow
unavailable (`_native_available()`/`_vienna_available()` below) -- the
default/shipped build needs no ViennaRNA header or runtime dependency.

HANG SAFETY (read before touching `PRODUCTION_CLEARANCE_LADDER`): a
diagnosis of the 17 timeouts a capped `(1.0, 1.25, 1.5)` ladder produced on
the hard set found NO uninterruptible C hang at clearance <= 1.5 -- the
puzzler C resolver's base layout finishes in ~0.1s even on >900nt
structures at that ladder. All 17 were either:

- **Empty-loop structures** (a bare `"()"` substring): puzzler hangs
  unconditionally on these, at ANY clearance. `EscalatingClearanceEngine`
  therefore raises `EngineError` on `has_empty_loop(structure)` BEFORE any
  C call (`layout_guaranteed`/`pipeline._try_primary` only guards
  pseudoknots, so without this in-process guard these structures would
  hang the real pipeline's caller). This is a required, load-bearing
  guard, not an optional optimization.
- **Large (>900nt), non-empty-loop structures**, where the SLOW part is
  entirely the Python post-pass (`remove_overlaps` grinding a large,
  many-overlap structure), not the C base layout. Python is finite and
  interruptible, so this is bounded by `PostPassConfig.time_budget_s`
  (`PostPassEngine.__init__`'s `time_budget_s` kwarg) rather than any
  further C-side change.

This is why the production engine is safe to run fully in-process (no
subprocess, no per-render process overhead): the clearance ladder is
capped at 1.5, and the only unbounded cost (the Python post-pass) is
wall-clock-bounded by construction.
"""

from __future__ import annotations

from rna_draw import _layout_core
from rna_draw.overlap import OverlapParams, check_overlaps, rescale_coords
from rna_draw.parameters import DrawParameters
from rna_draw.render_rna import get_pairmap_from_secstruct

from .base import (
    EngineError,
    LayoutEngine,
    has_empty_loop,
    is_pseudoknot_free,
    iter_adaptive_params,
)
from .postpass import POSTPASS_PARAMS, PostPassConfig, remove_overlaps
from .vienna import EXPECTED_PUZZLER_OPTIONS_SIZEOF, EXPECTED_VIENNA_ABI_VERSION

try:
    from rna_draw import _vienna_layout
except ImportError:  # pragma: no cover - exercised by the default (oracle-off) build
    _vienna_layout = None  # type: ignore[assignment]


def _native_available() -> bool:
    """Whether the native `_layout_core` puzzler call can be used.

    `_layout_core` is a default build target (compiled from this
    repository's own sources, no external ABI to drift) -- this is a
    plain import check, kept as a function (rather than inlined) so
    `_layout_at_clearance` reads as an explicit native-else-oracle
    decision and tests can monkeypatch it.
    """
    return _layout_core is not None


def _vienna_available() -> bool:
    """Whether the vendored ViennaRNA oracle engine can be used as a fallback.

    True only when `_vienna_layout` is built (`RNA_DRAW_BUILD_ORACLE=ON`)
    AND its ABI matches what `rna_draw.layout.vienna` was compiled against
    -- the same comparison `vienna._assert_vienna_abi` makes, but returning
    a bool instead of raising, so `_layout_at_clearance`'s fallback branch
    can silently decline (routing to `EngineError` instead) rather than
    crashing the caller.
    """
    if _vienna_layout is None:
        return False
    return (
        _vienna_layout.abi_version() == EXPECTED_VIENNA_ABI_VERSION
        and _vienna_layout.sizeof_puzzler_options() == EXPECTED_PUZZLER_OPTIONS_SIZEOF
    )

# Cap on how dirty a base layout may be, *measured at POSTPASS_PARAMS's floor
# radius* (see that constant's docstring), before the post-pass even attempts
# it. Structures above this are deep-dirty (rRNA-scale) territory the local
# rigid post-pass cannot realistically fix; attempting them only burns
# candidate rechecks. See `benchmarks/engines.py`'s measurement note (this
# value was tuned there: raising 6->20 lifted the hard-set clean rate
# 72.9%->87.3%).
POSTPASS_MAX_WITNESSES = 20

# Clearance ladder for the production (in-process, wall-clock-bounded)
# engine, capped at 1.5 -- see the module docstring's HANG SAFETY section.
PRODUCTION_CLEARANCE_LADDER: tuple[float, ...] = (1.0, 1.25, 1.5)

# Wall-clock budget, in seconds, for the post-pass's move loop on the
# production path (`PostPassConfig.time_budget_s`). Large non-empty-loop
# structures' post-pass is the only unbounded cost at the capped ladder
# (see module docstring); this value was measured to let those finish.
PRODUCTION_POSTPASS_BUDGET_S = 2.0


def _min_witnesses(x: list[float], y: list[float], pair_map: list[int]) -> int:
    """Minimum overlap-witness count over the adaptive radius range.

    Args:
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.

    Returns:
        The smallest `check_overlaps(...).num_overlaps` seen over
        `iter_adaptive_params(OverlapParams())`, stopping early at 0.
    """
    best = None
    for params in iter_adaptive_params(OverlapParams()):
        count = check_overlaps(x, y, pair_map, params).num_overlaps
        best = count if best is None else min(best, count)
        if best == 0:
            break
    return best if best is not None else 0


class EscalatingClearanceEngine:
    """Per-structure escalating clearance (M5).

    Lays out at the cheapest clearance first; if the checker still finds
    overlaps, retries at higher clearance (up to `levels`'s last entry),
    keeping the fewest-overlap result. Spends expensive high clearance only
    on structures that need it.
    """

    name = "escalating_clearance"

    def __init__(self, levels: tuple[float, ...] = (1.0, 1.25, 1.5, 1.75, 2.0)) -> None:
        """Store the clearance ladder.

        No ABI-drift guard is needed at construction (Milestone A step 10):
        the primary puzzler call is the native `_layout_core` extension
        (`_layout_at_clearance` below), which links no external library
        that could ABI-drift -- unlike the vendored `_vienna_layout` this
        replaced. `_layout_at_clearance` falls back to `_vienna_layout`
        (oracle build only) only if `_layout_core` is unavailable, and
        that fallback's own ABI check happens per-call
        (`_vienna_available`), not here.

        Args:
            levels: Clearance values tried in order, cheapest first.
                Defaults to the full 5-step ladder (matches the
                benchmark's measured-best setting); the production engine
                uses the hang-safe `PRODUCTION_CLEARANCE_LADDER` instead
                (see module docstring).
        """
        self._levels = levels
        self._primary = DrawParameters().PRIMARY_SPACE

    def layout(self, structure: str) -> tuple[list[float], list[float]]:
        """Lay out `structure` at the fewest-overlap clearance in the ladder.

        Args:
            structure: Dot-bracket secondary structure.

        Returns:
            `(x, y)` rescaled coordinates from the ladder level with the
            fewest overlaps (over the adaptive radius range).

        Raises:
            EngineError: If `structure` is a pseudoknot, contains an empty
                loop (puzzler hangs unconditionally on these -- see module
                docstring), or the underlying C call fails.
        """
        n = len(structure)
        if n == 0:
            return [], []
        if n < 2:
            return [0.0] * n, [0.0] * n
        if not is_pseudoknot_free(structure):
            raise EngineError(f"{self.name} cannot lay out a pseudoknot")
        if has_empty_loop(structure):
            raise EngineError(f"{self.name} would hang puzzler on an empty loop")

        pair_map = get_pairmap_from_secstruct(structure)
        return self._best_over_ladder(structure, pair_map)

    def _best_over_ladder(
        self, structure: str, pair_map: list[int]
    ) -> tuple[list[float], list[float]]:
        """Try every ladder level, keeping the fewest-overlap coordinates."""
        best_coords: tuple[list[float], list[float]] | None = None
        best_count: int | None = None
        for level in self._levels:
            x, y = self._layout_at_clearance(structure, level)
            count = _min_witnesses(x, y, pair_map)
            if best_count is None or count < best_count:
                best_count, best_coords = count, (x, y)
            if best_count == 0:
                break
        assert best_coords is not None  # `_levels` is never empty
        return best_coords

    def _layout_at_clearance(
        self, structure: str, clearance: float
    ) -> tuple[list[float], list[float]]:
        """Call the raw puzzler binding at one clearance and rescale.

        Native `_layout_core.plot_coords_puzzler_full` PRIMARY (Milestone A
        step 10; verified to reproduce this call's vendored predecessor,
        `_vienna_layout.plot_coords_puzzler_opts`, to a max diff of
        3.6e-12); falls back to the vendored oracle only if `_layout_core`
        is unavailable (`_native_available()`), and raises if neither
        extension is built.

        Raises:
            EngineError: If neither the native nor the vendored-oracle
                extension is built, or the one that ran raises
                `RuntimeError`/`ValueError` (mirrors `vienna.py`'s
                `_ViennaEngine.layout`).
        """
        try:
            if _native_available():
                rx, ry = _layout_core.plot_coords_puzzler_full(structure, False, 0, clearance)
            elif _vienna_available():
                rx, ry = _vienna_layout.plot_coords_puzzler_opts(structure, False, 0, clearance)
            else:
                raise EngineError(
                    f"{self.name}: neither rna_draw._layout_core nor rna_draw._vienna_layout "
                    "is built"
                )
        except (RuntimeError, ValueError) as exc:
            raise EngineError(f"{self.name} failed on {structure!r}") from exc
        return rescale_coords(rx, ry, self._primary)


class PostPassEngine:
    """Wrap a base engine with the rigid local post-pass (M5.2).

    Lays out with `base`, then -- only when the result is dirty and not
    too dirty to be worth it (`report_before.num_overlaps <=
    POSTPASS_MAX_WITNESSES`; see that constant's docstring) -- hands the
    coordinates to `remove_overlaps`, which is guaranteed to never increase
    the overlap count. Deep dirty structures short-circuit straight to the
    base coordinates.

    The gate/skip decision is made at `POSTPASS_PARAMS`'s floor radius, the
    same radius the pass's own acceptance checks use, because the hard gate
    scores a structure by the minimum overlap count over its adaptive
    radius range (see `rna_draw.layout.postpass`'s module docstring).
    """

    def __init__(
        self,
        base: LayoutEngine | None = None,
        name: str = "postpass",
        time_budget_s: float | None = None,
    ) -> None:
        """Store the base engine and post-pass wall-clock budget.

        Args:
            base: The engine to lay out with first; defaults to a fresh
                `EscalatingClearanceEngine()` (full ladder).
            name: Reported as `LayoutResult.engine_name`; the production
                factory (`production_engine`) passes `"production"`.
            time_budget_s: Forwarded into the `PostPassConfig` built for
                `remove_overlaps`; `None` (default) is unbounded.
        """
        self._base = base if base is not None else EscalatingClearanceEngine()
        self.name = name
        self._time_budget_s = time_budget_s

    def layout(self, structure: str) -> tuple[list[float], list[float]]:
        """Lay out `structure`, post-passing it if dirty-but-worth-it.

        Args:
            structure: Dot-bracket secondary structure.

        Returns:
            `(x, y)`: `base`'s coordinates, unchanged if already clean or
            too dirty to be worth post-passing, else the post-pass result
            (never worse than `base`'s coordinates).
        """
        x, y = self._base.layout(structure)
        if len(structure) < 2:
            return x, y
        pair_map = get_pairmap_from_secstruct(structure)
        report_before = check_overlaps(x, y, pair_map, POSTPASS_PARAMS)
        if report_before.passed or report_before.num_overlaps > POSTPASS_MAX_WITNESSES:
            return x, y
        config = PostPassConfig(time_budget_s=self._time_budget_s)
        result = remove_overlaps(x, y, pair_map, config)
        return result.x, result.y


def production_engine() -> LayoutEngine:
    """Build the production layout engine: escalation + bounded post-pass.

    Returns:
        `PostPassEngine(base=EscalatingClearanceEngine(PRODUCTION_CLEARANCE_LADDER),
        time_budget_s=PRODUCTION_POSTPASS_BUDGET_S, name="production")` --
        the hang-safe, wall-clock-bounded composition wired into
        `rna_draw.layout.pipeline.default_engine` (see module docstring).
    """
    return PostPassEngine(
        base=EscalatingClearanceEngine(PRODUCTION_CLEARANCE_LADDER),
        time_budget_s=PRODUCTION_POSTPASS_BUDGET_S,
        name="production",
    )


__all__ = [
    "POSTPASS_MAX_WITNESSES",
    "PRODUCTION_CLEARANCE_LADDER",
    "PRODUCTION_POSTPASS_BUDGET_S",
    "EscalatingClearanceEngine",
    "PostPassEngine",
    "production_engine",
]
