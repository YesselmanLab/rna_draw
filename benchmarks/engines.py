"""Pluggable engines for the hard-structure gate.

`build_engine(name)` returns an object with `layout(structure) -> (x, y)`.
Successive M5 levers are added here so each is measured by `hard_gate`
against the same frozen structure set and checker.

`EscalatingClearanceEngine`, `PostPassEngine`, and `POSTPASS_MAX_WITNESSES`
live in `rna_draw.layout.production` (the real pipeline's default engine
composition) and are re-imported here rather than duplicated, so the
benchmark and the production pipeline measure/run the exact same code (see
that module's docstring for the composition and hang-safety rationale).
"""

from __future__ import annotations

import rna_draw._vienna_layout as _vienna_layout
from rna_draw.layout.base import iter_adaptive_params
from rna_draw.layout.constructive import ConstructiveEngine
from rna_draw.layout.production import (
    POSTPASS_MAX_WITNESSES,
    EscalatingClearanceEngine,
    PostPassEngine,
)
from rna_draw.layout.vienna import (
    ViennaPuzzlerEngine,
    ViennaTurtleEngine,
)
from rna_draw.overlap import OverlapParams, check_overlaps, rescale_coords
from rna_draw.parameters import DrawParameters
from rna_draw.render_rna import get_pairmap_from_secstruct

_SIMPLE = {
    "puzzler": ViennaPuzzlerEngine,
    "turtle": ViennaTurtleEngine,
}


def _best_witnesses(x, y, pair_map) -> int:
    """Minimum overlap-witness count over the adaptive radius range.

    Delegates to `iter_adaptive_params(OverlapParams())` (the same [8, 10]
    range at 0.25 steps this used to compute inline) so the portfolio's
    engine ranking uses the one shared ladder implementation.
    """
    best = None
    for params in iter_adaptive_params(OverlapParams()):
        count = check_overlaps(x, y, pair_map, params).num_overlaps
        best = count if best is None else min(best, count)
        if best == 0:
            break
    return best if best is not None else 0


class PortfolioEngine:
    """Try each sub-engine, keep the layout with the fewest overlaps.

    The user's "different potential algorithms": puzzler/turtle fail on
    different structures, so per-structure best-of-two should beat any
    single engine. Ranks candidates by best-achievable witness count over
    the adaptive radius range; ties keep the earlier (more conventional)
    engine.
    """

    name = "portfolio"

    def __init__(self) -> None:
        self._subs = [ViennaPuzzlerEngine(), ViennaTurtleEngine()]

    def layout(self, structure: str):
        pair_map = get_pairmap_from_secstruct(structure)
        best_coords = None
        best_count = None
        for sub in self._subs:
            try:
                x, y = sub.layout(structure)
            except Exception:
                continue
            count = _best_witnesses(x, y, pair_map)
            if best_count is None or count < best_count:
                best_count, best_coords = count, (x, y)
            if best_count == 0:
                break
        if best_coords is None:
            raise RuntimeError(f"all portfolio engines failed on {structure!r}")
        return best_coords


class PuzzlerOptsEngine:
    """Puzzler with the resolver levers exposed by `plot_coords_puzzler_opts`
    (M5): `allow_flipping` and a per-call `max_config_changes` budget that
    replaces the source's hardcoded 25000 cap. Rescaled identically to
    `ViennaPuzzlerEngine` so the checker geometry matches the other engines.
    """

    def __init__(
        self, allow_flipping: bool, max_config_changes: int, clearance: float = 0.0
    ) -> None:
        self.name = (
            f"puzzler_opts(flip={int(allow_flipping)},"
            f"budget={max_config_changes},clearance={clearance})"
        )
        self._flip = allow_flipping
        self._budget = max_config_changes
        self._clearance = clearance
        self._primary = DrawParameters().PRIMARY_SPACE

    def layout(self, structure: str):
        n = len(structure)
        if n == 0:
            return [], []
        if n < 2:
            return [0.0] * n, [0.0] * n
        x, y = _vienna_layout.plot_coords_puzzler_opts(
            structure, self._flip, self._budget, self._clearance
        )
        return rescale_coords(x, y, self._primary)


def _parse_opts(name: str) -> tuple[bool, int, float]:
    """`puzzler_opts:flip,budget[,clearance]` e.g. `puzzler_opts:0,0,1.5`."""
    spec = name.split(":", 1)[1] if ":" in name else "0,0,0"
    parts = (spec.split(",") + ["0", "0", "0"])[:3]
    return bool(int(parts[0] or 0)), int(parts[1] or 0), float(parts[2] or 0.0)


__all__ = [
    "POSTPASS_MAX_WITNESSES",
    "EscalatingClearanceEngine",
    "PostPassEngine",
    "PortfolioEngine",
    "PuzzlerOptsEngine",
    "build_engine",
]


def build_engine(name: str):
    if name in _SIMPLE:
        return _SIMPLE[name]()
    if name == "constructive":
        return ConstructiveEngine()
    if name == "portfolio":
        return PortfolioEngine()
    if name == "escalating_clearance":
        return EscalatingClearanceEngine()
    if name.startswith("escalating_clearance:"):
        levels = tuple(float(x) for x in name.split(":", 1)[1].split(","))
        return EscalatingClearanceEngine(levels)
    if name == "postpass":
        return PostPassEngine()
    if name.startswith("postpass:"):
        levels = tuple(float(x) for x in name.split(":", 1)[1].split(","))
        return PostPassEngine(base=EscalatingClearanceEngine(levels))
    if name.startswith("puzzler_opts"):
        flip, budget, clearance = _parse_opts(name)
        return PuzzlerOptsEngine(flip, budget, clearance)
    raise ValueError(f"unknown gate engine: {name!r}")
