"""Pluggable engines for the hard-structure gate.

`build_engine(name)` returns an object with `layout(structure) -> (x, y)`.
Successive M5 levers are added here so each is measured by `hard_gate`
against the same frozen structure set and checker.
"""

from __future__ import annotations

import rna_draw._vienna_layout as _vienna_layout
from rna_draw.layout.vienna import (
    ViennaNaviewEngine,
    ViennaPuzzlerEngine,
    ViennaTurtleEngine,
)
from rna_draw.overlap import OverlapParams, check_overlaps, rescale_coords
from rna_draw.parameters import DrawParameters
from rna_draw.render_rna import get_pairmap_from_secstruct

_SIMPLE = {
    "puzzler": ViennaPuzzlerEngine,
    "naview": ViennaNaviewEngine,
    "turtle": ViennaTurtleEngine,
}

# Same adaptive radius range the gate scores on, used by the portfolio to
# rank candidate engines by their best achievable witness count.
_TARGET_R = 10.0
_FLOOR_R = 8.0
_STEP_R = 0.25


def _best_witnesses(x, y, pair_map) -> int:
    floor = _FLOOR_R
    radius = _TARGET_R
    best = None
    while radius >= floor - 1e-9:
        params = OverlapParams(
            node_r=radius,
            backbone_half_width=0.75 * radius,
            pair_half_width=0.75 * radius,
        )
        count = check_overlaps(x, y, pair_map, params).num_overlaps
        best = count if best is None else min(best, count)
        if best == 0:
            break
        radius -= _STEP_R
    return best if best is not None else 0


class PortfolioEngine:
    """Try each sub-engine, keep the layout with the fewest overlaps.

    The user's "different potential algorithms": puzzler/naview/turtle fail
    on different structures, so per-structure best-of-three should beat any
    single engine. Ranks candidates by best-achievable witness count over
    the adaptive radius range; ties keep the earlier (more conventional)
    engine.
    """

    name = "portfolio"

    def __init__(self) -> None:
        self._subs = [ViennaPuzzlerEngine(), ViennaNaviewEngine(), ViennaTurtleEngine()]

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

    def __init__(self, allow_flipping: bool, max_config_changes: int) -> None:
        self.name = f"puzzler_opts(flip={int(allow_flipping)},budget={max_config_changes})"
        self._flip = allow_flipping
        self._budget = max_config_changes
        self._primary = DrawParameters().PRIMARY_SPACE

    def layout(self, structure: str):
        n = len(structure)
        if n == 0:
            return [], []
        if n < 2:
            return [0.0] * n, [0.0] * n
        x, y = _vienna_layout.plot_coords_puzzler_opts(structure, self._flip, self._budget)
        return rescale_coords(x, y, self._primary)


def _parse_opts(name: str) -> tuple[bool, int]:
    """`puzzler_opts:flip,budget` e.g. `puzzler_opts:1,1000000`."""
    spec = name.split(":", 1)[1] if ":" in name else "0,0"
    flip_s, _, budget_s = spec.partition(",")
    return bool(int(flip_s or 0)), int(budget_s or 0)


def build_engine(name: str):
    if name in _SIMPLE:
        return _SIMPLE[name]()
    if name == "portfolio":
        return PortfolioEngine()
    if name.startswith("puzzler_opts"):
        flip, budget = _parse_opts(name)
        return PuzzlerOptsEngine(flip, budget)
    raise ValueError(f"unknown gate engine: {name!r}")
