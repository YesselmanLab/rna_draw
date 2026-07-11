"""Pluggable engines for the hard-structure gate.

`build_engine(name)` returns an object with `layout(structure) -> (x, y)`.
Successive M5 levers are added here so each is measured by `hard_gate`
against the same frozen structure set and checker.
"""

from __future__ import annotations

import rna_draw._vienna_layout as _vienna_layout
from rna_draw.layout.postpass import POSTPASS_PARAMS, remove_overlaps
from rna_draw.layout.vienna import (
    ViennaNaviewEngine,
    ViennaPuzzlerEngine,
    ViennaTurtleEngine,
)
from rna_draw.overlap import OverlapParams, check_overlaps, rescale_coords
from rna_draw.parameters import DrawParameters
from rna_draw.render_rna import get_pairmap_from_secstruct

# Cap on how dirty a base layout may be, *measured at POSTPASS_PARAMS's floor
# radius* (see that constant's docstring), before the post-pass even attempts
# it. Deep (>6-overlap) rRNA-scale structures are out of scope (constructive
# engine territory, not local rigid moves) and would otherwise burn many
# candidate rechecks for no realistic payoff; see the post-pass plan's R1.
POSTPASS_MAX_WITNESSES = 6

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


class EscalatingClearanceEngine:
    """Per-structure escalating clearance (M5). Lay out at the cheapest
    clearance first; if the checker still finds overlaps, retry at higher
    clearance, keeping the fewest-overlap result. Spends expensive high
    clearance only on structures that need it -- most of the cost of a flat
    high clearance is on large structures that stay dirty anyway, so
    escalating cleans the many small/mid dirty structures without paying
    that cost everywhere.
    """

    name = "escalating_clearance"

    # Finer ladder (measured best on the 450 hard set: 52.2% clean vs 48.0%
    # for the 3-step ladder). The 2.0x top can send puzzler's C resolver
    # into a long loop, so this engine REQUIRES the gate's hard worker-kill
    # timeout (hard_gate.py runs one process per structure and terminate()s
    # a stuck child) -- do not use it under a plain thread/pool without one.
    def __init__(self, levels: tuple[float, ...] = (1.0, 1.25, 1.5, 1.75, 2.0)) -> None:
        self._levels = levels
        self._primary = DrawParameters().PRIMARY_SPACE

    def layout(self, structure: str):
        n = len(structure)
        if n == 0:
            return [], []
        if n < 2:
            return [0.0] * n, [0.0] * n
        pair_map = get_pairmap_from_secstruct(structure)
        best_coords = None
        best_count = None
        for level in self._levels:
            rx, ry = _vienna_layout.plot_coords_puzzler_opts(structure, False, 0, level)
            x, y = rescale_coords(rx, ry, self._primary)
            count = _best_witnesses(x, y, pair_map)
            if best_count is None or count < best_count:
                best_count, best_coords = count, (x, y)
            if best_count == 0:
                break
        assert best_coords is not None
        return best_coords


class PostPassEngine:
    """Wrap a base engine with the rigid local post-pass (M5.2).

    Lays out with `base`, then -- only when the result is dirty and not
    too dirty to be worth it (`report_before.num_overlaps <=
    POSTPASS_MAX_WITNESSES`; see that constant's docstring) -- hands the
    coordinates to `rna_draw.layout.postpass.remove_overlaps`, which is
    guaranteed to never increase the overlap count. Deep dirty structures
    (rRNA-scale, >`POSTPASS_MAX_WITNESSES` overlaps) short-circuit straight
    to the base coordinates: they are out of scope for local rigid moves
    and would only spend candidate rechecks for no realistic payoff.

    The gate/skip decision is made at `POSTPASS_PARAMS`'s floor radius, the
    same radius the pass's own acceptance checks use (`PostPassConfig`'s
    default), because the hard gate scores a structure by the minimum
    overlap count over its adaptive radius range -- checking at the
    renderer's default radius instead would see a near-clean structure as
    far dirtier than the gate does and wrongly skip it (see
    `rna_draw.layout.postpass`'s module docstring).
    """

    name = "postpass"

    def __init__(self, base=None) -> None:
        self._base = base if base is not None else EscalatingClearanceEngine()

    def layout(self, structure: str):
        x, y = self._base.layout(structure)
        if len(structure) < 2:
            return x, y
        pair_map = get_pairmap_from_secstruct(structure)
        report_before = check_overlaps(x, y, pair_map, POSTPASS_PARAMS)
        if report_before.passed or report_before.num_overlaps > POSTPASS_MAX_WITNESSES:
            return x, y
        result = remove_overlaps(x, y, pair_map)
        return result.x, result.y


def _parse_opts(name: str) -> tuple[bool, int, float]:
    """`puzzler_opts:flip,budget[,clearance]` e.g. `puzzler_opts:0,0,1.5`."""
    spec = name.split(":", 1)[1] if ":" in name else "0,0,0"
    parts = (spec.split(",") + ["0", "0", "0"])[:3]
    return bool(int(parts[0] or 0)), int(parts[1] or 0), float(parts[2] or 0.0)


def build_engine(name: str):
    if name in _SIMPLE:
        return _SIMPLE[name]()
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
