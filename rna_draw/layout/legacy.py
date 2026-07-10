"""`LegacyEngine`: wraps the original tree-recursion layout unchanged.

Preserves today's coordinates exactly (same `DrawParameters`, same
post-shift arrays `RNARenderer.setup_tree` produces), so the M1 coordinate
baseline (`tests/test_layout_baseline.py`) still holds for this engine.
"""

from __future__ import annotations

from rna_draw.parameters import DrawParameters
from rna_draw.render_rna import RNARenderer

from .base import EngineUnavailableError, has_empty_loop


class LegacyEngine:
    """Lays out a structure with the original recursive tree algorithm."""

    name = "legacy"

    def __init__(self, params: DrawParameters | None = None) -> None:
        """Store the drawing parameters this engine lays out with.

        Args:
            params: `NODE_R`/`PRIMARY_SPACE`/`PAIR_SPACE` to use; defaults
                to `DrawParameters()`.
        """
        self._params = params or DrawParameters()

    def layout(self, secstruct: str) -> tuple[list[float], list[float]]:
        """Lay out `secstruct` via `RNARenderer.setup_tree`.

        Args:
            secstruct: Dot-bracket secondary structure.

        Returns:
            `(x, y)`, one entry per nucleotide in sequence order.

        Raises:
            EngineUnavailableError: If `secstruct` has an empty loop (a
                bare `"()"`) or otherwise crashes `setup_tree`'s recursion
                (a pre-existing `render_rna.py` `sys.exit(0)`; see the M3
                plan's Risk R7).
        """
        if len(secstruct) == 0:
            return [], []
        if has_empty_loop(secstruct):
            raise EngineUnavailableError(f"LegacyEngine cannot draw an empty loop: {secstruct!r}")

        renderer = RNARenderer()
        try:
            renderer.setup_tree(
                secstruct,
                self._params.NODE_R,
                self._params.PRIMARY_SPACE,
                self._params.PAIR_SPACE,
            )
        except SystemExit as exc:
            raise EngineUnavailableError(f"LegacyEngine failed to lay out {secstruct!r}") from exc
        return list(renderer.xarray_), list(renderer.yarray_)


__all__ = ["LegacyEngine"]
