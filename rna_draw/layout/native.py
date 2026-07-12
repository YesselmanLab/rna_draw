"""`LayoutEngine` for the owned native `rna_layout` core (`rna_draw._layout_core`).

`NativeTurtleEngine` calls the compiled `_layout_core` extension
(`src/layout_core/bindings.cpp`), a fresh modern-C++ port of RNAturtle's
base layout (`include/rna_layout/`, `src/layout_core/`) validated against
the vendored `_vienna_layout.plot_coords_turtle` as a parity oracle (see
`tests/test_native_parity.py`). No ViennaRNA header/runtime dependency.

SCOPE (Milestone A step 3): turtle-base only. The resolver-backed
`layout_puzzler` (and a matching `NativePuzzlerEngine`) is a later port
step; `rna_draw.layout.pipeline.resolve_engine`/`default_engine` are NOT
touched here -- this engine is importable/testable standalone but not yet
wired into the production selection chain (see `.claude/plans/
current-plan.md`, port order steps 10-11).
"""

from __future__ import annotations

from rna_draw import _layout_core
from rna_draw.overlap import rescale_coords
from rna_draw.parameters import DrawParameters

from .base import EngineError, EngineUnavailableError, is_pseudoknot_free


class NativeTurtleEngine:
    """Lays out a structure with the native `rna_layout` RNAturtle port.

    Mirrors `rna_draw.layout.vienna._ViennaEngine`'s guard/call/rescale
    shape (no ABI-drift guard is needed here: unlike `_vienna_layout`,
    `_layout_core` links no external library whose header could drift out
    from under a stale compiled extension -- it is built from this
    repository's own sources in the same CMake configure).
    """

    name = "native_turtle"

    def __init__(self, params: DrawParameters | None = None) -> None:
        """Store the drawing parameters this engine rescales output to.

        Args:
            params: Only `PRIMARY_SPACE` (target median backbone step) is
                used; defaults to `DrawParameters()`.
        """
        self._params = params or DrawParameters()

    def layout(self, secstruct: str) -> tuple[list[float], list[float]]:
        """Lay out `secstruct` with the native RNAturtle port.

        Args:
            secstruct: Dot-bracket secondary structure.

        Returns:
            `(x, y)`, one entry per nucleotide, rescaled so the median
            backbone step equals `DrawParameters.PRIMARY_SPACE`.

        Raises:
            EngineUnavailableError: If `secstruct` is a pseudoknot (this
                engine only supports well-nested structures, matching
                `rna_layout::make_pair_table`).
            EngineError: If the underlying C++ call fails (malformed
                input, e.g. an empty `"()"` loop).
        """
        n = len(secstruct)
        if n == 0:
            return [], []
        if n < 2:
            return [0.0] * n, [0.0] * n
        if not is_pseudoknot_free(secstruct):
            raise EngineUnavailableError(f"{self.name} cannot lay out a pseudoknot")

        try:
            x, y = _layout_core.plot_coords_turtle(secstruct)
        except (RuntimeError, ValueError) as exc:
            raise EngineError(f"{self.name} failed on {secstruct!r}") from exc
        return rescale_coords(x, y, self._params.PRIMARY_SPACE)


__all__ = ["NativeTurtleEngine"]
