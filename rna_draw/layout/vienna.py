"""In-process `LayoutEngine`s for ViennaRNA's puzzler/turtle.

Each engine calls the compiled `rna_draw._vienna_layout` extension
(`src/vienna_layout/bindings.cpp`) directly -- no subprocess, no EPS
parsing -- unlike `PuzzlerEngine` (`rna_draw/layout/puzzler.py`), which is
KEPT UNCHANGED as the parity oracle these were verified against (see
`tests/test_vienna_binding.py`).

RETIRED FROM PRODUCTION (Milestone A step 11): the native `rna_layout`
port (`rna_draw.layout.native`) is now the shipped production engine
(`rna_draw.layout.production`); `_vienna_layout` is compiled only under
the `RNA_DRAW_BUILD_ORACLE` CMake option (ON for the parity/oracle build,
OFF for the default/shipped build -- see `CMakeLists.txt`), so it is no
longer a runtime dependency of the default build. The
`from rna_draw import _vienna_layout` import below is therefore GUARDED:
importing this module always succeeds (so `rna_draw.layout`'s package
import, which pulls this module in, never breaks in a default build),
but constructing `ViennaPuzzlerEngine`/`ViennaTurtleEngine` raises
`EngineUnavailableError` when the oracle extension is not built.

Both `ViennaPuzzlerEngine` and `ViennaTurtleEngine` are reentrant
(RNApuzzler/RNAturtle keep no mutable file-scope state), so they are safe
to call concurrently. The extension links no `libRNA.a`: the two vendored
layout translation units are compiled standalone against a ~120 LOC
compat shim (`src/vienna_layout/vendor/vrna_compat.c`).
"""

from __future__ import annotations

from collections.abc import Callable

from rna_draw.overlap import rescale_coords
from rna_draw.parameters import DrawParameters

from .base import EngineError, EngineUnavailableError, is_pseudoknot_free

try:
    from rna_draw import _vienna_layout
except ImportError:  # pragma: no cover - exercised by the default (oracle-off) build
    _vienna_layout = None  # type: ignore[assignment]

CoordFn = Callable[[str], tuple[list[float], list[float]]]

# ABI drift guard for the in-process ViennaRNA bindings (M5.1 Step 2): a
# header bump that reorders/extends `vrna_plot_options_puzzler_t` or
# changes ViennaRNA's version must fail loudly here, not silently
# mis-layout downstream. Values captured from the first successful build
# against ViennaRNA 2.7.0 (`micromamba run -n py3 python -c "import
# rna_draw._vienna_layout as m; print(m.abi_version(),
# m.sizeof_puzzler_options())"` printed `(2, 7, 0) 64`).
EXPECTED_VIENNA_ABI_VERSION = (2, 7, 0)
EXPECTED_PUZZLER_OPTIONS_SIZEOF = 64


def _assert_vienna_abi() -> None:
    """Fail loudly on ViennaRNA header/ABI drift before trusting a binding.

    Called from `_ViennaEngine.__init__` so *every* construction path is
    guarded, not just `resolve_engine` -- M5.2 constructs these engines
    directly to drive `vrna_plot_options_puzzler_t`, and that struct is
    exactly what a silent header bump would reorder. Compares the
    compiled-against `abi_version()`/`sizeof_puzzler_options()` to the
    values recorded when M5.1's binding was first built.

    Raises:
        EngineUnavailableError: If `_vienna_layout` was not built (the
            default, oracle-off build -- Milestone A step 11) or either
            ABI value has drifted.
    """
    if _vienna_layout is None:
        raise EngineUnavailableError(
            "rna_draw._vienna_layout is not built (RNA_DRAW_BUILD_ORACLE=OFF): "
            "the vendored ViennaRNA engines are oracle/parity-only in the "
            "default build -- use rna_draw.layout.native's NativePuzzlerEngine/"
            "NativeTurtleEngine instead, or rebuild with -DRNA_DRAW_BUILD_ORACLE=ON"
        )
    actual_version = _vienna_layout.abi_version()
    if actual_version != EXPECTED_VIENNA_ABI_VERSION:
        raise EngineUnavailableError(
            f"ViennaRNA ABI drift: compiled against {actual_version}, "
            f"expected {EXPECTED_VIENNA_ABI_VERSION}"
        )
    actual_sizeof = _vienna_layout.sizeof_puzzler_options()
    if actual_sizeof != EXPECTED_PUZZLER_OPTIONS_SIZEOF:
        raise EngineUnavailableError(
            f"ViennaRNA ABI drift: sizeof(vrna_plot_options_puzzler_t) is "
            f"{actual_sizeof}, expected {EXPECTED_PUZZLER_OPTIONS_SIZEOF}"
        )


class _ViennaEngine:
    """Shared body for the three in-process ViennaRNA engines.

    Factored out once a third near-identical engine (`ViennaTurtleEngine`)
    would otherwise duplicate `ViennaPuzzlerEngine`'s guard/call/rescale
    body a third time.
    """

    name = "vienna_base"  # overridden by each subclass

    def __init__(self, params: DrawParameters | None = None) -> None:
        """Store the drawing parameters this engine rescales output to.

        Args:
            params: Only `PRIMARY_SPACE` (target median backbone step) is
                used; defaults to `DrawParameters()`.

        Raises:
            EngineUnavailableError: If the compiled binding has drifted
                from the ViennaRNA ABI it was built against.
        """
        _assert_vienna_abi()
        self._params = params or DrawParameters()

    def _coord_fn(self) -> CoordFn:
        """The `_vienna_layout` binding function this engine calls."""
        raise NotImplementedError

    def layout(self, secstruct: str) -> tuple[list[float], list[float]]:
        """Lay out `secstruct` with this engine's ViennaRNA algorithm.

        Args:
            secstruct: Dot-bracket secondary structure.

        Returns:
            `(x, y)`, one entry per nucleotide, rescaled so the median
            backbone step equals `DrawParameters.PRIMARY_SPACE`.

        Raises:
            EngineUnavailableError: If `secstruct` is a pseudoknot or not
                well-nested (this in-process path segfaults on unbalanced
                brackets otherwise -- see `bindings.cpp`'s
                `validate_well_nested`).
            EngineError: If the underlying C call fails.
        """
        n = len(secstruct)
        if n == 0:
            return [], []
        if n < 2:
            return [0.0] * n, [0.0] * n
        if not is_pseudoknot_free(secstruct):
            raise EngineUnavailableError(f"{self.name} cannot lay out a pseudoknot")

        try:
            x, y = self._coord_fn()(secstruct)
        except (RuntimeError, ValueError) as exc:
            raise EngineError(f"{self.name} failed on {secstruct!r}") from exc
        return rescale_coords(x, y, self._params.PRIMARY_SPACE)


class ViennaPuzzlerEngine(_ViennaEngine):
    """Lays out a structure with ViennaRNA's puzzler, in-process."""

    name = "vienna_puzzler"

    def _coord_fn(self) -> CoordFn:
        return _vienna_layout.plot_coords_puzzler


class ViennaTurtleEngine(_ViennaEngine):
    """Lays out a structure with ViennaRNA's turtle, in-process."""

    name = "turtle"

    def _coord_fn(self) -> CoordFn:
        return _vienna_layout.plot_coords_turtle


__all__ = [
    "EXPECTED_PUZZLER_OPTIONS_SIZEOF",
    "EXPECTED_VIENNA_ABI_VERSION",
    "ViennaPuzzlerEngine",
    "ViennaTurtleEngine",
]
