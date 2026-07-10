"""`PuzzlerEngine`: layout via ViennaRNA's puzzler (`RNAplot -t 4`).

The subprocess call and EPS parsing are isolated into module-level private
functions (`_run_rnaplot`, `_parse_coor_block`) so a future in-process
C++ binding (M5) can replace the guts without touching `PuzzlerEngine`'s
public contract.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from rna_draw.overlap import rescale_coords
from rna_draw.parameters import DrawParameters

from .base import EngineError, EngineUnavailableError, is_pseudoknot_free

RNAPLOT_BINARY = "RNAplot"
_COOR_BLOCK_RE = re.compile(r"/coor\s*\[(.*?)\]\s*def", re.S)
_COOR_POINT_RE = re.compile(r"\[\s*([-\d.eE]+)\s+([-\d.eE]+)\s*\]")


def _run_rnaplot(secstruct: str) -> str:
    """Run `RNAplot -t 4` on `secstruct` and return the generated EPS text.

    Args:
        secstruct: Dot-bracket secondary structure.

    Returns:
        The full text of the `rna.eps` file `RNAplot` writes.

    Raises:
        EngineUnavailableError: If the `RNAplot` binary is not on `PATH`.
    """
    # Puzzler layout (`-t 4`) is structure-only: nucleotide identity does not
    # affect the emitted coordinates, so any fixed-length placeholder
    # sequence folds to the same layout as the real one.
    placeholder_seq = "A" * len(secstruct)
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            subprocess.run(  # pragma: no cover -- exercised by the skip-guarded integration test
                [RNAPLOT_BINARY, "-t", "4"],
                input=f"{placeholder_seq}\n{secstruct}\n",
                capture_output=True,
                text=True,
                cwd=tmpdir,
                check=True,
            )
            return (Path(tmpdir) / "rna.eps").read_text()
    except FileNotFoundError as exc:
        raise EngineUnavailableError(f"{RNAPLOT_BINARY} not found on PATH") from exc


def _parse_coor_block(eps_text: str) -> tuple[list[float], list[float]]:
    """Parse the `/coor [...] def` block of a puzzler EPS file.

    Args:
        eps_text: Full text of an `rna.eps` file produced by `RNAplot -t 4`.

    Returns:
        `(x, y)` per-nucleotide coordinate lists, in sequence order.

    Raises:
        EngineError: If the EPS text has no `/coor` block.
    """
    match = _COOR_BLOCK_RE.search(eps_text)
    if match is None:
        raise EngineError("RNAplot EPS output missing /coor block")
    points = _COOR_POINT_RE.findall(match.group(1))
    return [float(px) for px, _ in points], [float(py) for _, py in points]


class PuzzlerEngine:
    """Lays out a structure with ViennaRNA's puzzler (`RNAplot -t 4`)."""

    name = "puzzler"

    def __init__(self, params: DrawParameters | None = None) -> None:
        """Store the drawing parameters this engine rescales output to.

        Args:
            params: Only `PRIMARY_SPACE` (target median backbone step) is
                used; defaults to `DrawParameters()`.
        """
        self._params = params or DrawParameters()

    @staticmethod
    def is_available() -> bool:
        """Check whether the `RNAplot` binary is on `PATH`.

        Returns:
            True if `RNAplot` can be invoked.
        """
        return shutil.which(RNAPLOT_BINARY) is not None

    def layout(self, secstruct: str) -> tuple[list[float], list[float]]:
        """Lay out `secstruct` with ViennaRNA's puzzler.

        Args:
            secstruct: Dot-bracket secondary structure.

        Returns:
            `(x, y)`, one entry per nucleotide, rescaled so the median
            backbone step equals `DrawParameters.PRIMARY_SPACE`.

        Raises:
            EngineUnavailableError: If `RNAplot` is absent or `secstruct`
                is a pseudoknot (puzzler cannot faithfully lay one out via
                this path).
            EngineError: If the parsed coordinate count does not match
                `len(secstruct)`.
        """
        n = len(secstruct)
        if n == 0:
            return [], []
        if n < 2:
            return [0.0] * n, [0.0] * n
        if not is_pseudoknot_free(secstruct):
            raise EngineUnavailableError("PuzzlerEngine cannot lay out a pseudoknot")

        x, y = _parse_coor_block(_run_rnaplot(secstruct))
        if len(x) != n:
            raise EngineError(
                f"RNAplot returned {len(x)} coordinates for a {n}-nucleotide structure"
            )
        return rescale_coords(x, y, self._params.PRIMARY_SPACE)


__all__ = ["PuzzlerEngine", "RNAPLOT_BINARY"]
