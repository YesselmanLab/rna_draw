"""The `LayoutEngine` seam: interface, shared result type, and helpers.

`LayoutEngine` is the single, stable contract every layout backend
implements -- today `LegacyEngine` (the original tree-recursion layout) and
`PuzzlerEngine` (ViennaRNA's puzzler), tomorrow an in-process C++ binding
(M5) or an axis-aligned engine (M6). Nothing engine-specific leaks into the
signature, so swapping backends never touches callers (`rna_draw/draw.py`,
`rna_draw/layout/pipeline.py`).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

from rna_draw.overlap import OverlapKind, OverlapReport

PSEUDOKNOT_BRACKETS = "()."


@runtime_checkable
class LayoutEngine(Protocol):
    """Maps a dot-bracket structure to per-nucleotide 2D coordinates.

    Attributes:
        name: Short stable id (`"legacy"`, `"puzzler"`, `"fallback"`) used
            to report which engine produced a `LayoutResult`.
    """

    name: str

    def layout(self, secstruct: str) -> tuple[list[float], list[float]]:
        """Lay out a secondary structure.

        Args:
            secstruct: Dot-bracket secondary structure.

        Returns:
            `(x, y)`, one entry per nucleotide in sequence order, in
            rna_draw layout units with median backbone step approximately
            equal to `DrawParameters.PRIMARY_SPACE`.

        Raises:
            EngineError: If this engine cannot lay out `secstruct` (the
                caller should try another engine or the safe fallback).
        """
        ...


@dataclass
class LayoutResult:
    """The outcome of `rna_draw.layout.pipeline.layout_guaranteed`.

    Honest contract: either `report.passed and not flagged` (a real engine
    produced a checker-clean layout), or `flagged is True` (the safe
    fallback was used, or the input is a pseudoknot/degenerate case) --
    never a silent overlap.

    Args:
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.
        engine_name: Which engine actually produced `x`/`y` (`"legacy"`,
            `"puzzler"`, `"fallback"`, or `"empty"` for a zero-length
            structure).
        report: The overlap report computed at `node_r`.
        flagged: True if the result is the safe fallback (or an
            empty/pseudoknot special case) rather than a checker-clean
            primary-engine layout.
        node_r: The disk radius the report was computed at -- the renderer
            MUST draw disks at this same radius (gate radius == render
            radius), so a clean report can never describe a differently
            rendered layout.
    """

    x: list[float]
    y: list[float]
    engine_name: str
    report: OverlapReport
    flagged: bool
    node_r: float


class EngineError(Exception):
    """A layout engine failed to lay out a structure."""


class EngineUnavailableError(EngineError):
    """A layout engine's backend is not usable right now.

    Examples: the `RNAplot` binary is absent, or the input is a shape
    `LegacyEngine` cannot draw (e.g. an empty hairpin loop).
    """


def is_pseudoknot_free(secstruct: str) -> bool:
    """Check that `secstruct` uses only `().` and is well-nested.

    Args:
        secstruct: Dot-bracket secondary structure.

    Returns:
        True iff every character is in `"()."` and parentheses balance
        (never go negative, end at zero). Pseudoknot bracket types
        (`[]{}<>`) or unbalanced `()` make this False.
    """
    depth = 0
    for char in secstruct:
        if char not in PSEUDOKNOT_BRACKETS:
            return False
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
            if depth < 0:
                return False
    return depth == 0


def has_empty_loop(secstruct: str) -> bool:
    """Check whether `secstruct` contains a pair with no loop nucleotides.

    A bare `"()"` substring crashes `RNARenderer.setup_tree` (see
    `rna_draw/layout/legacy.py`); `LegacyEngine` uses this to guard against
    that pre-existing `render_rna.py` behavior.

    Args:
        secstruct: Dot-bracket secondary structure.

    Returns:
        True if `"()"` occurs anywhere in `secstruct`.
    """
    return "()" in secstruct


def empty_report() -> OverlapReport:
    """Build a zero-witness `OverlapReport` for the length-0 short-circuit.

    `check_overlaps` rejects empty coordinate input, so the pipeline's
    `n == 0` case needs a report it can hand back without calling it.

    Returns:
        An `OverlapReport` with no witnesses and every `OverlapKind` count
        at zero.
    """
    return OverlapReport(witnesses=[], counts_by_kind={kind: 0 for kind in OverlapKind})


__all__ = [
    "LayoutEngine",
    "LayoutResult",
    "EngineError",
    "EngineUnavailableError",
    "is_pseudoknot_free",
    "has_empty_loop",
    "empty_report",
]
