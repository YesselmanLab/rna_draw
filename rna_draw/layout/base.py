"""The `LayoutEngine` seam: interface, shared result type, and helpers.

`LayoutEngine` is the single, stable contract every layout backend
implements -- today `LegacyEngine` (the original tree-recursion layout) and
`PuzzlerEngine` (ViennaRNA's puzzler), tomorrow an in-process C++ binding
(M5) or an axis-aligned engine (M6). Nothing engine-specific leaks into the
signature, so swapping backends never touches callers (`rna_draw/draw.py`,
`rna_draw/layout/pipeline.py`).
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from typing import Protocol, runtime_checkable

from rna_draw.overlap import OverlapKind, OverlapParams, OverlapReport

PSEUDOKNOT_BRACKETS = "()."

# Adaptive render-radius ladder shared by `pipeline._largest_clean_node_r`
# and the benchmark's `_best_witnesses` (`iter_adaptive_params`, below):
# never shrink readable disks below 80% of the target radius, searched in
# 0.25-unit steps.
MIN_NODE_R_FRACTION = 0.8
NODE_R_STEP = 0.25


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
            `"puzzler"`, `"production"`, `"constructive"`, `"fallback"`, or
            `"empty"` for a zero-length structure).
        report: The overlap report computed at `node_r`.
        flagged: True if the result is a fallback tier (the compact but
            non-primary `ConstructiveEngine`, the circle `SafeFallbackEngine`,
            or an empty/pseudoknot special case) rather than the checker-clean
            primary engine's own layout.
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


def iter_adaptive_params(params: OverlapParams) -> Iterator[OverlapParams]:
    """Yield `OverlapParams` at a descending node-radius ladder.

    The single source of truth for the adaptive-radius search duplicated
    across `pipeline._largest_clean_node_r`, the benchmark's
    `_best_witnesses`, and `hard_gate._min_witnesses`: start at
    `params.node_r` (the target) and step down by `NODE_R_STEP` to
    `MIN_NODE_R_FRACTION * params.node_r` (the floor), never shrinking
    readable disks below 80% of the target. `backbone_half_width` and
    `pair_half_width` scale with the radius, preserving their ratio to
    `params.node_r`, and `params.tol` is carried through unchanged --
    dropping `tol` would silently change the render/gate decision (see
    callers).

    Args:
        params: Target geometry; `params.node_r` is the search ceiling.

    Yields:
        `OverlapParams` at each radius in `[floor, target]`, descending
        from `target`, each carrying `params.tol` unchanged.
    """
    target = params.node_r
    floor = target * MIN_NODE_R_FRACTION
    bb_ratio = params.backbone_half_width / target if target else 0.75
    pr_ratio = params.pair_half_width / target if target else 0.75

    radius = target
    while radius >= floor - 1e-9:
        yield OverlapParams(
            node_r=radius,
            backbone_half_width=bb_ratio * radius,
            pair_half_width=pr_ratio * radius,
            tol=params.tol,
        )
        radius -= NODE_R_STEP


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
    "MIN_NODE_R_FRACTION",
    "NODE_R_STEP",
    "is_pseudoknot_free",
    "has_empty_loop",
    "empty_report",
    "iter_adaptive_params",
]
