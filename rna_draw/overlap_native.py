"""Opt-in Python wrappers around the native C++ overlap-checker twin.

`rna_draw._layout_core`'s `check_overlaps_count`/`check_overlaps_report`/
`check_overlaps_batch` (`src/layout_core/bindings.cpp`) are a modern-C++
port of the frozen `rna_draw.overlap.check_overlaps` (+ `geometry.py` +
`spatial_hash.py`), proven to AGREE WITH IT EXACTLY by the differential
parity harness (`tests/test_overlap_native_parity.py`, `.claude/plans/
current-plan-checker.md`). `rna_draw.overlap` itself is UNTOUCHED and
remains the definitional arbiter everywhere (`benchmarks/pipeline_qc.py`'s
0-silent-overlap gate stays on the Python checker) -- this module is only
the fast BATCH-QC path (`benchmarks/qc_throughput.py`), used nowhere in
the shipped rendering pipeline.

Mirrors `rna_draw.layout.native`'s import-guard shape: a plain
`try/except ImportError`, so callers on a build without `_layout_core`
(should not happen for the default build, but mirrors the existing
precedent) get a clear `RuntimeError` instead of an `AttributeError` deep
inside a pybind11 call.
"""

from __future__ import annotations

from collections.abc import Sequence

from .overlap import OverlapParams

try:
    from rna_draw import _layout_core
except ImportError:  # pragma: no cover - exercised only on a broken build
    _layout_core = None  # type: ignore[assignment]


def _native_available() -> bool:
    """Whether the native `_layout_core` overlap-checker twin can be used."""
    return _layout_core is not None


def _require_native() -> None:
    if not _native_available():
        raise RuntimeError(
            "rna_draw._layout_core is not importable; the native overlap-checker "
            "twin is unavailable (rebuild the C++ extension, or use "
            "rna_draw.overlap.check_overlaps directly)."
        )


def check_overlaps_native(
    x: Sequence[float],
    y: Sequence[float],
    pair_map: Sequence[int],
    params: OverlapParams | None = None,
) -> int:
    """Count overlaps for one layout with the native C++ checker twin.

    Verified to agree EXACTLY with `rna_draw.overlap.check_overlaps(...).
    num_overlaps` (`tests/test_overlap_native_parity.py`) -- use this only
    where the native call's speed matters (batch QC); `rna_draw.overlap`
    remains the arbiter for any pass/fail decision.

    Args:
        x: Nucleotide x-coordinates.
        y: Nucleotide y-coordinates.
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.
        params: Geometry parameters; defaults to `OverlapParams()`.

    Returns:
        Total overlap-witness count.

    Raises:
        RuntimeError: If `_layout_core` is not importable.
        ValueError: If `x`/`y`/`pair_map` lengths mismatch, are empty, or
            `pair_map` is not symmetric (mirrors `check_overlaps`'s own
            `_validate_inputs` contract).
    """
    _require_native()
    params = params or OverlapParams()
    return _layout_core.check_overlaps_count(
        list(x),
        list(y),
        list(pair_map),
        params.node_r,
        params.backbone_half_width,
        params.pair_half_width,
        params.tol,
    )


def check_overlaps_batch_native(
    coords_list: Sequence[tuple[Sequence[float], Sequence[float]]],
    pair_maps: Sequence[Sequence[int]],
    node_rs: Sequence[float],
    half_width_factor: float = 0.75,
    tol: float = 1e-6,
    num_threads: int = 0,
) -> list[int]:
    """Count overlaps for many layouts in parallel (the QC-throughput lever).

    A GIL-released `std::thread` pool over `_layout_core.check_overlaps_batch`
    -- for QC-scale (~1M structure) corpora where looping
    `check_overlaps_native` would itself be dominated by per-call
    Python/pybind11 overhead. `backbone_half_width`/`pair_half_width` are
    derived per-structure as `half_width_factor * node_rs[i]`, reproducing
    `benchmarks/pipeline_qc.py`'s `0.75 * node_r` convention.

    Args:
        coords_list: One `(x, y)` coordinate pair per structure.
        pair_maps: One pair map per structure, same order as `coords_list`.
        node_rs: One disk radius per structure, same order.
        half_width_factor: Multiplied by each structure's own `node_r` to
            get that structure's `backbone_half_width`/`pair_half_width`.
        tol: Shared tolerance applied to every structure.
        num_threads: Worker thread count; `<= 0` uses
            `std::thread::hardware_concurrency()`.

    Returns:
        One overlap-witness count per input structure, same order; `-1`
        for a structure whose inputs are malformed (never raises for a
        single bad element).

    Raises:
        RuntimeError: If `_layout_core` is not importable.
        ValueError: If `coords_list`/`pair_maps`/`node_rs` have differing
            outer lengths (i.e. a differing number of structures) -- a
            per-structure `x`/`y`/`pair_map` length mismatch instead yields
            `-1` for that element only, per this function's own contract.
    """
    _require_native()
    xs = [list(x) for x, _ in coords_list]
    ys = [list(y) for _, y in coords_list]
    return _layout_core.check_overlaps_batch(
        xs, ys, [list(pm) for pm in pair_maps], list(node_rs), half_width_factor, tol, num_threads
    )


__all__ = ["check_overlaps_batch_native", "check_overlaps_native"]
