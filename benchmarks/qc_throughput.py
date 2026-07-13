"""Batch-QC throughput benchmark: Python-loop `check_overlaps` vs the
native C++ batch overlap-checker twin.

Measures the actual throughput win the C++ overlap-checker twin
(`.claude/plans/current-plan-checker.md`) delivers at batch-QC scale: lay
out a size-stratified corpus sample once with the native puzzler engine's
parallel batch API (`NativePuzzlerEngine.layout_batch`), then time (i) a
plain Python loop of `rna_draw.overlap.check_overlaps(...).num_overlaps`
against (ii) one `check_overlaps_batch_native` call (GIL released,
`num_threads=0` -> `hardware_concurrency()`). Asserts per-structure counts
are IDENTICAL between the two before reporting wall-time + speedup -- the
throughput number is meaningful only because the two checkers already
proved to agree EXACTLY (`tests/test_overlap_native_parity.py`).

`benchmarks/pipeline_qc.py`'s own 0-silent-overlap gate is UNCHANGED and
stays on the Python `check_overlaps` -- this script only measures the
batch-QC throughput lever, never a production decision.

Usage (repo root, py3 env, `build/` already configured):
    python -m benchmarks.qc_throughput [N_PER_BUCKET] [--full]
"""

from __future__ import annotations

import sys
import time

from benchmarks.pipeline_qc import _sample
from rna_draw.layout.native import NativePuzzlerEngine
from rna_draw.overlap import OverlapParams, check_overlaps
from rna_draw.overlap_native import _native_available, check_overlaps_batch_native
from rna_draw.render_rna import get_pairmap_from_secstruct

NODE_R = 10.0
HALF_WIDTH_FACTOR = 0.75  # mirrors pipeline_qc.py's 0.75 * node_r convention

CoordsAndPairMap = tuple[list[float], list[float], list[int]]


def _lay_out_and_pair(structures: list[str]) -> list[CoordsAndPairMap]:
    """Lay out every structure with the native puzzler batch API.

    Args:
        structures: Dot-bracket secondary structures.

    Returns:
        `(x, y, pair_map)` for every structure the native engine laid out
        successfully (a pseudoknot or engine failure drops that structure
        -- this benchmark measures the checker, not the layout engine).
    """
    layouts = NativePuzzlerEngine().layout_batch(structures, num_threads=0)
    result: list[CoordsAndPairMap] = []
    for secstruct, layout in zip(structures, layouts):
        if layout is None:
            continue
        x, y = layout
        result.append((x, y, get_pairmap_from_secstruct(secstruct)))
    return result


def _python_loop_counts(cases: list[CoordsAndPairMap]) -> tuple[list[int], float]:
    """Time a plain Python loop over `rna_draw.overlap.check_overlaps`."""
    params = OverlapParams(
        node_r=NODE_R,
        backbone_half_width=HALF_WIDTH_FACTOR * NODE_R,
        pair_half_width=HALF_WIDTH_FACTOR * NODE_R,
    )
    start = time.perf_counter()
    counts = [check_overlaps(x, y, pm, params).num_overlaps for x, y, pm in cases]
    elapsed = time.perf_counter() - start
    return counts, elapsed


def _native_batch_counts(cases: list[CoordsAndPairMap]) -> tuple[list[int], float]:
    """Time one `check_overlaps_batch_native` call over the whole sample."""
    coords = [(x, y) for x, y, _pm in cases]
    pair_maps = [pm for _x, _y, pm in cases]
    node_rs = [NODE_R] * len(cases)
    start = time.perf_counter()
    counts = check_overlaps_batch_native(
        coords, pair_maps, node_rs, half_width_factor=HALF_WIDTH_FACTOR, num_threads=0
    )
    elapsed = time.perf_counter() - start
    return counts, elapsed


def _report(cases: list[CoordsAndPairMap], py_elapsed: float, cpp_elapsed: float,
           mismatches: int) -> None:
    n = len(cases)
    speedup = py_elapsed / cpp_elapsed if cpp_elapsed > 0 else float("inf")
    py_per_struct = py_elapsed / n if n else 0.0
    cpp_per_struct = cpp_elapsed / n if n else 0.0
    print("\n" + "=" * 60)
    print(f"QC THROUGHPUT: {n} structures")
    print(
        f"Python loop (check_overlaps):     {py_elapsed:.4f}s "
        f"({py_per_struct * 1e6:.1f} us/structure)"
    )
    print(
        f"C++ batch (check_overlaps_batch): {cpp_elapsed:.4f}s "
        f"({cpp_per_struct * 1e6:.1f} us/structure)"
    )
    print(f"speedup: {speedup:.1f}x")
    print(f"count mismatches: {mismatches}  (MUST be 0)")
    print(
        f"1M-structure projection: python {py_per_struct * 1_000_000 / 60:.1f} min, "
        f"cpp batch {cpp_per_struct * 1_000_000:.2f}s"
    )
    print("=" * 60)


def run(n_per_bucket: int) -> None:
    if not _native_available():
        raise SystemExit("rna_draw._layout_core is not importable; native checker unavailable")

    items = _sample(n_per_bucket)
    structures = [s for _bucket, s in items]
    print(f"Laying out {len(structures)} structures with the native puzzler batch API...")
    cases = _lay_out_and_pair(structures)
    print(f"{len(cases)}/{len(structures)} laid out successfully; checking overlaps...")

    py_counts, py_elapsed = _python_loop_counts(cases)
    cpp_counts, cpp_elapsed = _native_batch_counts(cases)
    mismatches = sum(1 for a, b in zip(py_counts, cpp_counts) if a != b)

    _report(cases, py_elapsed, cpp_elapsed, mismatches)
    assert mismatches == 0, (
        f"{mismatches} count mismatches between Python and native batch checkers"
    )


def main() -> None:
    full = "--full" in sys.argv[1:]
    positional = [a for a in sys.argv[1:] if not a.startswith("--")]
    n_per_bucket = int(positional[0]) if positional else (2000 if full else 50)
    run(n_per_bucket)


if __name__ == "__main__":
    main()
