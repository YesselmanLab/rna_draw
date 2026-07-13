"""Speed-optimization gate for the native `rna_layout` puzzler core.

Companion to `.claude/plans/current-plan-speed.md`: (a) a fast, oracle-free
"golden" self-check that an EXACT-parity optimization (A1/A2/A3/B1) did not
change native output at all, and (b) per-structure + batch timing so each
step's speedup claim is measured, not asserted.

Usage (repo root, py3 env, `build/` already configured -- see
`CMakeLists.txt`):
    python -m benchmarks.bench_speed golden          # (re)write the golden JSON
    python -m benchmarks.bench_speed check            # golden self-check (fast gate)
    python -m benchmarks.bench_speed time              # per-structure staged timing
    python -m benchmarks.bench_speed batch [threads]   # batch-API throughput
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

from benchmarks.hard_gate import CORPUS, HARD_SET_JSON, parse_dbn
from rna_draw.layout.base import has_empty_loop, is_pseudoknot_free

GOLDEN_JSON = Path(__file__).parent / "bench_speed_golden.json"

# Size-stratified sample buckets, mirroring the plan's xs/s/m/l/xl stage
# breakdown (`current-plan-speed.md` section 1a): (name, lo, hi, count).
BUCKETS: tuple[tuple[str, int, int, int], ...] = (
    ("xs", 0, 80, 40),
    ("s", 80, 250, 40),
    ("m", 250, 700, 40),
    ("l", 700, 1100, 28),
    ("xl", 1100, 10**9, 15),
)

# The slow tail: the largest structures in `benchmarks/hard_set.json` (up to
# ~4000 nt), reported as their own "tail" bucket -- the plan's ">=300 nt"
# optimize-bound regime, at the extreme end.
TAIL_COUNT = 10


def _bucket_of(length: int) -> str | None:
    for name, lo, hi, _ in BUCKETS:
        if lo <= length < hi:
            return name
    return None


def collect_tail(count: int = TAIL_COUNT) -> list[tuple[str, str]]:
    """The `count` longest pseudoknot-free, empty-loop-free hard-set structures."""
    entries = json.loads(HARD_SET_JSON.read_text())
    candidates = [
        (e["name"], e["structure"])
        for e in entries
        if is_pseudoknot_free(e["structure"]) and not has_empty_loop(e["structure"])
    ]
    candidates.sort(key=lambda pair: len(pair[1]))
    return candidates[-count:]


def collect_sample() -> list[tuple[str, str]]:
    """Deterministic (name-sorted), size-stratified, pseudoknot-free sample."""
    per_bucket: dict[str, list[tuple[str, str]]] = {b[0]: [] for b in BUCKETS}
    quota = {b[0]: b[3] for b in BUCKETS}
    for path in sorted(CORPUS.glob("*.dbn"), key=lambda p: p.name):
        if all(len(per_bucket[b[0]]) >= b[3] for b in BUCKETS):
            break
        parsed = parse_dbn(path)
        if parsed is None:
            continue
        _seq, structure = parsed
        bucket = _bucket_of(len(structure))
        if bucket is None or len(per_bucket[bucket]) >= quota[bucket]:
            continue
        if not is_pseudoknot_free(structure) or has_empty_loop(structure):
            continue
        per_bucket[bucket].append((path.name, structure))
    sample: list[tuple[str, str]] = []
    for name, *_rest in BUCKETS:
        sample.extend(per_bucket[name])
    return sample


def _layout(native, structure: str) -> tuple[list[float], list[float]]:
    return native.plot_coords_puzzler_full(structure, False, 0, 1.0)


def write_golden() -> None:
    from rna_draw import _layout_core as native

    sample = collect_sample() + collect_tail()
    golden = {}
    for fname, structure in sample:
        x, y = _layout(native, structure)
        golden[fname] = {"structure": structure, "x": x, "y": y}
    GOLDEN_JSON.write_text(json.dumps(golden))
    print(f"wrote {len(golden)} structures to {GOLDEN_JSON}")


def check_golden() -> bool:
    from rna_draw import _layout_core as native

    golden = json.loads(GOLDEN_JSON.read_text())
    mismatches = 0
    for fname, entry in golden.items():
        x, y = _layout(native, entry["structure"])
        if x != entry["x"] or y != entry["y"]:
            mismatches += 1
            print(f"MISMATCH: {fname} (len={len(entry['structure'])})")
    ok = mismatches == 0
    print(f"golden check: {len(golden) - mismatches}/{len(golden)} bit-identical")
    return ok


def _percentile(values: list[float], pct: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    idx = min(len(ordered) - 1, int(round(pct * (len(ordered) - 1))))
    return ordered[idx]


def time_sample() -> None:
    from rna_draw import _layout_core as native

    labeled = [
        (fname, structure, _bucket_of(len(structure))) for fname, structure in collect_sample()
    ]
    labeled += [(fname, structure, "tail") for fname, structure in collect_tail()]

    bucket_names = [b[0] for b in BUCKETS] + ["tail"]
    by_bucket: dict[str, list[float]] = {name: [] for name in bucket_names}
    all_times: list[float] = []
    for _fname, structure, bucket in labeled:
        start = time.perf_counter()
        _layout(native, structure)
        elapsed = time.perf_counter() - start
        by_bucket[bucket].append(elapsed)
        all_times.append(elapsed)

    print(f"{'bucket':<6}{'n':>5}{'p50 (us)':>12}{'p90 (us)':>12}{'max (us)':>12}")
    for name in bucket_names:
        times_us = [t * 1e6 for t in by_bucket[name]]
        print(
            f"{name:<6}{len(times_us):>5}{_percentile(times_us, 0.5):>12.1f}"
            f"{_percentile(times_us, 0.9):>12.1f}{(max(times_us) if times_us else 0.0):>12.1f}"
        )
    all_us = [t * 1e6 for t in all_times]
    print(
        f"{'ALL':<6}{len(all_us):>5}{_percentile(all_us, 0.5):>12.1f}"
        f"{_percentile(all_us, 0.9):>12.1f}{(max(all_us) if all_us else 0.0):>12.1f}"
    )


def bench_batch(num_threads: int) -> None:
    from rna_draw import _layout_core as native

    sample = collect_sample()
    structures = [s for _fname, s in sample]

    start = time.perf_counter()
    results = native.plot_coords_puzzler_batch(structures, False, 0, 1.0, num_threads)
    elapsed = time.perf_counter() - start
    ok = sum(1 for r in results if r[2])
    rate = len(structures) / elapsed if elapsed > 0 else float("inf")
    print(
        f"batch threads={num_threads}: {len(structures)} structures in {elapsed:.4f}s "
        f"({rate:.0f} structures/s, {ok}/{len(structures)} ok); "
        f"1M projection: {1_000_000 / rate:.1f}s"
    )


def main() -> None:
    mode = sys.argv[1] if len(sys.argv) > 1 else "time"
    if mode == "golden":
        write_golden()
    elif mode == "check":
        sys.exit(0 if check_golden() else 1)
    elif mode == "time":
        time_sample()
    elif mode == "batch":
        num_threads = int(sys.argv[2]) if len(sys.argv) > 2 else 0
        bench_batch(num_threads)
    else:
        raise SystemExit(f"unknown mode: {mode}")


if __name__ == "__main__":
    main()
