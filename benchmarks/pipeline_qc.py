"""Large-scale QC of the WHOLE `layout_guaranteed` pipeline on the bpRNA corpus.

Unlike `hard_gate.py` (which measures one engine on the frozen 450-hard set),
this runs the real production pipeline -- production primary -> constructive
fallback -> circle fallback -- on a large size-stratified sample of the full
corpus, and verifies the honest contract at scale: every returned layout is
either checker-clean-and-unflagged (compact primary) or flagged. It counts a
SILENT OVERLAP (a `flagged=False` result that the checker finds dirty at its
own render radius) as a hard failure -- that must be 0.

One process per structure with a hard per-structure kill (the circle fallback
can be very slow on the largest structures), so a slow structure is a bounded
timeout, not a stall. Reports the tier breakdown + clean rate by size bucket.

Usage (repo root, py3 env):
    python -m benchmarks.pipeline_qc [N_PER_BUCKET] [WORKERS]
"""

from __future__ import annotations

import multiprocessing as mp
import sys
import time
from collections import Counter, defaultdict

from benchmarks.hard_gate import CORPUS, parse_dbn
from rna_draw.layout.base import is_pseudoknot_free

_BUCKETS = (
    ("<300", 0, 300),
    ("300-600", 300, 600),
    ("600-1200", 600, 1200),
    (">1200", 1200, 10**9),
)
_TIMEOUT_S = 30


def _classify(secstruct: str) -> tuple[str, int, float]:
    """Child body: run the pipeline; return (tier, num_overlaps_if_leak, seconds)."""
    from rna_draw.layout.pipeline import layout_guaranteed
    from rna_draw.overlap import OverlapParams, check_overlaps
    from rna_draw.render_rna import get_pairmap_from_secstruct

    t = time.monotonic()
    result = layout_guaranteed(secstruct)
    elapsed = time.monotonic() - t
    leak = 0
    if not result.flagged:
        pair_map = get_pairmap_from_secstruct(secstruct)
        params = OverlapParams(
            node_r=result.node_r,
            backbone_half_width=0.75 * result.node_r,
            pair_half_width=0.75 * result.node_r,
        )
        leak = check_overlaps(result.x, result.y, pair_map, params).num_overlaps
    if not result.flagged:
        tier = "primary"
    elif result.engine_name == "constructive":
        tier = "constructive"
    else:
        tier = "circle"
    return tier, leak, elapsed


def _run_one(secstruct: str, q) -> None:
    try:
        q.put(_classify(secstruct))
    except Exception as exc:  # noqa: BLE001
        q.put(("error", -1, 0.0))
        del exc


def _sample(n_per_bucket: int) -> list[tuple[str, str]]:
    """Deterministic size-stratified sample of pseudoknot-free structures."""
    filled: dict[str, list[str]] = {b[0]: [] for b in _BUCKETS}
    for path in sorted(CORPUS.glob("*.dbn"), key=lambda p: p.name):
        parsed = parse_dbn(path)
        if parsed is None:
            continue
        _seq, struct = parsed
        if not is_pseudoknot_free(struct):
            continue
        for name, lo, hi in _BUCKETS:
            if lo <= len(struct) < hi and len(filled[name]) < n_per_bucket:
                filled[name].append(struct)
                break
        if all(len(v) >= n_per_bucket for v in filled.values()):
            break
    return [(b, s) for b in filled for s in filled[b]]


def run(n_per_bucket: int, workers: int) -> None:
    items = _sample(n_per_bucket)
    ctx = mp.get_context("fork")
    by_bucket: dict[str, Counter] = defaultdict(Counter)
    times: dict[str, list[float]] = defaultdict(list)
    leaks = 0
    active: dict = {}
    nxt = done = 0

    while nxt < len(items) or active:
        while len(active) < workers and nxt < len(items):
            q = ctx.Queue()
            proc = ctx.Process(target=_run_one, args=(items[nxt][1], q))
            proc.start()
            active[proc] = (nxt, q, time.monotonic() + _TIMEOUT_S)
            nxt += 1
        for proc in list(active):
            idx, q, deadline = active[proc]
            bucket = items[idx][0]
            if proc.is_alive() and time.monotonic() <= deadline:
                continue
            if proc.is_alive():
                proc.terminate()
                proc.join()
                by_bucket[bucket]["timeout"] += 1
            else:
                proc.join()
                try:
                    tier, leak, secs = q.get_nowait()
                except Exception:  # noqa: BLE001
                    tier, leak, secs = "error", -1, 0.0
                by_bucket[bucket][tier] += 1
                if leak > 0:
                    leaks += 1
                if secs:
                    times[bucket].append(secs)
            del active[proc]
            done += 1
            if done % 200 == 0:
                print(f"  {done}/{len(items)}")
        if active:
            time.sleep(0.02)

    _report(by_bucket, times, leaks, len(items))


def _report(by_bucket, times, leaks: int, total: int) -> None:
    print("\n" + "=" * 74)
    print(f"PIPELINE QC — {total} structures, real layout_guaranteed()")
    print("-" * 74)
    print(f"{'bucket':<10}{'n':>5}{'primary%':>10}{'constr':>8}{'circle':>7}{'t/out':>7}{'p90s':>7}")
    tot = Counter()
    for name, _lo, _hi in _BUCKETS:
        c = by_bucket[name]
        n = sum(c.values())
        tot.update(c)
        if not n:
            continue
        ts = sorted(times[name])
        p90 = ts[int(0.9 * len(ts))] if ts else 0.0
        print(
            f"{name:<10}{n:>5}{100*c['primary']/n:>9.1f}%{c['constructive']:>8}"
            f"{c['circle']:>8}{c['timeout']:>9}{p90:>8.2f}"
        )
    n = sum(tot.values())
    print("-" * 74)
    print(
        f"OVERALL: primary(compact-clean) {tot['primary']}/{n} = {100*tot['primary']/n:.1f}%; "
        f"constructive-fallback {tot['constructive']}; circle {tot['circle']}; "
        f"timeout {tot['timeout']}; error {tot['error']}"
    )
    print(f"SILENT OVERLAPS (flagged=False but checker-dirty): {leaks}  (MUST be 0)")
    print("=" * 74)


def main() -> None:
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    run(n, workers)


if __name__ == "__main__":
    main()
