"""Readability/quality harness for the real `layout_guaranteed` pipeline.

`hard_gate.py` and `pipeline_qc.py` answer "is the layout overlap-free?"; this
module answers the *next* question -- "is the overlap-free layout actually
readable, or has it sprawled into a big-rRNA hairball where every feature is
sub-pixel once fit to a page?". It runs the real production pipeline
(`layout_guaranteed`: production primary -> constructive -> circle) over a
size-stratified corpus sample and records, PER STRUCTURE:

  * length bucket, tier (primary/constructive/circle/timeout), checker pass/fail,
  * SPRAWL: bounding-box area per nucleotide, normalized by `node_r**2` so it
    is dimensionless and n-independent (a perfectly compact layout is roughly
    constant across sizes; sprawl grows it),
  * MIN FEATURE SIZE: the minimum nearest-neighbour nucleotide distance and the
    minimum helix rung (base-pair) length, each normalized by `node_r`, and
  * the headline READABILITY INDEX = min_feature_size / sqrt(bbox_area): the
    smallest feature as a fraction of the drawing's linear span. Fit the whole
    layout to a page of P pixels and the smallest feature spans
    `readability_index * P` pixels -- so a tiny index literally means an
    unreadable (sub-pixel-feature) drawing. This is the single number that
    captures "sprawl makes big rRNA visually unusable", and the yardstick an
    owned engine must beat and a resolver edit must move.

Same process-per-structure + hard-per-structure-kill pattern as `pipeline_qc.py`
(the circle fallback can be very slow on the largest structures), so a slow
structure is a bounded timeout, not a stall. Reuses `pipeline_qc._sample` /
`_BUCKETS` for the deterministic size-stratified corpus sample.

Usage (repo root, py3 env):
    python -m benchmarks.quality_metrics [N_PER_BUCKET] [WORKERS]
"""

from __future__ import annotations

import json
import multiprocessing as mp
import sys
import time
from collections import Counter, defaultdict
from math import hypot, isfinite, sqrt
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from benchmarks.pipeline_qc import _BUCKETS, _sample

_TIMEOUT_S = 30

# Metric keys aggregated per bucket, in report order. Kept in one place so the
# child, the aggregator, and the report never drift.
_METRIC_KEYS = ("area_per_nt_norm", "min_nn_norm", "min_rung_norm", "readability_index")


def compute_readability(
    x: list[float], y: list[float], pair_map: list[int], node_r: float
) -> dict:
    """Compute per-structure readability metrics for one laid-out structure.

    Pure and deterministic (no pipeline call) so it is unit-testable in
    isolation. All size metrics are normalized by `node_r` (or `node_r**2`
    for an area) so they are dimensionless and comparable across the fixed
    render geometry.

    Args:
        x: Nucleotide x-coordinates (layout units).
        y: Nucleotide y-coordinates (layout units).
        pair_map: Entry `i` holds the partner index of nucleotide `i`, or
            `-1` if unpaired.
        node_r: Disk radius the layout was gated/rendered at.

    Returns:
        A dict with `n`, `area_per_nt_norm` (bbox area / n / node_r**2),
        `min_nn_norm` (closest nucleotide-nucleotide distance / node_r),
        `min_rung_norm` (shortest base-pair rung / node_r, or `None` if
        there are no pairs), and `readability_index` (min feature size /
        sqrt(bbox area), or `None` if the layout is degenerate/collinear so
        the bbox area is 0).
    """
    n = len(x)
    ax = np.asarray(x, dtype=float)
    ay = np.asarray(y, dtype=float)
    width = float(ax.max() - ax.min()) if n else 0.0
    height = float(ay.max() - ay.min()) if n else 0.0
    area = width * height

    if n >= 2:
        pts = np.column_stack([ax, ay])
        # k=2: nearest neighbour is the 2nd result (the 1st is the point itself).
        dists, _ = cKDTree(pts).query(pts, k=2)
        min_nn = float(dists[:, 1].min())
    else:
        min_nn = 0.0

    rungs = [hypot(x[i] - x[j], y[i] - y[j]) for i, j in enumerate(pair_map) if j != -1 and i < j]
    min_rung = min(rungs) if rungs else None

    min_feature = min([min_nn] + ([min_rung] if min_rung is not None else []))
    extent = sqrt(area) if area > 0 else 0.0
    readability = min_feature / extent if extent > 0 else None

    return {
        "n": n,
        "area_per_nt_norm": (area / n) / node_r**2 if node_r > 0 else None,
        "min_nn_norm": min_nn / node_r if node_r > 0 else None,
        "min_rung_norm": (min_rung / node_r) if (min_rung is not None and node_r > 0) else None,
        "readability_index": readability,
    }


def _measure(secstruct: str) -> tuple[str, bool, dict, float]:
    """Child body: run the pipeline and measure readability of the result."""
    from rna_draw.layout.pipeline import layout_guaranteed
    from rna_draw.render_rna import get_pairmap_from_secstruct

    t = time.monotonic()
    result = layout_guaranteed(secstruct)
    elapsed = time.monotonic() - t

    if not result.flagged:
        tier = "primary"
    elif result.engine_name == "constructive":
        tier = "constructive"
    else:
        tier = "circle"

    pair_map = get_pairmap_from_secstruct(secstruct)
    metrics = compute_readability(result.x, result.y, pair_map, result.node_r)
    return tier, result.report.passed, metrics, elapsed


def _run_one(secstruct: str, q) -> None:
    try:
        q.put(_measure(secstruct))
    except Exception:  # noqa: BLE001
        q.put(("error", False, {}, 0.0))


def _percentile(values: list[float], q: float) -> float:
    """Nearest-rank percentile of a non-empty numeric list."""
    ordered = sorted(values)
    idx = min(len(ordered) - 1, int(q * len(ordered)))
    return ordered[idx]


def run(n_per_bucket: int, workers: int) -> None:
    """Measure a size-stratified corpus sample with a hard per-structure kill."""
    items = _sample(n_per_bucket)  # [(bucket_name, structure), ...]
    ctx = mp.get_context("fork")
    tiers: dict[str, Counter] = defaultdict(Counter)
    metric_vals: dict[str, dict[str, list[float]]] = defaultdict(
        lambda: {k: [] for k in _METRIC_KEYS}
    )
    times: dict[str, list[float]] = defaultdict(list)
    dirty = 0
    per_struct: list[dict] = []
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
                tiers[bucket]["timeout"] += 1
                per_struct.append({"bucket": bucket, "tier": "timeout", "passed": False})
            else:
                proc.join()
                try:
                    tier, passed, metrics, secs = q.get_nowait()
                except Exception:  # noqa: BLE001
                    tier, passed, metrics, secs = "error", False, {}, 0.0
                tiers[bucket][tier] += 1
                if tier not in ("error", "timeout"):
                    if not passed:
                        dirty += 1
                    for k in _METRIC_KEYS:
                        v = metrics.get(k)
                        if v is not None and isfinite(v):
                            metric_vals[bucket][k].append(v)
                    times[bucket].append(secs)
                    per_struct.append({"bucket": bucket, "tier": tier, "passed": passed, **metrics})
                else:
                    per_struct.append({"bucket": bucket, "tier": tier, "passed": False})
            del active[proc]
            done += 1
            if done % 200 == 0:
                print(f"  {done}/{len(items)}")
        if active:
            time.sleep(0.02)

    _report(tiers, metric_vals, times, dirty, len(items))
    out = Path(__file__).parent / "quality_result.json"
    out.write_text(json.dumps(per_struct, indent=0))
    print(f"per-structure results -> {out}")


def _report(tiers, metric_vals, times, dirty: int, total: int) -> None:
    print("\n" + "=" * 92)
    print(f"QUALITY / READABILITY — {total} structures, real layout_guaranteed()")
    print("-" * 92)
    hdr = (
        f"{'bucket':<10}{'n':>5}{'prim%':>7}"
        f"{'area/nt (med/p90)':>20}{'minNN (med/p90)':>18}{'readIdx (med/p90)':>20}"
    )
    print(hdr)

    all_read: list[float] = []
    for name, _lo, _hi in _BUCKETS:
        c = tiers[name]
        n = sum(c.values())
        if not n:
            continue
        mv = metric_vals[name]
        all_read.extend(mv["readability_index"])

        def med_p90(key: str, mv=mv) -> str:
            vals = mv[key]
            if not vals:
                return f"{'--':>9}{'--':>9}"
            return f"{_percentile(vals, 0.5):>9.3g}{_percentile(vals, 0.9):>9.3g}"

        prim = 100 * c["primary"] / n
        # readability median/p90: report p10 as the "worst 10%" is the risk tail,
        # but keep median/p90 columns consistent -- readability is higher==better,
        # so p90 here is the *good* tail and the median is the honest headline.
        print(
            f"{name:<10}{n:>5}{prim:>6.1f}%"
            f"  {med_p90('area_per_nt_norm')}  {med_p90('min_nn_norm')}  {med_p90('readability_index')}"
        )

    print("-" * 92)
    tot = Counter()
    for c in tiers.values():
        tot.update(c)
    n = sum(tot.values())
    print(
        f"TIERS: primary {tot['primary']}/{n} = {100 * tot['primary'] / n:.1f}%; "
        f"constructive {tot['constructive']}; circle {tot['circle']}; "
        f"timeout {tot['timeout']}; error {tot['error']}"
    )
    print(f"CHECKER-DIRTY (drawn but not clean at render radius): {dirty}")
    headline = _percentile(all_read, 0.5) if all_read else float("nan")
    print("-" * 92)
    print(f"HEADLINE READABILITY INDEX (overall median min_feature / sqrt(bbox_area)): {headline:.4f}")
    print("  higher == more readable; on a P-pixel page the smallest feature spans index*P px.")
    print("=" * 92)


def main() -> None:
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    run(n, workers)


if __name__ == "__main__":
    main()
