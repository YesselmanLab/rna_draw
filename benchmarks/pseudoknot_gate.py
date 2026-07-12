"""Pseudoknot layout gate (M3): the corpus-scale honesty check for
`rna_draw.layout.pseudoknot.layout_pseudoknot`.

Mirrors `benchmarks/hard_gate.py`'s process-per-structure hard-kill
harness (`hard_gate.py:196`), but selects the PSEUDOKNOTTED corpus subset
(the complement of `hard_gate.py`'s own `not is_pseudoknot_free` skip) and
reports M3-specific metrics: clean %, **0 silent overlaps** (the pass/fail
gate), the % of crossing stems' pairs placed as in-plane PK-B connectors
vs routed PK-A lines vs left unplaced, mean PK-A polyline length, mean
crossing-endpoint distance (the un-biased Phase 2 objective this milestone
does NOT optimize -- Phase 2b is deferred), and readability
(`benchmarks.quality_metrics.compute_readability`).

Usage (repo root, py3 env):
    python -m benchmarks.pseudoknot_gate select [--limit 400] [--max-length 600]
    python -m benchmarks.pseudoknot_gate run [--workers 10] [--timeout 30]
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import time
from dataclasses import asdict, dataclass
from math import hypot
from pathlib import Path

from benchmarks.hard_gate import CORPUS, parse_dbn
from rna_draw.layout.base import is_pseudoknot_free

PSEUDOKNOT_SET_JSON = Path(__file__).parent / "pseudoknot_set.json"
DEFAULT_LIMIT = 400
DEFAULT_MAX_LENGTH = 600
PER_STRUCT_TIMEOUT_S = 30


def select_pseudoknot_set(limit: int, max_length: int) -> None:
    """Scan the corpus once and freeze a deterministic pseudoknot set.

    Args:
        limit: Stop once this many structures are selected.
        max_length: Skip structures longer than this (keeps a default gate
            run fast; M3's own `ConstructiveEngine` guard already caps
            nested-subset length at 4000nt independently, see
            `pseudoknot.engine._layout_nested`).
    """
    files = sorted(CORPUS.glob("*.dbn"), key=lambda p: p.name)
    selected: list[dict] = []
    for path in files:
        parsed = parse_dbn(path)
        if parsed is None:
            continue
        _seq, struct = parsed
        if is_pseudoknot_free(struct):
            continue
        if set(struct) - set(".()[]{}<>"):
            continue  # M3 scope is the four bracket types; skip bpRNA's rarer page letters
        if len(struct) > max_length:
            continue
        selected.append({"name": path.name, "structure": struct})
        if len(selected) >= limit:
            break
    PSEUDOKNOT_SET_JSON.write_text(json.dumps(selected, indent=0))
    print(f"Wrote {len(selected)} pseudoknotted structures to {PSEUDOKNOT_SET_JSON}")


@dataclass
class PkResult:
    """One structure's M3 gate measurement."""

    name: str
    length: int
    flagged: bool
    checker_clean: bool
    n_crossing_pairs: int
    n_pk_b: int
    n_pk_a: int
    n_unplaced: int
    mean_pk_a_length: float
    mean_crossing_endpoint_distance: float
    readability_index: float | None
    error: str


def _mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def _polyline_length(points: list[tuple[float, float]]) -> float:
    return sum(hypot(x1 - x0, y1 - y0) for (x0, y0), (x1, y1) in zip(points, points[1:]))


def _measure_one(name: str, struct: str) -> PkResult:
    """Lay `struct` out via `layout_pseudoknot` and compute the gate metrics."""
    from benchmarks.quality_metrics import compute_readability
    from rna_draw.layout.pseudoknot.engine import layout_pseudoknot
    from rna_draw.layout.pseudoknot.extraction import max_nested_subset
    from rna_draw.layout.pseudoknot.parsing import group_stems, parse_all_pairs
    from rna_draw.overlap import OverlapParams

    n = len(struct)
    result = layout_pseudoknot(struct, OverlapParams())

    stems = group_stems(parse_all_pairs(struct))
    _retained, crossing = max_nested_subset(stems)
    n_crossing_pairs = sum(stem.length for stem in crossing)

    pk_a_lengths = [_polyline_length(line.points) for line in result.crossing_lines]
    endpoint_dists = [
        hypot(result.x[i] - result.x[j], result.y[i] - result.y[j])
        for i, j in result.crossing_pairs
    ] + [
        hypot(result.x[line.i] - result.x[line.j], result.y[line.i] - result.y[line.j])
        for line in result.crossing_lines
    ]
    readability = compute_readability(result.x, result.y, result.pair_map, result.node_r)

    return PkResult(
        name=name,
        length=n,
        flagged=result.flagged,
        checker_clean=result.report.passed,
        n_crossing_pairs=n_crossing_pairs,
        n_pk_b=len(result.crossing_pairs),
        n_pk_a=len(result.crossing_lines),
        n_unplaced=n_crossing_pairs - len(result.crossing_pairs) - len(result.crossing_lines),
        mean_pk_a_length=_mean(pk_a_lengths),
        mean_crossing_endpoint_distance=_mean(endpoint_dists),
        readability_index=readability["readability_index"],
        error="",
    )


def _run_one_in_proc(name: str, struct: str, q: mp.Queue) -> None:
    """Child-process body: measure one structure, put the result (or an error)."""
    try:
        q.put(asdict(_measure_one(name, struct)))
    except Exception as exc:  # never let one bad structure kill the gate
        q.put({"name": name, "length": len(struct), "error": repr(exc)})


def run_gate(workers: int, timeout_s: int) -> None:
    """Measure every structure in `PSEUDOKNOT_SET_JSON` with a hard per-structure kill."""
    items = json.loads(PSEUDOKNOT_SET_JSON.read_text())
    ctx = mp.get_context("fork")
    results: list[dict] = []
    active: dict = {}
    next_i = 0

    while next_i < len(items) or active:
        while len(active) < workers and next_i < len(items):
            q = ctx.Queue()
            item = items[next_i]
            proc = ctx.Process(target=_run_one_in_proc, args=(item["name"], item["structure"], q))
            proc.start()
            active[proc] = (item, q, time.monotonic() + timeout_s)
            next_i += 1

        for proc in list(active):
            item, q, deadline = active[proc]
            finished = not proc.is_alive()
            timed_out = time.monotonic() > deadline
            if not finished and not timed_out:
                continue
            if finished:
                proc.join()
                try:
                    results.append(q.get_nowait())
                except Exception:
                    results.append(_error_result(item, "no result"))
            else:
                proc.terminate()
                proc.join()
                results.append(_error_result(item, "timeout"))
            del active[proc]
            if len(results) % 50 == 0:
                print(f"  {len(results)}/{len(items)}")

        if active:
            time.sleep(0.02)

    _report(results)


def _error_result(item: dict, reason: str) -> dict:
    """A placeholder result for a structure that errored/timed out."""
    return {"name": item["name"], "length": len(item["structure"]), "error": reason}


def _report(results: list[dict]) -> None:
    errors = [r for r in results if r.get("error")]
    ok = [r for r in results if not r.get("error")]
    clean_unflagged = sum(1 for r in ok if not r["flagged"])
    silent_overlap = [r for r in ok if not r["flagged"] and not r["checker_clean"]]

    total_pairs = sum(r["n_crossing_pairs"] for r in ok)
    total_pk_b = sum(r["n_pk_b"] for r in ok)
    total_pk_a = sum(r["n_pk_a"] for r in ok)
    total_unplaced = sum(r["n_unplaced"] for r in ok)

    print("\n" + "=" * 60)
    print(f"structures measured: {len(results)}   errors/timeouts: {len(errors)}")
    clean_pct = f"{100 * clean_unflagged / len(ok):.1f}%" if ok else "n/a"
    print(f"clean (flagged=False): {clean_unflagged}/{len(ok)} ({clean_pct})")
    print(f"SILENT OVERLAPS (flagged=False, checker-dirty): {len(silent_overlap)}  <-- MUST be 0")
    if total_pairs:
        pk_b_pct = 100 * total_pk_b / total_pairs
        pk_a_pct = 100 * total_pk_a / total_pairs
        unplaced_pct = 100 * total_unplaced / total_pairs
        print(f"crossing pairs total: {total_pairs}")
        print(f"  in-plane PK-B connectors: {total_pk_b} ({pk_b_pct:.1f}%)")
        print(f"  routed PK-A lines:        {total_pk_a} ({pk_a_pct:.1f}%)")
        print(f"  unplaced:                 {total_unplaced} ({unplaced_pct:.1f}%)")
    print("=" * 60)

    out = Path(__file__).parent / "pseudoknot_gate_result.json"
    out.write_text(json.dumps(results, indent=0))
    print(f"per-structure results -> {out}")


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("select")
    s.add_argument("--limit", type=int, default=DEFAULT_LIMIT)
    s.add_argument("--max-length", type=int, default=DEFAULT_MAX_LENGTH)
    r = sub.add_parser("run")
    r.add_argument("--workers", type=int, default=10)
    r.add_argument("--timeout", type=int, default=PER_STRUCT_TIMEOUT_S)
    args = ap.parse_args()
    if args.cmd == "select":
        select_pseudoknot_set(args.limit, args.max_length)
    else:
        run_gate(args.workers, args.timeout)


if __name__ == "__main__":
    main()
