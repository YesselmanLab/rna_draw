"""Hard-structure overlap gate: the number we drive down.

M5 goal is a layout engine that is overlap-free on *all* non-pseudoknotted
RNAs. Stock puzzler is ~90% clean overall but collapses on large rRNAs
(bpRNA: ~46% clean at 300-600 nt, ~8% above 1200 nt). This gate fixes a
deterministic set of the *hard* (large) structures, lays each out, and
counts the total residual overlap witnesses across the set at the best
radius in the adaptive range. That single total is the optimization target:
lower is better, 0 is the goal.

The engine is pluggable (`--engine`) so successive levers -- puzzler option
sweeps, the turtle portfolio, and eventually an altered resolver -- are all
measured against the *same* frozen structure set and the *same* frozen M2
checker (`check_overlaps`), so the number is comparable across iterations.

PARALLELISM: the gate runs one process per structure (never threads). This
bounds turtle's known upstream per-call memory leak via worker recycling
and lets a stuck C resolver be hard-killed on timeout. See
`rna_draw/layout/vienna.py`.

Usage:
    python -m benchmarks.hard_gate select      # build the frozen hard set (once)
    python -m benchmarks.hard_gate run --engine puzzler [--workers 10]
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import time
from dataclasses import dataclass
from pathlib import Path

from rna_draw.layout.base import is_pseudoknot_free
from rna_draw.overlap import OverlapParams, check_overlaps
from rna_draw.render_rna import get_pairmap_from_secstruct

CORPUS = Path("/Users/jyesselman2/Downloads/dbnFiles")
HARD_SET_JSON = Path(__file__).parent / "hard_set.json"
WORST_SET_JSON = Path(__file__).parent / "worst_set.json"

# Size buckets (nt) and how many structures to sample from each. Deterministic
# even-spaced sampling over the name-sorted bucket, so the set is reproducible
# and stable across engine iterations.
BUCKETS: tuple[tuple[str, int, int, int], ...] = (
    ("300-600", 300, 600, 200),
    ("600-1200", 600, 1200, 150),
    ("1200-4000", 1200, 4000, 100),
)

# Adaptive radius search, mirroring pipeline._largest_clean_node_r: render at
# the largest radius that is clean, but for the metric we take the *minimum*
# witness count over the range (best achievable in-conventional-look before
# any fallback), so a structure that never goes clean still contributes its
# smallest honest overlap count rather than being masked by the safe circle.
TARGET_NODE_R = 10.0
MIN_NODE_R_FRACTION = 0.8
NODE_R_STEP = 0.25
PER_STRUCT_TIMEOUT_S = 30


def parse_dbn(path: Path) -> tuple[str, str] | None:
    """Return `(sequence, structure)` from a bpRNA `.dbn`, or None if unusable.

    Format: some `#`-comment lines, then the sequence, then the dot-bracket
    structure. The structure is the first non-comment line made solely of
    dot-bracket characters.
    """
    seq = None
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if seq is None:
            seq = line
            continue
        if set(line) <= set(".()[]{}<>anAB"):
            return seq, line
    return None


def _bucket_of(length: int) -> str | None:
    for name, lo, hi, _ in BUCKETS:
        if lo <= length < hi:
            return name
    return None


def select_hard_set() -> None:
    """Scan the corpus once and freeze a deterministic hard-structure set."""
    per_bucket: dict[str, list[tuple[str, str]]] = {b[0]: [] for b in BUCKETS}
    files = sorted(CORPUS.glob("*.dbn"), key=lambda p: p.name)
    scanned = 0
    for path in files:
        scanned += 1
        if scanned % 10000 == 0:
            print(f"  scanned {scanned}/{len(files)}")
        parsed = parse_dbn(path)
        if parsed is None:
            continue
        _seq, struct = parsed
        bucket = _bucket_of(len(struct))
        if bucket is None:
            continue
        if not is_pseudoknot_free(struct):
            continue
        per_bucket[bucket].append((path.name, struct))

    selected: list[dict] = []
    for name, _lo, _hi, want in BUCKETS:
        pool = per_bucket[name]
        if not pool:
            continue
        if len(pool) <= want:
            picks = pool
        else:
            step = len(pool) / want
            picks = [pool[int(i * step)] for i in range(want)]
        for fname, struct in picks:
            selected.append({"name": fname, "bucket": name, "structure": struct})
        print(f"  bucket {name}: {len(pool)} available -> {len(picks)} selected")

    HARD_SET_JSON.write_text(json.dumps(selected, indent=0))
    print(f"Wrote {len(selected)} structures to {HARD_SET_JSON}")


@dataclass
class StructResult:
    name: str
    bucket: str
    length: int
    witnesses: int  # min over the adaptive radius range; -1 == engine error
    best_radius: float
    engine: str


def _min_witnesses(x, y, pair_map) -> tuple[int, float]:
    """Minimum overlap-witness count over the adaptive radius range."""
    floor = TARGET_NODE_R * MIN_NODE_R_FRACTION
    best_count = None
    best_radius = TARGET_NODE_R
    radius = TARGET_NODE_R
    while radius >= floor - 1e-9:
        params = OverlapParams(
            node_r=radius,
            backbone_half_width=0.75 * radius,
            pair_half_width=0.75 * radius,
        )
        count = check_overlaps(x, y, pair_map, params).num_overlaps
        if best_count is None or count < best_count:
            best_count, best_radius = count, radius
        if best_count == 0:
            break
        radius -= NODE_R_STEP
    return (best_count if best_count is not None else 0), best_radius


def _run_one_in_proc(engine_name: str, struct: str, q) -> None:
    """Child-process body: lay out one structure, put (witnesses, radius).

    Runs in its own process so the parent can HARD-KILL it on timeout --
    puzzler's C resolver can enter a long loop that a Python SIGALRM cannot
    interrupt (the signal only fires between C calls), so the previous
    in-worker alarm approach let one structure stall the whole run. Here a
    stuck child is `terminate()`d by the parent and recorded as an error.
    """
    from benchmarks.engines import build_engine

    try:
        engine = build_engine(engine_name)
        x, y = engine.layout(struct)
        witnesses, best_radius = _min_witnesses(x, y, get_pairmap_from_secstruct(struct))
        q.put((witnesses, best_radius))
    except Exception:
        q.put((-1, 0.0))


def build_worst_set(from_result: Path, n: int) -> None:
    """Freeze the N worst (most-overlapping) structures for fast iteration.

    Joins a prior `gate_result_*.json` (witness counts) back to the full
    hard set (structures) and writes the top-N by witness count. Errors
    (witnesses < 0) sort last -- they are failures a lever can't measure by
    overlap count, kept only if fewer than N dirty structures exist.
    """
    by_name = {s["name"]: s for s in json.loads(HARD_SET_JSON.read_text())}
    results = json.loads(from_result.read_text())
    ranked = sorted(results, key=lambda r: (r["witnesses"] < 0, -r["witnesses"]))
    worst = [by_name[r["name"]] for r in ranked[:n] if r["name"] in by_name]
    WORST_SET_JSON.write_text(json.dumps(worst, indent=0))
    total = sum(max(r["witnesses"], 0) for r in ranked[:n])
    print(f"Wrote {len(worst)} worst structures ({total} overlaps) to {WORST_SET_JSON}")


def run_gate(engine_name: str, workers: int, set_path: Path) -> None:
    """Measure every structure in `set_path` with a HARD per-structure kill.

    One process per structure, at most `workers` concurrent. A child still
    running after `PER_STRUCT_TIMEOUT_S` is terminated and recorded as an
    error -- the only way to bound puzzler's uninterruptible C hangs at high
    clearance. Process-based (never threads) also caps turtle's known
    upstream per-call memory leak via worker recycling.
    """
    items = json.loads(set_path.read_text())
    ctx = mp.get_context("fork")
    results: list[StructResult | None] = [None] * len(items)
    active: dict = {}  # proc -> (index, queue, deadline)
    next_i = 0
    done = 0

    while next_i < len(items) or active:
        while len(active) < workers and next_i < len(items):
            q = ctx.Queue()
            proc = ctx.Process(
                target=_run_one_in_proc, args=(engine_name, items[next_i]["structure"], q)
            )
            proc.start()
            active[proc] = (next_i, q, time.monotonic() + PER_STRUCT_TIMEOUT_S)
            next_i += 1

        for proc in list(active):
            idx, q, deadline = active[proc]
            item = items[idx]
            n = len(item["structure"])
            finished = not proc.is_alive()
            timed_out = time.monotonic() > deadline
            if not finished and not timed_out:
                continue
            if finished:
                proc.join()
                try:
                    witnesses, best_radius = q.get_nowait()
                except Exception:
                    witnesses, best_radius = -1, 0.0
            else:  # hard-kill the stuck child
                proc.terminate()
                proc.join()
                witnesses, best_radius = -1, 0.0
            results[idx] = StructResult(
                item["name"], item["bucket"], n, witnesses, best_radius, engine_name
            )
            del active[proc]
            done += 1
            if done % 50 == 0:
                print(f"  {done}/{len(items)}")

        if active:
            time.sleep(0.02)

    _report(engine_name, [r for r in results if r is not None])


def _report(engine_name: str, results: list[StructResult]) -> None:
    errors = [r for r in results if r.witnesses < 0]
    ok = [r for r in results if r.witnesses >= 0]
    total = sum(r.witnesses for r in ok)
    clean = sum(1 for r in ok if r.witnesses == 0)
    dirty = sum(1 for r in ok if r.witnesses > 0)

    # Primary metric = NON-CLEAN structures (dirty + errors). Robust to the
    # error-handling artifact where a structure an engine cannot lay out
    # silently drops from the overlap sum; a failure is not an improvement.
    non_clean = dirty + len(errors)
    print("\n" + "=" * 60)
    print(f"ENGINE: {engine_name}")
    print(f"structures: {len(results)}")
    print(f"NON-CLEAN (dirty+errors): {non_clean}   "
          f"[clean {clean} / {len(results)} = {100 * clean / len(results):.1f}%]")
    print(f"  dirty(drawn w/ overlaps): {dirty}   errors/timeouts: {len(errors)}")
    print(f"TOTAL OVERLAPS (over drawn structures only): {total}")
    print("-" * 60)
    print(f"{'bucket':<12}{'n':>5}{'clean':>7}{'dirty':>7}{'total_ov':>10}{'worst':>7}")
    for name, _lo, _hi, _w in BUCKETS:
        b = [r for r in ok if r.bucket == name]
        if not b:
            continue
        bt = sum(r.witnesses for r in b)
        bc = sum(1 for r in b if r.witnesses == 0)
        worst = max((r.witnesses for r in b), default=0)
        print(f"{name:<12}{len(b):>5}{bc:>7}{len(b) - bc:>7}{bt:>10}{worst:>7}")
    print("=" * 60)

    out = Path(__file__).parent / f"gate_result_{engine_name.replace(':', '_')}.json"
    out.write_text(json.dumps([r.__dict__ for r in results], indent=0))
    print(f"per-structure results -> {out}")


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("select")
    w = sub.add_parser("worst")
    w.add_argument("--from", dest="from_result", default="benchmarks/gate_result_puzzler.json")
    w.add_argument("--n", type=int, default=40)
    r = sub.add_parser("run")
    r.add_argument("--engine", default="puzzler")
    r.add_argument("--workers", type=int, default=10)
    r.add_argument("--set", dest="set_name", choices=["hard", "worst"], default="hard")
    args = ap.parse_args()
    if args.cmd == "select":
        select_hard_set()
    elif args.cmd == "worst":
        build_worst_set(Path(args.from_result), args.n)
    else:
        set_path = WORST_SET_JSON if args.set_name == "worst" else HARD_SET_JSON
        run_gate(args.engine, args.workers, set_path)


if __name__ == "__main__":
    main()
