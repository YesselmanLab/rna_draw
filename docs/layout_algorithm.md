# RNA layout: the never-overlap rebuild

Status report for the layout-engine rebuild (branch `modernize-nonoverlap-layout`). The old
engine placed nucleotides on top of each other; this documents every lever tried, the algorithm
now shipping as the default, and what remains.

All clean-rate numbers are measured on a **frozen set of the 450 hardest large structures** from
the bpRNA corpus (`benchmarks/hard_set.json`, 300–4000 nt), scored by the independent M2 overlap
checker. A visual version with before/after renders is published as a Claude Artifact.

## Headline

| Metric | Result |
|---|---|
| Clean conventional layouts (hard set) | **30% → 88%** |
| Total overlaps across the set | **2,943 → 945 (−68%)** |
| Silent overlaps | **0** — every output is checker-clean or explicitly flagged |
| Layout compactness | median bounding box **1.00×** stock puzzler (no ballooning) |

## Current best algorithm (the default engine)

`rna_draw`'s `default_engine()` is `production_engine()` (in `rna_draw/layout/production.py`),
run behind the checker-gated `layout_guaranteed`. Three stages:

1. **Per-structure clearance escalation.** Lay out with ViennaRNA's RNApuzzler (vendored,
   in-process) at an increasing *intersection-clearance* factor (ladder `1.0, 1.25, 1.5×`). This
   scales puzzler's own overlap-detection margin so it resolves the near-touches our stricter
   checker flags — the key insight being that puzzler resolves every intersection *it* detects,
   but its clearance model is more lenient than ours. Keep the first clean result. Capped at 1.5×
   because the C resolver is fast there (~0.1 s) and never hangs.

2. **Rigid loop-inflation post-pass.** For the residual — which is small-loop crowding of unpaired
   bases, *not* branch collisions — inflate the offending loop: spread its unpaired members
   radially and translate each child branch rigidly so helices stay straight. A candidate move is
   applied only if the frozen checker reports **strictly fewer** overlaps, so the pass is monotone
   (never worse) and wall-clock bounded. Engine-agnostic; operates on coordinates + a pair map.

3. **Checker-gated adaptive radius + guaranteed fallback.** Render at the largest disk radius in
   `[0.8×, 1.0×]target` that passes the checker. If no conventional layout clears, fall back to a
   guaranteed-clean circle, `flagged=True`. The honest contract: every result is checker-clean or
   explicitly flagged — never a silent overlap.

## Every attempt

The wins only made sense after the dead ends — they proved the overlaps were not a resolver-power
problem but first a *clearance-model mismatch*, then a *small-loop-crowding* problem.

| # | Lever | Hypothesis | Result | Verdict |
|---|---|---|---|---|
| 1 | Portfolio: best of puzzler / naview / turtle | They fail on different structures | Clean count identical (137); naview/turtle never rescue a puzzler-dirty structure | dead end |
| 2 | Enable puzzler's `allowFlipping` | Flipping resolves what rotation can't | Worse (329 vs 313 non-clean) | dead end |
| 3 | Raise config-change budget 25k → 1M | Resolver runs out of iterations | Byte-identical output — it converges, not budget-limited | dead end |
| 4 | Global loop-radius inflation | Bigger loops = more room | NaN coords past 1.5× + catastrophic slowdown | dead end |
| 5 | **Scale puzzler's intersection clearance** | Puzzler under-detects vs our checker | 2943 → 2008 overlaps; clean 30% → 43% | **win** |
| 6 | Per-structure escalating clearance ladder | Different structures need different clearance | Clean 43% → 52% | win |
| 7 | **Rigid loop-inflation post-pass (monotone)** | Residual = small-loop crowding, not collisions | Clean 52% → 73% — the biggest jump | **big win** |
| 8 | Raise post-pass cap 6 → 20, moves 12 → 40 | Most 7–20-overlap structures are the same crowding | Clean 73% → 87% | win |
| 9 | Productionize as the default engine | Make the stack rna_draw's real output | Default renders stock-dirty structures clean; honest contract | shipped |
| 10 | Spatial-hash O(L²) → O(L) for long chords | Fallback renders took minutes | Checker ~14× faster on the fallback circle; output byte-identical (proven vs brute force) | shipped |

## Clean-rate progression (450 hard structures)

```
Stock puzzler                 30.4%  (137/450)   2943 overlaps
+ intersection clearance 1.5×  42.9%  (193/450)   2008
+ escalating finer ladder      52.2%  (235/450)   1896
+ rigid loop-inflation pass    72.9%  (328/450)   1739
+ cap / move-budget tuning     87.3%  (393/450)    945
```

## What's left

- **Deep-rRNA tail** (~37 giant structures): mix local crowding with distant branch collisions that
  local rigid moves can't coordinate away. They get the guaranteed-clean circle fallback today. The
  planned fix is a bottom-up **constructive envelope engine** (overlap-free by construction), which
  would also replace the slow circle fallback with a compact one.
- **Empty-loop inputs** (14 structures): contain a degenerate `()` loop that hangs RNApuzzler
  unconditionally; guarded and routed to the safe fallback. Fix is to handle the zero-length loop in
  the vendored resolver.

## Infrastructure

- In-process pybind11 binding to an **editable vendored copy** of RNApuzzler
  (`src/vienna_layout/`), so the algorithm can be modified in-tree.
- Independent, frozen overlap checker (`rna_draw/overlap.py`) is the arbiter for every engine and
  every post-pass move.
- A frozen benchmark gate (`benchmarks/hard_gate.py`) with a hard per-structure kill, so a slow or
  hanging structure is a bounded error, not a stall.
- 596 tests, 95% coverage; ruff + mypy clean.
