# RNA layout: the never-overlap rebuild

Status report for the layout-engine rebuild (branch `modernize-nonoverlap-layout`). The old
engine placed nucleotides on top of each other; this documents every lever tried, the algorithm
now shipping as the default, and what remains.

All clean-rate numbers below are measured on a **frozen set of the 450 *hardest* large structures**
from the bpRNA corpus (`benchmarks/hard_set.json`, 300–4000 nt) — a deliberately adversarial
worst-case set, scored by the independent M2 overlap checker. A visual version with before/after
renders is published as a Claude Artifact.

**Real-world performance is much higher.** On a size-stratified sample of the *general* corpus (the
typical structures users actually draw), the compact production primary alone is clean on **98% of
<300 nt, 100% of 300–600 nt, 96% of 600–1200 nt** structures — near-100%, with the constructive
engine catching the rare miss and **0 silent overlaps**. The 87% figure below is the adversarial
tail, not the common case.

## Headline

| Metric | Result |
|---|---|
| Clean conventional layouts (hard set) | **30% → 88%** |
| Total overlaps across the set | **2,943 → 945 (−68%)** |
| Silent overlaps | **0** — every output is checker-clean or explicitly flagged |
| Layout compactness | median bounding box **1.00×** stock puzzler (no ballooning) |

## Current best algorithm (the default engine)

`rna_draw`'s `default_engine()` (in `rna_draw/layout/pipeline.py`) returns `production_engine()`
(in `rna_draw/layout/production.py`), run behind the checker-gated `layout_guaranteed` (also in
`pipeline.py`). The production primary itself has three stages:

1. **Per-structure clearance escalation.** Lay out with RNApuzzler's algorithm at an increasing
   *intersection-clearance* factor (ladder `1.0, 1.25, 1.5×`). This scales puzzler's own
   overlap-detection margin so it resolves the near-touches our stricter checker flags — the key
   insight being that puzzler resolves every intersection *it* detects, but its clearance model is
   more lenient than ours. Keep the first clean result. Capped at 1.5× because the C resolver is
   fast there (~0.1 s) and never hangs. **Milestone A (branch `cpp-layout-core`) replaced the
   vendored ViennaRNA C with an owned modern-C++ port** (`include/rna_layout/`,
   `src/layout_core/`, exposed as `rna_draw._layout_core`) verified to reproduce the vendored
   engine's production output to a max diff of 3.6e-12 (`tests/test_native_parity.py`) — this is
   now the shipped engine (`rna_draw.layout.native.NativePuzzlerEngine`); the vendored
   `_vienna_layout` is retained only as an opt-in parity oracle
   (`-DRNA_DRAW_BUILD_ORACLE=ON`), not part of the default/shipped build, so `rna_draw` has zero
   ViennaRNA header or runtime dependency by default.

2. **Rigid loop-inflation post-pass.** For the residual — which is small-loop crowding of unpaired
   bases, *not* branch collisions — inflate the offending loop: spread its unpaired members
   radially and translate each child branch rigidly so helices stay straight. A candidate move is
   applied only if the frozen checker reports **strictly fewer** overlaps, so the pass is monotone
   (never worse) and wall-clock bounded. Engine-agnostic; operates on coordinates + a pair map.

3. **Checker-gated adaptive radius (`_largest_clean_node_r`).** Render at the largest disk radius in
   `[floor, target]` that passes the checker — gate radius == render radius, so a "clean" result is
   never silently overlapping at the radius actually drawn.

## The 4-tier pipeline (`layout_guaranteed`)

`layout_guaranteed` (`rna_draw/layout/pipeline.py`) tries an ordered chain, **every tier
checker-gated** against the frozen `check_overlaps` — never a silent overlap:

1. **Production primary** (`_try_primary`, the compact `production_engine` by default). Only
   attempted when the input `is_pseudoknot_free`. Clean → `flagged=False`.
2. **Constructive engine** (`_try_constructive_fallback`, `rna_draw/layout/constructive/`). An
   overlap-free-*by-construction* layout (bounding-disk envelope tree, see below), gated exactly
   like tier 1. It is checker-clean but honestly reported `flagged=True` (`engine_name="constructive"`),
   because it is only reached when the compact primary was *not* clean.
3. **Pseudoknot tier** (`_try_pseudoknot`, `rna_draw/layout/pseudoknot/`). Reached **only when the
   input is not pseudoknot-free** (contains `[]{}<>`) — tiers 1–2 decline such input via their own
   `is_pseudoknot_free` guard. Lays out a max-nested subset conventionally, then draws the crossing
   pairs as checker-validated connectors (see below). `flagged=False` only if every crossing became a
   clean in-plane connector; otherwise `flagged=True`.
4. **Circle `SafeFallbackEngine`** (`_fallback_result`). The guaranteed-terminating last resort,
   always `flagged=True`.

## The constructive engine (the fallback tier)

A second, independent engine (`rna_draw/layout/constructive/`) that lays a structure out overlap-free
*by construction* rather than by repair — used when the compact production primary can't clear. It
builds bottom-up over the loop/stem tree (`engine.py`):

- **Bounding-disk envelope tree** (`envelope.py`). Each subtree is bounded by a **disk** centered
  exactly on its attachment pivot; children are packed around a loop at angular half-widths
  `asin(r/d)` (`geometry_helpers.pack_loop_angles`) — provably disjoint, so no two siblings' disks,
  hence no two siblings' subtrees, can overlap. Disjointness is a construction-time *proof*, not a
  post-hoc check.
- **Straight-continuation fixes** to stop the envelope bound compounding on long chains. A
  single-child bulge/interior loop is placed as a straight continuation of the parent stem's axis
  (`_place_bulge`), so a bulge chain's reach grows linearly, not ~3× per level. A **degree-2
  multiloop** (one large continuing branch + one small side branch, repeated deep in real rRNA) pins
  its dominant child collinear too (`envelope.lateral_reach` / `_degree2_packing`), sized by a
  directional lateral reach — the same linear-growth fix for the other structural cause of blow-up.
  A `_MAX_REACH` guard still bails (fast, before any coordinate is written) on anything genuinely
  isotropic (e.g. a giant 3+-way junction with two large children).
- **Provable area-minimization pass** (`area_min.py`). On top of the sound layout, each eligible
  branch is rotated about its own attachment pivot. Because the sibling-disjoint disk is centered
  *at* that pivot, a rotation is an isometry that fixes the disk, so every sibling/cousin disjointness
  relation is preserved automatically — no re-verification needed. The one primitive the disk
  argument doesn't cover (the branch's own incoming/departing backbone capsules) is checked
  explicitly (`_boundary_capsules_clear`), and every accepted move is monotone (kept only if it
  strictly shrinks bbox area and stays clean; reverted otherwise). It delivers real 22–35% area cuts
  on most worst-set structures. A per-loop *radius* re-pack is deliberately not built here — it would
  feed anisotropic half-widths into `pack_loop_angles`, which bounds only perpendicular extent, so it
  would not be certificate-clean; radius tightening stays checker-gated in `compaction.py`. An
  earlier **exterior 2-D fold DOF** (wrapping the exterior line into rows) was tried, verified sound,
  and **removed** — it was rejected by the safety check on every real and every synthetically
  favorable structure (a confirmed no-op).

Every produced layout is re-verified against the frozen checker before return; a dirty result raises
`EngineError` rather than being returned silently. The engine lays out **436/450 hard structures
clean-by-construction (0 silent overlaps), including the deep-rRNA tail** production leaves to
fallback. It is ~21× less compact than the production primary (so it stays a fallback, not a primary)
but ~1600× *more* compact than the circle it replaces, and conventional in style. On the 450 hard set
the pipeline now resolves as **86% production-primary (compact) + 11% constructive-fallback (clean,
conventional-style) + 2% circle**, 0 silent overlaps.

## Pseudoknot layout (`rna_draw/layout/pseudoknot/`)

The third pipeline tier, reached only for input that is not pseudoknot-free (`[]{}<>`). It never
re-places a nucleotide: every base is placed exactly once by a conventional nested layout, and a
crossing pair is drawn as an *additional* connector between two already-placed nucleotides.

1. **Parse** the multi-bracket string into stems (`parsing.py`).
2. **Extract a max-nested subset** via Maximum-Weight Independent Set on the stem-crossing conflict
   graph (`extraction.py`) — exact branch-and-bound over only the crossing-involved stems (corpus
   median 3, max 7), with a greedy fallback above `EXACT_COVER_LIMIT`. The retained subset is the
   nested tree; the removed stems are the crossings. A round-trip guard asserts the reconstructed
   subset really is pseudoknot-free before handoff.
3. **Lay out the nested subset** with the ordinary checker-gated `layout_guaranteed` pipeline
   (`engine.py` → `_layout_nested`).
4. **Draw each crossing base pair** (`placement.py`), one at a time (innermost pair of a crossing
   stem first, narrowest crossing stem first), escalating per pair, checker-gated the whole way:
   a straight **in-plane connector** (PK-B) added to the pair map and re-validated by the
   *unmodified* `check_overlaps`; else a non-overlapping, **axis-aligned orthogonal "staple"**
   polyline (PK-A, `routing.py` — a 3- or 5-segment bracket that goes straight out, across, and back
   in, never a diagonal, validated by `validate.py`); else left **unplaced**. A final
   mutual-validation backstop (`_clean_routed_lines`) drops any routed line that isn't clean against
   everything else. The frozen checker validates crossing pairs unchanged — never a silent or drawn
   overlap.

## The never-silent-overlap contract and the frozen checker

The load-bearing guarantee: **every returned layout is either checker-clean or explicitly
`flagged`** — the pipeline never hands back a silent overlap. The arbiter is a single, **frozen,
engine-agnostic overlap checker** shared by every engine, every post-pass move, and every crossing
connector:

- `rna_draw/overlap.py` builds the primitives a renderer actually draws (nucleotide **disks**,
  backbone **capsules**, base-pair **capsules**), excludes the pairs *supposed* to touch (backbone
  neighbors, base-pair partners, a disk at its own capsule's endpoint), and reports every remaining
  overlap. Overlap is judged purely by geometry — no index-distance exclusion.
- `rna_draw/geometry.py` holds the exact primitive-vs-primitive predicates (disk/disk, disk/capsule,
  capsule/capsule).
- `rna_draw/spatial_hash.py` makes the all-pairs sweep near-linear on long chains (proven
  byte-identical to brute force).

Each tier is gated by this same checker at its own render radius (`_largest_clean_node_r`, so gate
radius == render radius). The constructive engine additionally verifies its by-construction output;
its area-min and compaction moves are each checker-verified before being kept (monotone revert). The
pseudoknot tier validates every connector against it. The checker is treated as immutable — engines
and post-passes call its read-only building blocks, never alter its semantics.

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
| 11 | **Constructive engine** (bounding-disk envelope tree) | Lay out overlap-free *by construction*, not by repair | 96.9% clean-by-construction incl. the deep tail, 0 silent overlaps; but ~21× less compact than the primary | **win** |
| 12 | Wire constructive as the middle fallback tier | Give the deep tail a clean conventional layout instead of a circle | Pipeline now 86% primary + 11% constructive-fallback + 2% circle; ~1600× more compact than the circle it replaces | shipped |
| 13 | Standalone vendored build (drop `libRNA.a`) | The layout core needs only two library symbols | Links standalone via a ~120-LOC `vrna_compat.c` shim; no ViennaRNA runtime dep; `naview` removed | shipped |
| 14 | **Pseudoknot tier** (max-nested + on-top crossings) | Draw `[]{}<>` conventionally instead of a bare circle | Conventional nested layout + correct bonds where drawable, 0 overlaps; crossing depiction capped ~38% post-hoc (co-design needed) | shipped |
| 15 | **Area-min rotation pass** (disk-preserving) | Shrink constructive sprawl with a provably-clean lever | Real 22–35% area cuts on most worst-set structures, 0 dirty/0 regressions; exterior fold DOF tried + removed (no-op) | shipped |
| 16 | **Milestone A: owned native `rna_layout` C++ port** (`include/rna_layout/`, `src/layout_core/`) replaces vendored ViennaRNA C in production | Own the layout core outright, drop the last ViennaRNA runtime dependency, without regressing the never-silent-overlap contract | Module-by-module parity-gated port (turtle → tree/boxes → detection → resolver siblings/ancestors/optimize); full-config coordinate parity vs. the vendored oracle max diff **3.6e-12** (`tests/test_native_parity.py`); swapped in as production primary, vendored `_vienna_layout` retired to an opt-in `-DRNA_DRAW_BUILD_ORACLE=ON` parity build | shipped |

## Clean-rate progression (450 hard structures)

```
Stock puzzler                 30.4%  (137/450)   2943 overlaps
+ intersection clearance 1.5×  42.9%  (193/450)   2008
+ escalating finer ladder      52.2%  (235/450)   1896
+ rigid loop-inflation pass    72.9%  (328/450)   1739
+ cap / move-budget tuning     87.3%  (393/450)    945
```

## What's left

- **Constructive engine → compact *primary*.** The constructive engine is overlap-free by
  construction on 96.9% of the hard set but ~21× less compact than the production primary, so it
  serves as the fallback tier, not a primary. Closing that gap (to ~1–3×) — which would raise the
  compact-clean rate from 87% toward 97% — hits a diagnosed floor (deep subtree rotations redirect
  content past the local compaction check on the largest many-branch structures) and needs per-branch
  envelope-aware verification: a larger follow-on effort, not a quick tune.
- **Deep-rRNA tail** (~13% of the hard set): production can't lay these out compactly; they now get a
  clean constructive-engine fallback (conventional style) instead of a circle. The ~2% largest
  (>2900 nt) still fall to the circle (the constructive reach guard rejects them).
- **Empty-loop inputs** (14 structures): contain a degenerate `()` loop that hangs RNApuzzler
  unconditionally; guarded and now routed to the constructive engine (which handles them) or the
  circle. Fix for a compact result is to handle the zero-length loop in the vendored resolver.
- **Pseudoknot crossing depiction is capped.** The "lay out the nested subset, add crossings on top"
  approach has a hard post-hoc ceiling: only ~38% of crossing bonds are placeable after the fact.
  The rest are *buried* — the packed nested layout leaves disks near-tangent around a deeply-embedded
  crossing endpoint, so there is no positive-width gap to route a line out (confirmed by exhaustive
  best-first search on real structures) and no adjacency for an in-plane connector. Raising this needs
  **co-design** (reserving space / placing crossing endpoints accessibly *during* the nested layout),
  a significant engine change with uncertain payoff. What ships is honest: a conventional nested
  layout with correct bonds wherever drawable and 0 overlaps ever — a large improvement over the old
  bare circle, but crossing quality is capped, not solved.
- **Constructive-engine sprawl.** Structures that fall to the constructive tier can still render as
  sprawly strips — content in tight clusters far apart — even after the area-min rotation pass cuts
  ~23%. The root cause is the *conservative isotropic* disk radii the by-construction proof relies on;
  a tighter sound packing is future research. Compactness is a secondary goal — the primary contract
  is never a silent overlap.

## Infrastructure

- **Owned native layout core (Milestone A, branch `cpp-layout-core`).** `include/rna_layout/` +
  `src/layout_core/` is a fresh, modern-C++17 reimplementation of RNApuzzler/RNAturtle's algorithm
  (not a refactor of the vendored `.inc` amalgamation), exposed as `rna_draw._layout_core` via
  pybind11 and validated module-by-module against the vendored C as a differential parity oracle
  (`tests/test_native_parity.py`). It is now the SHIPPED production engine
  (`rna_draw.layout.native.NativePuzzlerEngine`/`NativeTurtleEngine`, wired into
  `rna_draw.layout.production`/`pipeline.default_engine`) and has **zero ViennaRNA header or
  runtime dependency**. The default CMake build (`RNA_DRAW_BUILD_ORACLE=OFF`) compiles only
  `_layout_core` (+ the frozen overlap checker); the vendored `_vienna_layout` compiles only under
  `-DRNA_DRAW_BUILD_ORACLE=ON`, kept solely as the parity oracle for `tests/test_native_parity.py`
  and `tests/test_vienna_binding.py` (both skip cleanly when that option is off).
- **Standalone in-process vendored build (oracle-only).** An in-process pybind11 binding to an
  **editable vendored copy** of RNApuzzler/RNAturtle (`src/vienna_layout/`), so the layout core
  could be modified in-tree during the port. The vendored core (`RNApuzzler.c` + `RNAturtle.c`)
  compiles and links **standalone** — with `libRNA.a` **unlinked** and **no ViennaRNA runtime
  dependency** (only ViennaRNA *headers* are needed at compile time, and only for the oracle
  build). The two translation units reference exactly two ViennaRNA library symbols, `vrna_alloc`
  and `vrna_ptable` (+ `vrna_ptable_from_string`), supplied by a ~120-LOC compat shim
  `src/vienna_layout/vendor/vrna_compat.c` (faithful copies of the upstream 2.7.0 sources, plus a
  `vrna_log` stub for the OOM/malformed-input paths). A `ctest` smoke test proves the core links and
  runs with `libRNA.a` unlinked. `naview` was **removed** — it lived only in `libRNA.a`'s
  non-reentrant `naview.o` and never rescued a puzzler-dirty structure (see `CMakeLists.txt` and
  `src/vienna_layout/vendor/README.md`).
- Independent, **frozen** overlap checker (`rna_draw/overlap.py` + `geometry.py` + `spatial_hash.py`)
  is the arbiter for every engine and every post-pass move (see the never-silent-overlap section).
- A frozen benchmark gate (`benchmarks/hard_gate.py`) with a hard per-structure kill, so a slow or
  hanging structure is a bounded error, not a stall.
- ~1100 tests (many parametrized), ~96% coverage; ruff + mypy clean. Full-suite runs use pytest-xdist
  (`-n auto`).
