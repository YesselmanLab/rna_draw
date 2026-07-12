/*
 * vendor_instrument.c -- oracle-side parity instrumentation for the
 * vendored RNApuzzler/RNAturtle layout core (vendor/RNApuzzler/), compiled
 * ONLY when RNA_DRAW_BUILD_ORACLE is on (see CMakeLists.txt).
 *
 * MECHANISM (named per the plan review's refinement #2 -- see
 * .claude/plans/current-plan.md, "cpp-reviewer APPROVED" section): the
 * vendored layout files are `#include`d `.inc` amalgams, each compiled
 * into ONE translation unit per engine (RNApuzzler.c / RNAturtle.c already
 * do this). Their tree/box/config-change internals (configtree.inc,
 * handleConfigChanges.inc, ...) are PRIVATE (`static`) functions, invisible
 * from any other TU -- a plain wrapper TU (like this one) CANNOT call them
 * directly. Rather than editing the vendored .inc files to expose them
 * (forbidden -- see this repo's SAFETY rules), this TU `#include`s the SAME
 * .inc amalgam a SECOND time, into ITS OWN translation unit: every
 * `PRIVATE`/`static` function becomes locally callable here, with the exact
 * same logic RNApuzzler.c's copy runs (C header guards are per-TU, so this
 * compiles a second, independent copy -- exactly how RNApuzzler.c and
 * RNAturtle.c already each get their own copy of the shared `.inc`s today).
 * `MACRO INTERPOSITION` proper (`#define <name> <shim>` before an include,
 * to redirect calls INSIDE an *existing* TU) is reserved for a later step
 * that needs to observe calls RNApuzzler.c's own resolver makes internally
 * (`dump_change_trace`, Milestone A step 7) -- this dump does not need
 * that, since it drives the tree-build pipeline itself rather than
 * observing another TU's private call graph. Only one piece of shared
 * state crosses TUs: `rnadraw_clearance_value` (declared `extern` by
 * `definitions.inc`, defined once in RNApuzzler.c) -- this dump does not
 * touch it (uses the stock `epsilonRecognize`/`epsilonFix` via clearance
 * 0 => `_rnadraw_clearance()` falls back to 1.0).
 *
 * THIS SLICE (Milestone A step 4): `rnadraw_oracle_dump_tree` runs the
 * turtle pass + `buildConfigtree` + `updateBoundingBoxes` (mirroring
 * `RNApuzzler.c:421-476`'s setup through the resolver -- the resolver
 * itself is NOT run) and serializes the resulting T1 tree to JSON text, in
 * the same field-name shape `bindings.cpp`'s `dump_config_tree_binding`
 * (native `_layout_core.dump_tree`) returns as a `list[dict]`, so
 * `tests/test_native_parity.py` can compare both sides after
 * `json.loads()`.
 *
 * MILESTONE A STEP 5 addition: `rnadraw_oracle_dump_detections` runs the
 * SAME setup (through `updateBoundingBoxes`) and then computes the FULL
 * intersection detection set over the resulting tree -- every intersecting
 * non-root node pair (`intersectNodeNode`, called for all `n*(n-1)/2` pairs
 * in the SAME DFS pre-order `id` numbering `dump_tree`'s JSON already uses)
 * plus every direct root child that intersects the exterior baseline
 * (`intersectNodeExterior`, forced on via `checkExteriorIntersections = 1`).
 * This mirrors `include/rna_layout/intersect_tree.hpp`'s
 * `detect_intersections` exactly -- see that file's header for why this
 * "check every pair" detection set has no single vendored counterpart (the
 * resolver only ever queries specific pairs on demand): calling the
 * vendored `intersectNodeNode`/`intersectNodeExterior` directly, in a loop
 * this file adds, is NOT an edit to vendored logic (same "second
 * independent copy of the .inc amalgam, in this TU" mechanism the rest of
 * this file already uses -- see the header above).
 *
 * MILESTONE A STEP 7 addition: `rnadraw_oracle_dump_change_trace` runs the
 * T0 tree, then `checkAndFixIntersections` with `checkSiblingIntersections
 * = 1` / `checkAncestorIntersections = 0` / `optimize = 0` (the SIBLING-only
 * resolver path), and returns the ORDERED sequence of config-change
 * decisions it made. This is the ONE dump that genuinely needs
 * MACRO INTERPOSITION proper (not just "a second independent copy of the
 * .inc amalgam" like every dump above): the private, `static`
 * `checkAndApplyConfigChanges` (`handleConfigChanges.inc:55`) is called
 * from INSIDE this TU's own copy of `handleSiblingIntersections.inc`'s
 * private call graph, so there is no external hook to attach to -- instead,
 * `#define checkAndApplyConfigChanges <renamed>` before including
 * `handleConfigChanges.inc` renames the REAL implementation, then this file
 * defines its OWN `checkAndApplyConfigChanges` (under the original name,
 * the only symbol left with that name once the rename's header guard
 * prevents a second, unrenamed definition) that records a trace entry
 * around a call to the renamed real one. Every later `#include` in this
 * TU's copy of the amalgam (`handleSiblingIntersections.inc`,
 * `handleAncestorIntersections.inc`, `resolveIntersections.inc`) then calls
 * THIS wrapper for every `checkAndApplyConfigChanges(...)` call site,
 * because C's header guards make `handleConfigChanges.inc` a no-op the
 * second time any of them `#include`s it. This is the mechanism the plan's
 * cpp-reviewer refinement #2 named ("MACRO INTERPOSITION ... consistent
 * with the existing rna_draw fork practice") -- it edits NO vendored `.inc`
 * logic (the renamed function's BODY, copied verbatim by the preprocessor,
 * is untouched).
 */

#include <ViennaRNA/structures/pairtable.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "includes/boundingBoxes.inc"
#include "includes/boundingWedge.inc"
#include "includes/calcDeltas.inc"
#include "includes/configtree.inc"
#include "includes/coordinates.inc"
#include "includes/definitions.inc"
#include "includes/drawingconfig.inc"
#include "includes/intersectLevelTreeNodes.inc"
#include "includes/vector_math.inc"

/*---------------------------------------------------------------------------
 *  Growable string buffer (JSON text builder)
 *--------------------------------------------------------------------------*/

typedef struct {
  char* data;
  size_t length;
  size_t capacity;
} strbuf_t;

static void strbuf_init(strbuf_t* buf) {
  buf->capacity = 4096;
  buf->length = 0;
  buf->data = (char*)vrna_alloc(buf->capacity);
  buf->data[0] = '\0';
}

static void strbuf_reserve(strbuf_t* buf, size_t extra) {
  if (buf->length + extra + 1 <= buf->capacity) return;

  while (buf->length + extra + 1 > buf->capacity) buf->capacity *= 2;

  char* grown = (char*)realloc(buf->data, buf->capacity);
  if (grown == NULL) {
    free(buf->data);
    // `fprintf_s` is a Windows/Annex-K-only extension, unavailable on this
    // platform; this is a fixed literal (no untrusted format/length input).
    // NOLINTNEXTLINE(clang-analyzer-security.insecureAPI.DeprecatedOrUnsafeBufferHandling)
    fprintf(stderr, "vendor_instrument: strbuf realloc failed\n");
    exit(EXIT_FAILURE);
  }
  buf->data = grown;
}

static void strbuf_append(strbuf_t* buf, const char* text) {
  size_t len = strlen(text);

  strbuf_reserve(buf, len);
  memcpy(buf->data + buf->length, text, len + 1);
  buf->length += len;
}

/* Appends a formatted double with full round-trip precision. */
static void strbuf_append_double(strbuf_t* buf, double value) {
  char formatted[64];

  snprintf(formatted, sizeof(formatted), "%.17g", value);
  strbuf_append(buf, formatted);
}

static void strbuf_append_int(strbuf_t* buf, int value) {
  char formatted[32];

  snprintf(formatted, sizeof(formatted), "%d", value);
  strbuf_append(buf, formatted);
}

/*---------------------------------------------------------------------------
 *  Tree -> JSON
 *--------------------------------------------------------------------------*/

static void dump_config_json(strbuf_t* buf, const config* cfg) {
  if (cfg == NULL) {
    strbuf_append(buf, "null");
    return;
  }

  strbuf_append(buf, "{\"radius\":");
  strbuf_append_double(buf, cfg->radius);
  strbuf_append(buf, ",\"min_radius\":");
  strbuf_append_double(buf, cfg->minRadius);
  strbuf_append(buf, ",\"default_radius\":");
  strbuf_append_double(buf, cfg->defaultRadius);
  strbuf_append(buf, ",\"arcs\":[");
  for (int i = 0; i < cfg->numberOfArcs; i++) {
    if (i > 0) strbuf_append(buf, ",");

    strbuf_append(buf, "{\"segments\":");
    strbuf_append_int(buf, cfg->cfgArcs[i].numberOfArcSegments);
    strbuf_append(buf, ",\"angle\":");
    strbuf_append_double(buf, cfg->cfgArcs[i].arcAngle);
    strbuf_append(buf, "}");
  }
  strbuf_append(buf, "]}");
}

static void dump_lbox_json(strbuf_t* buf, const loopBox* lbox) {
  if (lbox == NULL) {
    strbuf_append(buf, "null");
    return;
  }

  strbuf_append(buf, "{\"cx\":");
  strbuf_append_double(buf, lbox->c[0]);
  strbuf_append(buf, ",\"cy\":");
  strbuf_append_double(buf, lbox->c[1]);
  strbuf_append(buf, ",\"r\":");
  strbuf_append_double(buf, lbox->r);
  strbuf_append(buf, "}");
}

static void dump_sbox_json(strbuf_t* buf, const stemBox* sbox) {
  if (sbox == NULL) {
    strbuf_append(buf, "null");
    return;
  }

  strbuf_append(buf, "{\"ax\":");
  strbuf_append_double(buf, sbox->a[0]);
  strbuf_append(buf, ",\"ay\":");
  strbuf_append_double(buf, sbox->a[1]);
  strbuf_append(buf, ",\"bx\":");
  strbuf_append_double(buf, sbox->b[0]);
  strbuf_append(buf, ",\"by\":");
  strbuf_append_double(buf, sbox->b[1]);
  strbuf_append(buf, ",\"cx\":");
  strbuf_append_double(buf, sbox->c[0]);
  strbuf_append(buf, ",\"cy\":");
  strbuf_append_double(buf, sbox->c[1]);
  strbuf_append(buf, ",\"ex\":");
  strbuf_append_double(buf, sbox->e[0]);
  strbuf_append(buf, ",\"ey\":");
  strbuf_append_double(buf, sbox->e[1]);
  strbuf_append(buf, ",\"bulge_count\":");
  strbuf_append_int(buf, sbox->bulgeCount);
  strbuf_append(buf, ",\"bulge_dist\":");
  strbuf_append_double(buf, sbox->bulgeDist);
  strbuf_append(buf, "}");
}

/* DFS pre-order, matching `id`'s own assignment order (`configtree.inc`'s
 * `treeHandleStem`, `++(*nodeID)` immediately before recursing) -- see
 * `include/rna_layout/debug_dump.hpp`'s matching note on the native side. */
static void dump_node_json(strbuf_t* buf, const treeNode* node, short is_first) {
  if (!is_first) strbuf_append(buf, ",");

  strbuf_append(buf, "{\"id\":");
  strbuf_append_int(buf, getNodeID(node));
  strbuf_append(buf, ",\"parent_id\":");
  strbuf_append_int(buf, getNodeID(getParent(node)));
  strbuf_append(buf, ",\"loop_start\":");
  strbuf_append_int(buf, node->loop_start);
  strbuf_append(buf, ",\"stem_start\":");
  strbuf_append_int(buf, node->stem_start);
  strbuf_append(buf, ",\"cfg\":");
  dump_config_json(buf, node->cfg);
  strbuf_append(buf, ",\"lbox\":");
  dump_lbox_json(buf, node->lBox);
  strbuf_append(buf, ",\"sbox\":");
  dump_sbox_json(buf, node->sBox);
  strbuf_append(buf, "}");

  for (int i = 0; i < node->childCount; i++) dump_node_json(buf, getChild(node, i), 0);
}

/*---------------------------------------------------------------------------
 *  Shared T1-tree build (dump_tree + dump_detections setup)
 *--------------------------------------------------------------------------*/

/* Everything `rnadraw_oracle_dump_tree`/`rnadraw_oracle_dump_detections`
 * need to free once they are done with the tree. */
typedef struct {
  short* pair_table;
  tBaseInformation* base_information;
  double* x;
  double* y;
  treeNode* tree;
} t1_tree_t;

static void t1_tree_free(t1_tree_t* built) {
  if (built->tree != NULL) freeTree(built->tree);
  free(built->x);
  free(built->y);
  free(built->base_information);
  free(built->pair_table);
}

/* Runs the turtle pass + buildConfigtree + updateBoundingBoxes on
 * `structure` (mirroring `RNApuzzler.c:421-476` through, but not
 * including, the resolver). Returns 1 on success (`*out` populated;
 * caller must `t1_tree_free(out)` when done) or 0 on a malformed/
 * degenerate structure (`*out`'s fields are NULL/zeroed either way). */
static short build_t1_tree(const char* structure, double paired, double unpaired,
                           t1_tree_t* out) {
  out->pair_table = NULL;
  out->base_information = NULL;
  out->x = NULL;
  out->y = NULL;
  out->tree = NULL;

  short* pair_table = vrna_ptable(structure);

  if (pair_table == NULL) return 0;

  int length = pair_table[0];

  if (length <= 0) {
    free(pair_table);
    return 0;
  }
  out->pair_table = pair_table;

  tBaseInformation* baseInformation = vrna_alloc((length + 1) * sizeof(tBaseInformation));

  for (int i = 0; i <= length; i++) {
    baseInformation[i].baseType = TYPE_BASE_NONE;
    baseInformation[i].distance = unpaired;
    baseInformation[i].angle = 0.0;
    baseInformation[i].config = NULL;
  }
  out->base_information = baseInformation;

  cfgGenerateConfig(pair_table, baseInformation, unpaired, paired);
  computeAffineCoordinates(pair_table, paired, unpaired, baseInformation);

  double* x = (double*)vrna_alloc(length * sizeof(double));
  double* y = (double*)vrna_alloc(length * sizeof(double));

  affineToCartesianCoordinates(baseInformation, length, x, y);
  out->x = x;
  out->y = y;

  double distBulge = sqrt(unpaired * unpaired - 0.25 * unpaired * unpaired);

  treeNode* tree = buildConfigtree(pair_table, baseInformation, x, y, distBulge);
  out->tree = tree;

  /* Fully-initialized defaults (same constructor `bindings.cpp`'s
   * `PuzzlerOptions` RAII wrapper uses), overriding only the two fields
   * `updateBoundingBoxes` reads. */
  vrna_plot_options_puzzler_t* puzzler_options = vrna_plot_options_puzzler();

  puzzler_options->paired = paired;
  puzzler_options->unpaired = unpaired;
  updateBoundingBoxes(tree, puzzler_options);
  vrna_plot_options_puzzler_free(puzzler_options);

  return 1;
}

/*---------------------------------------------------------------------------
 *  Public entry points
 *--------------------------------------------------------------------------*/

/*
 * Returns a malloc'd JSON array of `structure`'s T1 tree (see
 * `build_t1_tree`), one object per node (see this file's header for the
 * schema). Caller owns the returned buffer; free() it. Returns NULL on a
 * malformed/degenerate structure (caller should treat that as an error,
 * same contract as `vrna_plot_coords_puzzler` returning 0).
 */
char* rnadraw_oracle_dump_tree(const char* structure, double paired, double unpaired) {
  t1_tree_t built;

  if (!build_t1_tree(structure, paired, unpaired, &built)) return NULL;

  strbuf_t buf;

  strbuf_init(&buf);
  strbuf_append(&buf, "[");
  dump_node_json(&buf, built.tree, 1);
  strbuf_append(&buf, "]");

  t1_tree_free(&built);

  return buf.data;
}

/*---------------------------------------------------------------------------
 *  Detection set (Milestone A step 5)
 *--------------------------------------------------------------------------*/

/* Growable flat array of `treeNode*`, filled by `flatten_tree` in the SAME
 * DFS pre-order `getNodeID` assignment uses (`configtree.inc`'s
 * `treeHandleStem`), so `nodes[i]` and `getNodeID(nodes[i])` agree -- the
 * same invariant `dump_node_json` already relies on. */
typedef struct {
  treeNode** items;
  int count;
  int capacity;
} node_list_t;

static void node_list_init(node_list_t* list) {
  list->capacity = 64;
  list->count = 0;
  list->items = (treeNode**)vrna_alloc((size_t)list->capacity * sizeof(treeNode*));
}

static void node_list_push(node_list_t* list, treeNode* node) {
  if (list->count >= list->capacity) {
    list->capacity *= 2;
    list->items = (treeNode**)realloc(list->items, (size_t)list->capacity * sizeof(treeNode*));
  }
  list->items[list->count++] = node;
}

static void flatten_tree(treeNode* node, node_list_t* list) {
  node_list_push(list, node);
  for (int i = 0; i < node->childCount; i++) flatten_tree(getChild(node, i), list);
}

/*
 * Returns a malloc'd JSON array of `[[node1_id, node2_id, "type"], ...]`
 * over the FULL detection set of `structure`'s T1 tree: every intersecting
 * non-root node pair (`intersectNodeNode`, all `n*(n-1)/2` pairs, `id_i <
 * id_j`, in that nested-loop discovery order) then every direct root child
 * that intersects the exterior baseline (`intersectNodeExterior`, forced on
 * via `checkExteriorIntersections = 1` -- `node2_id` is the root's own id,
 * `0`, for these) -- see `include/rna_layout/intersect_tree.hpp`'s
 * `detect_intersections`, which this exactly mirrors. `"type"` is
 * `intersectionTypeToString`'s own short code (`"LxL"`, `"SxS"`, ...,
 * `"EXT"`), so both sides serialize identically without a translation
 * table on the Python test side. Caller owns the returned buffer; free()
 * it. Returns NULL on a malformed/degenerate structure.
 */
char* rnadraw_oracle_dump_detections(const char* structure, double paired, double unpaired) {
  t1_tree_t built;

  if (!build_t1_tree(structure, paired, unpaired, &built)) return NULL;

  node_list_t nodes;

  node_list_init(&nodes);
  flatten_tree(built.tree, &nodes);

  vrna_plot_options_puzzler_t* puzzler_options = vrna_plot_options_puzzler();

  puzzler_options->checkExteriorIntersections = 1;

  strbuf_t buf;

  strbuf_init(&buf);
  strbuf_append(&buf, "[");
  short is_first = 1;

  for (int i = 1; i < nodes.count; i++) {
    for (int j = i + 1; j < nodes.count; j++) {
      intersectionType it = intersectNodeNode(nodes.items[i], nodes.items[j]);
      if (it == noIntersection) continue;

      if (!is_first) strbuf_append(&buf, ",");
      is_first = 0;
      strbuf_append(&buf, "[");
      strbuf_append_int(&buf, getNodeID(nodes.items[i]));
      strbuf_append(&buf, ",");
      strbuf_append_int(&buf, getNodeID(nodes.items[j]));
      strbuf_append(&buf, ",\"");
      strbuf_append(&buf, intersectionTypeToString(it));
      strbuf_append(&buf, "\"]");
    }
  }

  for (int i = 1; i < nodes.count; i++) {
    if (getParent(nodes.items[i]) != built.tree) continue;
    if (!intersectNodeExterior(nodes.items[i], puzzler_options)) continue;

    if (!is_first) strbuf_append(&buf, ",");
    is_first = 0;
    strbuf_append(&buf, "[");
    strbuf_append_int(&buf, getNodeID(nodes.items[i]));
    strbuf_append(&buf, ",");
    strbuf_append_int(&buf, getNodeID(built.tree));
    strbuf_append(&buf, ",\"EXT\"]");
  }

  strbuf_append(&buf, "]");

  vrna_plot_options_puzzler_free(puzzler_options);
  free(nodes.items);
  t1_tree_free(&built);

  return buf.data;
}

/*---------------------------------------------------------------------------
 *  Change trace (Milestone A step 7) -- MACRO INTERPOSITION, see this
 *  file's header for the full mechanism explanation.
 *--------------------------------------------------------------------------*/

/* Growable list of trace entries, filled by the `checkAndApplyConfigChanges`
 * wrapper below and read back out by `rnadraw_oracle_dump_change_trace`. */
typedef struct {
  int node_id;
  intersectionType type;
  double* deltas;
  int num_deltas;
  short accepted;
} trace_entry_t;

typedef struct {
  trace_entry_t* items;
  int count;
  int capacity;
} trace_list_t;

/* Non-NULL only for the duration of one `rnadraw_oracle_dump_change_trace`
 * call (single-threaded, same discipline `rnadraw_clearance_value`,
 * `definitions.inc:39`, already uses for its own call-scoped global). */
static trace_list_t* g_trace = NULL;

static void trace_list_init(trace_list_t* list) {
  list->capacity = 16;
  list->count = 0;
  list->items = (trace_entry_t*)vrna_alloc((size_t)list->capacity * sizeof(trace_entry_t));
}

static void trace_list_free(trace_list_t* list) {
  for (int i = 0; i < list->count; i++) free(list->items[i].deltas);
  free(list->items);
}

static void trace_list_push(trace_list_t* list, trace_entry_t entry) {
  if (list->count >= list->capacity) {
    list->capacity *= 2;
    list->items = (trace_entry_t*)realloc(list->items, (size_t)list->capacity * sizeof(trace_entry_t));
  }
  list->items[list->count++] = entry;
}

/* Renames the REAL `checkAndApplyConfigChanges` definition
 * (`handleConfigChanges.inc:55`) for the duration of this one `#include`;
 * the vendored BODY is copied by the preprocessor unmodified. */
#define checkAndApplyConfigChanges rnadraw_traced_checkAndApplyConfigChanges_real
#include "includes/handleConfigChanges.inc"
#undef checkAndApplyConfigChanges

/* The trace-recording wrapper: the ONLY definition left under the original
 * name in this TU (`handleConfigChanges.inc`'s header guard makes every
 * later `#include` of it, from `handleSiblingIntersections.inc`/
 * `handleAncestorIntersections.inc`, a no-op) -- so every
 * `checkAndApplyConfigChanges(...)` call site the resolver's private call
 * graph makes, from here on in this TU, resolves to this wrapper. */
PRIVATE short checkAndApplyConfigChanges(treeNode* tree, double* deltaCfg,
                                         const intersectionType it,
                                         vrna_plot_options_puzzler_t* puzzler) {
  short changed = rnadraw_traced_checkAndApplyConfigChanges_real(tree, deltaCfg, it, puzzler);

  if (g_trace != NULL) {
    int num_deltas = tree->cfg->numberOfArcs;
    trace_entry_t entry;

    entry.node_id = getNodeID(tree);
    entry.type = it;
    entry.num_deltas = num_deltas;
    /* `deltaCfg` was mutated IN PLACE by the "fix too small changes" step
     * inside the real function above -- read AFTER the call, matching the
     * native side's capture point (`config_changes.cpp`'s
     * `check_and_apply_config_changes`, which records post-adjustment). */
    entry.deltas = (double*)vrna_alloc((size_t)num_deltas * sizeof(double));
    for (int i = 0; i < num_deltas; i++) entry.deltas[i] = deltaCfg[i];
    entry.accepted = changed;

    trace_list_push(g_trace, entry);
  }

  return changed;
}

#include "includes/handleSiblingIntersections.inc"
#include "includes/handleAncestorIntersections.inc"
#include "includes/resolveIntersections.inc"

/*
 * Returns a malloc'd JSON array of `[{"node_id":.., "type":"BRA",
 * "deltas":[...], "accepted":true}, ...]`: the ORDERED sequence of
 * config-change decisions `checkAndFixIntersections` made on `structure`'s
 * T1 tree with `checkSiblingIntersections = 1`, `checkAncestorIntersections
 * = 0`, `optimize = 0` (the SIBLING-only resolver path, Milestone A step
 * 7) -- see `include/rna_layout/resolve.hpp`'s `ChangeTraceEntry`, which
 * this exactly mirrors. Caller owns the returned buffer; free() it.
 * Returns NULL on a malformed/degenerate structure.
 *
 * WARNING (see `benchmarks/hard_set_sibling_only_oracle_hangs.json`): the
 * vendored `checkAndFixIntersections` does not terminate on every
 * structure under this option combination -- callers MUST NOT invoke this
 * on a structure from that exclusion list (or any other not already known
 * to terminate) without an external timeout; this file adds none, per its
 * "expose dumps, don't change vendored logic" scope.
 */
char* rnadraw_oracle_dump_change_trace(const char* structure, double paired, double unpaired,
                                       int max_config_changes) {
  t1_tree_t built;

  if (!build_t1_tree(structure, paired, unpaired, &built)) return NULL;

  vrna_plot_options_puzzler_t* puzzler_options = vrna_plot_options_puzzler();

  puzzler_options->paired = paired;
  puzzler_options->unpaired = unpaired;
  puzzler_options->checkSiblingIntersections = 1;
  puzzler_options->checkAncestorIntersections = 0;
  puzzler_options->optimize = 0;
  puzzler_options->numberOfChangesAppliedToConfig = 0;
  puzzler_options->maximumNumberOfConfigChangesAllowed =
      (max_config_changes <= 0) ? 25000 : max_config_changes;

  trace_list_t trace;

  trace_list_init(&trace);
  g_trace = &trace;
  checkAndFixIntersections(built.tree, 0, puzzler_options);
  g_trace = NULL;

  strbuf_t buf;

  strbuf_init(&buf);
  strbuf_append(&buf, "[");
  for (int i = 0; i < trace.count; i++) {
    if (i > 0) strbuf_append(&buf, ",");

    trace_entry_t* entry = &trace.items[i];
    strbuf_append(&buf, "{\"node_id\":");
    strbuf_append_int(&buf, entry->node_id);
    strbuf_append(&buf, ",\"type\":\"");
    strbuf_append(&buf, intersectionTypeToString(entry->type));
    strbuf_append(&buf, "\",\"deltas\":[");
    for (int j = 0; j < entry->num_deltas; j++) {
      if (j > 0) strbuf_append(&buf, ",");
      strbuf_append_double(&buf, entry->deltas[j]);
    }
    strbuf_append(&buf, "],\"accepted\":");
    strbuf_append(&buf, entry->accepted ? "true" : "false");
    strbuf_append(&buf, "}");
  }
  strbuf_append(&buf, "]");

  trace_list_free(&trace);
  vrna_plot_options_puzzler_free(puzzler_options);
  t1_tree_free(&built);

  return buf.data;
}

const char* rnadraw_oracle_instrumentation_version(void) {
  return "vendor_instrument v3: turtle dump + dump_tree (config tree + "
         "bounding boxes, Milestone A step 4) + dump_detections "
         "(intersection detection set, Milestone A step 5) + "
         "dump_change_trace (SIBLING resolver config-change trace, "
         "Milestone A step 7, via macro-interposition on "
         "checkAndApplyConfigChanges -- see this file's header).";
}
