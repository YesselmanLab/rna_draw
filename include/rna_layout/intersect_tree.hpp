#pragma once

/**
 * @file intersect_tree.hpp
 * @brief Intersection DETECTION over `StemBox`/`LoopBox`/`TreeNode`, ported
 *        from `intersectLevelBoundingBoxes.inc` + `intersectLevelTreeNodes.inc`
 *        + `intersectionType.inc` (Milestone A step 5).
 *
 * This is DETECTION only: every function here is a pure predicate over an
 * already-built (T0) config tree -- it reads `StemBox`/`LoopBox`/`Aabb`
 * fields but never mutates a `TreeNode`. Tree MUTATION (the resolver
 * proper -- `checkAndApplyConfigChanges`, `handleAncestorIntersections.inc`,
 * `handleSiblingIntersections.inc`) is Milestone A steps 7-9, not built
 * here.
 *
 * `detect_intersections` is the one function with NO vendored counterpart:
 * it is the exhaustive, order-independent "detection set" the plan's
 * `dump_detections` parity seam needs (`.claude/plans/current-plan.md`'s
 * "Parity oracle harness"). The vendored resolver never computes this set
 * directly -- it only ever queries SPECIFIC pairs on demand (a node against
 * its ancestors, a node's children against each other, one subtree against
 * another) via the primitives this file also exposes
 * (`intersect_node_node`, `intersect_trees`, `intersect_node_lists`), which
 * Milestone A steps 7-9 will call unchanged. `detect_intersections` composes
 * those same primitives over ALL node pairs instead, giving a
 * traversal-order-independent ground truth to diff against a matching new
 * oracle-side dump (`vendor_instrument.c`'s `rnadraw_oracle_dump_detections`,
 * which enumerates the same pairs by calling the vendored `intersectNodeNode`/
 * `intersectNodeExterior` directly -- itself not an edit to vendored logic).
 */

#include <cstdint>
#include <vector>

#include "rna_layout/tree.hpp"
#include "rna_layout/types.hpp"

namespace rna_layout {

/**
 * Mirrors the vendored `intersectionType` enum (`intersectionType.inc`)
 * value-for-value: downstream parity comparisons (native vs. the oracle's
 * JSON dump) match on these exact integers. `XxY` reads as "node1 is an X,
 * node2 is a Y" (e.g. `stem_loop` is node1's STEM intersecting node2's
 * LOOP) -- see `intersect_node_node`'s doc comment for how the vendored
 * `intersectNodeNode` establishes that convention.
 */
enum class IntersectionType : std::uint8_t {
  none = 0,
  loop_loop = 1,    // LxL
  loop_stem = 2,    // LxS
  stem_loop = 3,    // SxL
  stem_stem = 4,    // SxS
  loop_bulge = 5,   // LxB
  bulge_loop = 6,   // BxL
  stem_bulge = 7,   // SxB
  bulge_stem = 8,   // BxS
  bulge_bulge = 9,  // BxB
  siblings = 10,    // "BRA" (resolver change-context tag; not returned by
                    // any detection predicate here -- Milestone A step 7+)
  exterior = 11,    // "EXT"
};

/// `intersectionTypeToString` (`intersectionType.inc:27`), reproduced so
/// both the native and oracle-side parity dumps serialize the identical
/// short code for each type.
[[nodiscard]] const char* intersection_type_to_string(IntersectionType type);

/*=============================================================================
 *  Box-level predicates (`intersectLevelBoundingBoxes.inc`)
 *============================================================================*/

/// Whether two stems' oriented rectangles overlap. Mirrors `intersectStemStem`
/// (`intersectLevelBoundingBoxes.inc:375`): brute-force segment intersection
/// of each stem's two long (`AB`/`CD`) sides against the other's.
[[nodiscard]] bool intersect_stem_stem(const StemBox& stem1, const StemBox& stem2);

/// Whether two loops' bounding circles overlap. Mirrors `intersectLoopLoop`
/// (`intersectLevelBoundingBoxes.inc:94`); each radius is grown by
/// `0.5 * epsilon_recognize(clearance)` before the circle-circle test.
[[nodiscard]] bool intersect_loop_loop(const LoopBox& loop1, const LoopBox& loop2,
                                       double clearance);

/// Whether @p loop's bounding circle (grown by `epsilon_recognize(clearance)`)
/// contains the point on @p stem's rectangle closest to @p loop's center.
/// Mirrors `intersectStemLoop` (`intersectLevelBoundingBoxes.inc:353`).
[[nodiscard]] bool intersect_stem_loop(const StemBox& stem, const LoopBox& loop, double clearance);

/// The bulge index a `*x_bulges`-family predicate found intersecting, or
/// `-1` if none did -- mirrors those predicates' `int *bulge` out-param
/// (`intersectLevelBoundingBoxes.inc`).
struct BulgeHit {
  bool intersects = false;
  int bulge = -1;
};

/// Whether any of @p stem's bulge triangles intersects @p loop's bounding
/// circle (grown by `epsilon_recognize(clearance)`). Mirrors
/// `intersectLoopBulges` (`intersectLevelBoundingBoxes.inc:523`); returns the
/// FIRST intersecting bulge found (matching the vendored early-return loop).
[[nodiscard]] BulgeHit intersect_loop_bulges(const LoopBox& loop, const StemBox& stem,
                                             double clearance);

/// The pair of bulge indices (one per stem) a `*x_bulges`-family predicate
/// found intersecting, or `-1`/`-1` if none did.
struct BulgeBulgeHit {
  bool intersects = false;
  int bulge1 = -1;
  int bulge2 = -1;
};

/// Whether any bulge triangle of @p stem1 intersects any bulge triangle of
/// @p stem2 (each bulge's edges grown `0.5 * epsilon_recognize(clearance)`
/// further out). Mirrors `intersectBulgesBulges`
/// (`intersectLevelBoundingBoxes.inc:551`); returns the FIRST intersecting
/// pair found.
[[nodiscard]] BulgeBulgeHit intersect_bulges_bulges(const StemBox& stem1, const StemBox& stem2,
                                                    double clearance);

/// Whether any bulge triangle of @p stem2 crosses one of @p stem1's two long
/// sides (each bulge's edges grown `epsilon_recognize(clearance)` further
/// out). Mirrors `intersectStemBulges` (`intersectLevelBoundingBoxes.inc:586`).
[[nodiscard]] BulgeHit intersect_stem_bulges(const StemBox& stem1, const StemBox& stem2,
                                             double clearance);

/*=============================================================================
 *  Node-level predicates (`intersectLevelTreeNodes.inc`)
 *============================================================================*/

/// The classification + (if applicable) intersecting bulge indices
/// `intersect_node_node` found; mirrors `intersectNodeNode`'s
/// `intersectionType` return plus its `int *bulge1`/`int *bulge2` out-params.
struct NodeIntersection {
  IntersectionType type = IntersectionType::none;
  int bulge1 = -1;
  int bulge2 = -1;
};

/**
 * Classify the intersection (if any) between two DISTINCT, non-root
 * `TreeNode`s. Mirrors `intersectNodeNode` (`intersectLevelTreeNodes.inc:175`):
 * an AABB early-reject, then the SxS/LxL/SxL/LxS/LxB/BxL/SxB/BxS/BxB checks
 * IN THAT EXACT ORDER (each gated by whether `node1`/`node2` is the other's
 * parent, since adjacent stem/loop pairs in a valid config never truly
 * intersect and would otherwise false-positive), returning the FIRST type
 * that matches.
 *
 * @param node1 First node; must have `sbox`/`lbox`/`aabb` set (i.e. not the
 *     tree's true exterior root -- see `intersect_node_exterior` for that
 *     case).
 * @param node2 Second node; same requirement.
 * @param clearance `PuzzlerOptions::clearance`.
 */
[[nodiscard]] NodeIntersection intersect_node_node(const TreeNode& node1, const TreeNode& node2,
                                                   double clearance);

/**
 * Whether @p node -- a DIRECT CHILD of the exterior root -- dips below the
 * exterior baseline (`geom::kExteriorY`). Mirrors `intersectNodeExterior`
 * (`intersectLevelTreeNodes.inc:88`): `false` for the root itself or any
 * node whose parent is not the root, and `false` whenever
 * @p check_exterior_intersections is `false` (mirroring the vendored
 * `puzzler->checkExteriorIntersections` gate).
 */
[[nodiscard]] bool intersect_node_exterior(const TreeNode& node, bool check_exterior_intersections,
                                           double clearance);

/**
 * The first node in @p tree's subtree (including @p tree itself) that
 * intersects @p node, or `nullptr` if none does. Mirrors `intersectNodeTree`
 * (`intersectLevelTreeNodes.inc:269`): checks @p tree itself before
 * recursing into its children, left to right.
 */
[[nodiscard]] const TreeNode* intersect_node_tree(const TreeNode& node, const TreeNode& tree,
                                                  double clearance);

/**
 * Whether ANY node of @p tree1's subtree intersects ANY node of @p tree2's
 * subtree (early-exit). Mirrors `intersectTrees`/`intersect_iterateTree`
 * (`intersectLevelTreeNodes.inc:293-335`): iterates @p tree1's subtree
 * (pre-order), calling `intersect_node_tree` against the whole of @p tree2
 * for each.
 */
[[nodiscard]] bool intersect_trees(const TreeNode& tree1, const TreeNode& tree2, double clearance);

/**
 * Whether any node in @p list1 intersects any node in @p list2 (early-exit).
 * Mirrors `intersectNodeLists` (`intersectLevelTreeNodes.inc:338`): a node
 * that IS the exterior root is compared via `intersect_node_exterior`
 * against its counterpart (which must then be a direct root child); two
 * non-root nodes are compared via `intersect_node_node`.
 */
[[nodiscard]] bool intersect_node_lists(const std::vector<const TreeNode*>& list1,
                                        const std::vector<const TreeNode*>& list2,
                                        bool check_exterior_intersections, double clearance);

/**
 * `check_optimize_intersections`'s (`optimize.cpp:54`) EXACT-parity, faster
 * replacement for `intersect_node_lists(subtree, subtree, ...) ||
 * intersect_node_lists(subtree, ancestor_list, ...)` -- SPEED lever A1,
 * `.claude/plans/current-plan-speed.md`. Implemented in `broad_phase.cpp`.
 *
 * MECHANISM: for `subtree.size() + ancestor_list.size()` above a small
 * brute-force threshold, this builds a uniform-grid broad-phase index over
 * every non-exterior node's `Aabb`, expanded by a per-call, per-node
 * CONSERVATIVE upper bound on `intersect_nodes_bounding_boxes`'s
 * (`intersect_tree.cpp:20`) `extra_distance` (`epsilon_recognize(clearance)
 * + max_bulge_dist`, where `max_bulge_dist` is the largest `StemBox::
 * bulge_dist` among every node this call considers). Only pairs whose
 * expanded boxes overlap -- a PROVABLE SUPERSET of every pair
 * `intersect_nodes_bounding_boxes` would itself accept -- are tested with
 * the exact, unmodified `intersect_node_node`; ancestor-vs-ancestor pairs
 * (never tested by the original two `intersect_node_lists` calls) are never
 * generated. The candidate set differs from the brute-force O(m^2) scan,
 * but the OR-over-candidate-pairs boolean it computes does not: every
 * non-candidate pair is provably `none`, so the returned boolean -- and
 * therefore every coordinate downstream of it -- is BIT-IDENTICAL. Below
 * the threshold, this literally calls `intersect_node_lists` (the original,
 * unmodified function), so the small-subtree case is trivially exact too.
 *
 * @param subtree `optimize_tree`'s collected subtree-node list.
 * @param ancestor_list `optimize_tree`'s collected ancestor-chain list.
 * @param check_exterior_intersections `PuzzlerOptions::check_exterior`.
 * @param clearance `PuzzlerOptions::clearance`.
 */
[[nodiscard]] bool any_intersection(const std::vector<const TreeNode*>& subtree,
                                    const std::vector<const TreeNode*>& ancestor_list,
                                    bool check_exterior_intersections, double clearance);

/*=============================================================================
 *  Detection set (parity-only; no vendored counterpart -- see file header)
 *============================================================================*/

/// One entry of `detect_intersections`'s result: an intersecting node pair
/// (identified by the DFS pre-order `id` `debug_dump.hpp`'s `dump_tree`
/// assigns -- `node2_id == 0` (the root's id) marks an `exterior` detection,
/// where `node1_id` is the exterior-intersecting child, not a true pair).
struct Detection {
  int node1_id = -1;
  int node2_id = -1;
  IntersectionType type = IntersectionType::none;
};

/**
 * The FULL intersection detection set over @p root's whole tree (T0: an
 * already `update_bounding_boxes`-d tree, pre-resolver): every intersecting
 * NON-ROOT node pair (via `intersect_node_node`, over all `n*(n-1)/2` pairs
 * in `dump_tree`'s DFS pre-order id numbering), plus every direct root child
 * that intersects the exterior baseline (via `intersect_node_exterior`,
 * always checked -- i.e. as if `check_exterior_intersections = true`).
 *
 * Order: node-pair detections first (in `(id_i, id_j)` nested-loop
 * discovery order, `id_i < id_j`), then exterior detections (in id order) --
 * deterministic, matching `vendor_instrument.c`'s
 * `rnadraw_oracle_dump_detections` construction so both sides can be
 * compared as ordered sequences, not just as sets.
 *
 * @param root The tree's root (exterior loop; `parent == nullptr`).
 * @param clearance `PuzzlerOptions::clearance`.
 */
[[nodiscard]] std::vector<Detection> detect_intersections(const TreeNode& root, double clearance);

}  // namespace rna_layout
