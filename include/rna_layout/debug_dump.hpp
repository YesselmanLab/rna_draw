#pragma once

/**
 * @file debug_dump.hpp
 * @brief Native-side parity dumps (`.claude/plans/current-plan.md`'s
 *        "Parity oracle harness" instrumentation seams), PARITY-ONLY: never
 *        called by any production layout path.
 *
 * `dump_tree` is the T0/T1 config-tree dump (criterion 2): a flat,
 * DFS-pre-order list of every node's topology + `Config` + boxes, compared
 * against the vendored oracle's own dump
 * (`src/vienna_layout/vendor_instrument.c`) from Python
 * (`tests/test_native_parity.py`).
 */

#include <optional>
#include <vector>

#include "rna_layout/tree.hpp"

namespace rna_layout {

/// Mirrors `ConfigArc` (`types.hpp`); duplicated as a dump-only value type
/// so `debug_dump.hpp` has no pybind11/JSON dependency of its own -- the
/// binding layer (`bindings.cpp`) converts this to a Python object.
struct DumpConfigArc {
  int segments = 0;
  double angle = 0.0;
};

struct DumpConfig {
  double radius = 0.0;
  double min_radius = 0.0;
  double default_radius = 0.0;
  std::vector<DumpConfigArc> arcs;
};

struct DumpLoopBox {
  double cx = 0.0;
  double cy = 0.0;
  double r = 0.0;
};

struct DumpStemBox {
  double ax = 0.0;
  double ay = 0.0;
  double bx = 0.0;
  double by = 0.0;
  double cx = 0.0;
  double cy = 0.0;
  double ex = 0.0;
  double ey = 0.0;
  int bulge_count = 0;
  double bulge_dist = 0.0;
};

/// One flattened `TreeNode`, DFS-pre-order. `id` is this entry's own index
/// in the list `dump_tree` returns, which coincides with the vendored
/// `treeNode.id` -- see `dump_tree`'s doc comment for why.
struct DumpTreeNode {
  int id = -1;
  int parent_id = -1;
  int loop_start = -1;
  int stem_start = -1;
  std::optional<DumpConfig> cfg;
  std::optional<DumpLoopBox> lbox;
  std::optional<DumpStemBox> sbox;
};

/**
 * Flatten the subtree rooted at @p root into a DFS-pre-order list, one
 * entry per node.
 *
 * The vendored `treeNode.id` is assigned by a single counter incremented
 * once per node in exactly this pre-order (`treeHandleStem`'s `++(*nodeID)`
 * fires immediately before its recursive `treeHandleLoop` call,
 * `configtree.inc:814-821`) -- so a fresh DFS pre-order numbering of an
 * INDEPENDENTLY built native tree reproduces the same id sequence as long
 * as both trees discover children in the same left-to-right order, which
 * `build_config_tree` does (same `pair_table` stem-discovery loop). This
 * lets the T0/T1 parity gate compare topology by id/parent_id without the
 * native `TreeNode` needing its own `id` field (`tree.hpp` intentionally
 * has none -- see the plan's value-type list).
 *
 * @param root The tree (or subtree) to flatten.
 * @return One `DumpTreeNode` per node in @p root's subtree, in DFS
 *     pre-order.
 */
[[nodiscard]] std::vector<DumpTreeNode> dump_tree(const TreeNode& root);

}  // namespace rna_layout
