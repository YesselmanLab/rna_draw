#pragma once

/**
 * @file tree.hpp
 * @brief RAII replacement for the vendored `treeNode` (`configtree_struct.h`).
 *
 * `TreeNode` is foundation laid down for the config-tree/bounding-box/resolver
 * steps (Milestone A, steps 4+); the turtle-base pass (this slice) does not
 * build a tree at all -- `vrna_plot_coords_turtle_pt` returns before
 * `buildConfigtree` is ever called (`RNAturtle.c`, `coordinates.inc`). It is
 * included now so those later steps do not need a second foundation pass, and
 * so `unique_ptr` ownership (replacing the vendored `freeTree`/`freeBulges`)
 * is decided once, up front.
 */

#include <memory>
#include <optional>
#include <vector>

#include "rna_layout/types.hpp"

namespace rna_layout {

/**
 * One node of the config tree: one loop (or the exterior loop) plus the
 * stem that encloses it.
 *
 * Ownership mirrors the vendored `treeNode`'s parent/children shape, but
 * RAII replaces `configtree.inc`'s manual `freeTree`: destroying a
 * `TreeNode` recursively destroys its children via `std::unique_ptr`.
 * `parent` is intentionally a non-owning raw pointer (the parent already
 * owns this node through `children`); it is `nullptr` for the root.
 */
struct TreeNode {
  TreeNode* parent = nullptr;
  std::vector<std::unique_ptr<TreeNode>> children;

  std::optional<Config> cfg;
  int loop_start = -1;
  int stem_start = -1;

  std::optional<LoopBox> lbox;
  std::optional<StemBox> sbox;
  Aabb aabb;

  /**
   * Append a freshly constructed child, wire its `parent` back to `this`,
   * and return a non-owning pointer to it (the canonical way callers keep
   * building the tree without fighting `unique_ptr` ownership).
   *
   * @return Non-owning pointer to the newly appended child.
   */
  TreeNode* add_child() {
    children.push_back(std::make_unique<TreeNode>());
    TreeNode* child = children.back().get();
    child->parent = this;
    return child;
  }
};

}  // namespace rna_layout
