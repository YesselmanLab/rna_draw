/**
 * @file resolve.cpp
 * @brief Implements `resolve.hpp`'s `check_and_fix_intersections`, ported
 *        from `checkAndFixIntersections` (`resolveIntersections.inc:24`).
 */

#include "rna_layout/resolve.hpp"

#include "resolve_internal.hpp"
#include "rna_layout/config_tree.hpp"

namespace rna_layout {

namespace {

/// Whether `check_and_fix_intersections` should invoke `optimize_tree` at
/// @p node once its sibling/ancestor resolution has settled. Ported from
/// `checkAndFixIntersections`'s "----- OPTIMIZATIONS -----" gate
/// (`resolveIntersections.inc:120-136`): the exterior root never
/// optimizes; a node whose PARENT is the exterior root always does (its
/// loop was never itself intersection-tested against an ancestor); any
/// other node only optimizes once its own resolved radius has grown to
/// more than 10x its un-widened default -- a cheap proxy for "this loop's
/// resolution likely left slack `optimize_tree` can shrink back out".
bool should_optimize_node(const TreeNode& node) {
  if (is_exterior(node)) {
    return false;
  }
  if (is_exterior(*node.parent)) {
    return true;
  }
  const Config& cfg = *node.cfg;  // NOLINT(bugprone-unchecked-optional-access)
  return cfg.radius > 10 * cfg.default_radius;
}

}  // namespace

TreeNode* check_and_fix_intersections(TreeNode* node, const PuzzlerOptions& opts,
                                      ResolverState& state) {
  bool check_tree = true;
  while (check_tree) {
    check_tree = false;

    // On the way from root to leaves: resolve ancestor intersections
    // (Milestone A step 8). Ported from `resolveIntersections.inc:52-60`: if
    // `node` itself has an ancestor (or exterior) intersection that was
    // fixed by rotating some node ON THE PATH BETWEEN THEM, propagate that
    // node straight back up the RECURSION STACK (not the tree) -- see the
    // "propagation" note below, at this function's recursive-call site,
    // for how each enclosing level decides whether to keep propagating or
    // restart at its own level.
    if (opts.check_ancestor && !is_exterior(*node)) {
      TreeNode* changed_by_ancestor = check_node_against_ancestors(*node, opts, state);
      if (changed_by_ancestor != nullptr) {
        return changed_by_ancestor;
      }
    }

    // Recursive call for all children.
    if (!check_tree) {
      for (auto& child : node->children) {
        TreeNode* changed_by_recursion = check_and_fix_intersections(child.get(), opts, state);
        if (changed_by_recursion != nullptr) {
          // PROPAGATION (`resolveIntersections.inc:76-88`): a non-null
          // return only ever comes from the ancestor branch above, applied
          // somewhere at or below `child`, and that call's `rotationNode`
          // is ALWAYS a genuine tree-ancestor of the node it was checking
          // (`fix_intersection_with_ancestor`'s caller only ever offers
          // nodes strictly between an intersector and its own ancestor) --
          // so `changed_by_recursion`'s id is always `<= node->id` here,
          // with equality iff it IS `node`. There is no third case (the
          // vendored `if`/`else if` has no trailing `else` either).
          if (changed_by_recursion->id < node->id) {
            // The rotated node is a STRICT ancestor of `node` too: keep
            // propagating up unchanged (this level's own state is
            // untouched by that rotation, so there is nothing to restart
            // here).
            return changed_by_recursion;
          }
          if (changed_by_recursion == node) {
            // The rotated node IS `node`: its own Config/boxes changed, so
            // restart this level's `while` loop (re-run the ancestor check,
            // children recursion, and sibling check against the new state).
            check_tree = true;
            break;
          }
        }
        // `changed_by_recursion == nullptr`: that child's whole subtree
        // (and its own ancestor chain up to and including `node`) is
        // already intersection-free; move on to the next child.
      }
    }

    // On the way back from leaves to root: resolve sibling intersections.
    if (opts.check_sibling && !is_exterior(*node) && !check_tree) {
      const SiblingCheckResult result = check_siblings(*node, opts, state);
      if (result == SiblingCheckResult::aborted) {
        return nullptr;
      }
      if (result == SiblingCheckResult::restart) {
        check_tree = true;
        continue;
      }
    }
  }

  // ----- OPTIMIZATIONS ----- (Milestone A step 9). Ported from
  // `resolveIntersections.inc:120-140`: unlike the ancestor/sibling checks
  // above (each re-run every time this function's own `while` loop
  // restarts), this gate + `optimize_tree` call happens EXACTLY ONCE per
  // node, at the very end of this function's call for that node -- i.e.
  // once per node across the WHOLE RECURSION, not just at the tree's true
  // root.
  if (opts.optimize && should_optimize_node(*node)) {
    // The returned shrinking ratio is informational only (mirrors
    // `optimizeTree`'s vendored caller, which also discards it here --
    // `resolveIntersections.inc:136`); this driver only cares about the
    // tree it mutates in place.
    (void)optimize_tree(*node, opts, state);
  }

  return nullptr;
}

}  // namespace rna_layout
