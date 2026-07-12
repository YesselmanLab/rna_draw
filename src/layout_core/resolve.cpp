/**
 * @file resolve.cpp
 * @brief Implements `resolve.hpp`'s `check_and_fix_intersections`, ported
 *        from `checkAndFixIntersections` (`resolveIntersections.inc:24`).
 */

#include "rna_layout/resolve.hpp"

#include <stdexcept>

#include "resolve_internal.hpp"
#include "rna_layout/config_tree.hpp"

namespace rna_layout {

TreeNode* check_and_fix_intersections(TreeNode* node, const PuzzlerOptions& opts,
                                      ResolverState& state) {
  bool check_tree = true;
  while (check_tree) {
    check_tree = false;

    // On the way from root to leaves: resolve ancestor intersections
    // (Milestone A step 8, NOT YET PORTED). `layout_puzzler` throws before
    // this driver ever runs with `opts.check_ancestor == true`
    // (`puzzler.cpp`), so this branch's condition mirrors
    // `resolveIntersections.inc:52-60`'s gate faithfully but is PROVABLY
    // unreachable in this build -- the throw below is a defensive
    // trip-wire, not expected behavior.
    if (opts.check_ancestor && !is_exterior(*node)) {
      throw std::logic_error(
          "rna_layout::check_and_fix_intersections: ancestor intersection "
          "resolution is not implemented yet (Milestone A step 8)");
    }

    // Recursive call for all children.
    if (!check_tree) {
      for (auto& child : node->children) {
        TreeNode* changed_by_recursion = check_and_fix_intersections(child.get(), opts, state);
        if (changed_by_recursion != nullptr) {
          // `checkNodeAgainstAncestors` (the ancestor branch above) is the
          // ONLY producer of a non-null return in the reference, and it is
          // never called while `opts.check_ancestor == false` (see above) --
          // so this is unreachable too, for the same reason.
          throw std::logic_error(
              "rna_layout::check_and_fix_intersections: unexpected non-null "
              "ancestor-propagation return (Milestone A step 8)");
        }
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

  // ----- OPTIMIZATIONS ----- (Milestone A step 9, NOT YET PORTED).
  // `layout_puzzler` throws before this driver ever runs with
  // `opts.optimize == true` (`puzzler.cpp`); same defensive trip-wire as
  // above.
  if (opts.optimize) {
    throw std::logic_error(
        "rna_layout::check_and_fix_intersections: optimize is not "
        "implemented yet (Milestone A step 9)");
  }

  return nullptr;
}

}  // namespace rna_layout
