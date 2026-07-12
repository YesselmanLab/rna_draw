/**
 * @file resolve_siblings.cpp
 * @brief `check_siblings` + its two private helpers, ported from
 *        `handleSiblingIntersections.inc`: resolve intersections between
 *        `node`'s DIRECT children by widening the angular space between
 *        the intersecting pair (via `calc_deltas`) and applying the result
 *        (via `check_and_apply_config_changes`), one pair at a time, until
 *        the first pair that actually changes something.
 */

#include <cmath>
#include <optional>
#include <utility>

#include "resolve_internal.hpp"
#include "rna_layout/bounding_wedge.hpp"
#include "rna_layout/geometry.hpp"

namespace rna_layout {

namespace {

/// Try to fix the intersection between @p node's children @p left and
/// @p right by widening the wedge between them (`calc_deltas`) and applying
/// the result (`check_and_apply_config_changes`). Ported from
/// `fixIntersectionOfSiblings` (`handleSiblingIntersections.inc:48`).
///
/// @param node Common ancestor (parent) of the intersecting @p left/@p right
///     children.
/// @param delta_cfg Scratch buffer, reused (not re-zeroed) across every pair
///     `check_siblings` tries for @p node -- see `calc_deltas`'s doc
///     comment: it fully overwrites this buffer whenever it actually
///     computes anything, so stale content from a prior pair is never read.
bool fix_intersection_of_siblings(TreeNode& node, int left, int right,
                                  std::vector<double>& delta_cfg, const PuzzlerOptions& opts,
                                  ResolverState& state) {
  const double min_angle = bounding_wedge(node, right, opts.clearance).min_angle;
  const double max_angle = bounding_wedge(node, left, opts.clearance).max_angle;
  double target_angle = min_angle - max_angle;

  bool changed = false;
  if (target_angle < 0) {
    // Limit the angle to avoid malformed structures.
    target_angle = std::fmax(target_angle, -geom::kPiHalf);
    const double changed_angle = calc_deltas(node, node.parent, left, right, -target_angle,
                                             opts.paired, opts.clearance, delta_cfg);

    if (changed_angle != 0.0) {
      changed = check_and_apply_config_changes(node, delta_cfg, IntersectionType::siblings,
                                               opts.unpaired, opts.paired, state);
    }
  }
  return changed;
}

/// Resolve @p intersections (pairs of @p node's child indices) one at a
/// time, stopping at the first one that actually changes something. Ported
/// from `handleIntersectionOfSiblings` (`handleSiblingIntersections.inc:99`).
///
/// @return `std::nullopt` if `state.changes_applied` exceeds
///     `state.max_config_changes` (mirrors the vendored `return -1`);
///     otherwise whether any pair was actually fixed.
std::optional<bool> handle_intersection_of_siblings(
    TreeNode& node, const std::vector<std::pair<int, int>>& intersections,
    const PuzzlerOptions& opts, ResolverState& state) {
  if (state.changes_applied > state.max_config_changes) {
    return std::nullopt;
  }

  bool changed = false;
  const int config_size = static_cast<int>(node.children.size()) + 1;
  std::vector<double> delta_cfg(static_cast<std::size_t>(config_size), 0.0);

  for (const auto& [left, right] : intersections) {
    changed = fix_intersection_of_siblings(node, left, right, delta_cfg, opts, state);
    if (changed) {
      break;
    }
  }

  return changed;
}

}  // namespace

SiblingCheckResult check_siblings(TreeNode& node, const PuzzlerOptions& opts,
                                  ResolverState& state) {
  const int child_count = static_cast<int>(node.children.size());

  std::vector<std::pair<int, int>> intersections;
  for (int i = 0; i < child_count; ++i) {
    for (int j = i + 1; j < child_count; ++j) {
      if (intersect_trees(*node.children[static_cast<std::size_t>(i)],
                          *node.children[static_cast<std::size_t>(j)], opts.clearance)) {
        intersections.emplace_back(i, j);
      }
    }
  }

  if (intersections.empty()) {
    return SiblingCheckResult::none;
  }

  const std::optional<bool> outcome =
      handle_intersection_of_siblings(node, intersections, opts, state);
  if (!outcome.has_value()) {
    return SiblingCheckResult::aborted;
  }

  // DOCUMENTED DEVIATION FROM THE VENDORED BIT-FLAG SEMANTICS (see the
  // vendored `checkSiblings`, `handleSiblingIntersections.inc:149-223`:
  // `ret = _intersect | (retFix ? _changed : 0)`, TRUTHY -- i.e. "restart"
  // -- whenever ANY intersecting pair was found, whether or not
  // `handle_intersection_of_siblings` actually changed anything).
  //
  // Verified (this port's parity harness, hard_set): on a nontrivial
  // fraction of real structures, some intersecting sibling pair's bounding
  // wedges never overlap in angle-space (`target_angle >= 0` in
  // `fix_intersection_of_siblings`, every time), so NO config change is
  // ever attempted for it. Restarting on that outcome re-runs this exact
  // deterministic check against a byte-for-byte UNCHANGED tree (nothing
  // this call touched any `Config`/box when `outcome == false`) --
  // producing the identical `restart` result FOREVER. This is not a
  // porting bug: the VENDORED oracle hangs identically, on the identical
  // structures, under the equivalent options
  // (`checkSiblingIntersections=1, checkAncestorIntersections=0,
  // optimize=0` -- `plot_coords_puzzler_sibling_only`,
  // `src/vienna_layout/bindings.cpp`).
  //
  // RE-DERIVED AT MILESTONE A STEP 8 (ancestor resolution now real): this
  // shortcut REMAINS behavior-preserving with `check_ancestor == true`. See
  // `resolve_internal.hpp`'s `SiblingCheckResult` doc comment for the full
  // argument -- in short, `outcome == false` means NO `Config`/box anywhere
  // in the tree changed (not just at this node), so a restart's ancestor-
  // check/children-recursion/sibling-check all reproduce their identical
  // prior answers regardless of whether ancestor-checking is enabled; a real
  // restart driven by an ancestor fix flows through `resolve.cpp`'s
  // SEPARATE `checkTree = true` paths, never through this function's return
  // value. Empirically re-verified too: re-running the vendored oracle's
  // sibling-only hang list under `checkSibling=1, checkAncestor=1,
  // optimize=0` -- see `tests/test_native_parity.py`'s sibling+ancestor
  // suite and this port's handoff report for the (possibly narrower)
  // re-derived exclusion set.
  if (!outcome.value()) {
    return SiblingCheckResult::none;
  }
  return SiblingCheckResult::restart;
}

}  // namespace rna_layout
