#pragma once

/**
 * @file resolve.hpp
 * @brief The intersection-resolver driver, ported from
 *        `resolveIntersections.inc:24 checkAndFixIntersections` (Milestone A
 *        step 7: the SIBLING-intersection vertical only).
 *
 * SCOPE (Milestone A step 7): only the SIBLING-resolution path is real here
 * -- `check_and_fix_intersections` requires `opts.check_ancestor == false`
 * and `opts.optimize == false` (throws otherwise); `opts.check_sibling` may
 * be `true`. The ANCESTOR path (`handleAncestorIntersections.inc`, Milestone
 * A step 8) and `optimizeTree` (`optimize.inc`, Milestone A step 9) are not
 * ported -- their gates are preserved in the driver's control flow (matching
 * the reference's exact shape) but throw if ever reached, rather than being
 * silently skipped or approximated. See `puzzler.cpp`'s call site for how
 * `PuzzlerOptions` is validated before this ever runs.
 *
 * MUTATION MODEL: unlike every earlier Milestone A step (pure construction/
 * detection over an already-built tree), this driver MUTATES the tree in
 * place -- a resolved sibling intersection changes a node's `Config` (radius,
 * arc angles) and recomputes that node's subtree's boxes
 * (`apply_changes_to_config_and_bounding_boxes`, `config_tree.hpp`). No node
 * is ever reparented, added, or deleted (`TreeNode`'s `unique_ptr` ownership,
 * `tree.hpp`, needs no change for this) -- only `cfg`/`lbox`/`sbox`/`aabb`
 * fields already declared mutable-by-value are written through the SAME
 * `TreeNode*` addresses `puzzler.cpp`'s `build_config_tree` returned, which
 * is why this driver takes raw `TreeNode*` (never re-allocating), matching
 * the reference's own in-place-mutation model over `treeNode*`.
 */

#include <vector>

#include "rna_layout/intersect_tree.hpp"
#include "rna_layout/tree.hpp"
#include "rna_layout/types.hpp"

namespace rna_layout {

/**
 * One config-change DECISION the resolver made, in the order it made them --
 * the change-trace parity seam (`.claude/plans/current-plan.md`'s
 * "Parity oracle harness", criterion 4b): instruments
 * `check_and_apply_config_changes` (`config_changes.cpp`), the native
 * counterpart of the vendored `checkAndApplyConfigChanges`
 * (`handleConfigChanges.inc:55`), which `vendor_instrument.c`
 * macro-interposes on the oracle side so both traces can be compared.
 */
struct ChangeTraceEntry {
  /// `TreeNode::id` of the node whose `Config` this entry changed (or tried
  /// to).
  int node_id = -1;
  /// The context that triggered this change attempt; `siblings` for every
  /// entry this step's resolver produces (`IntersectionType::siblings`,
  /// mirroring the vendored `"BRA"` tag) -- other values are reserved for
  /// the ancestor path (Milestone A step 8).
  IntersectionType type = IntersectionType::none;
  /// Per-arc angle deltas (radians), one entry per the node's `Config::arcs`,
  /// AFTER `check_and_apply_config_changes`'s "fix too small changes"
  /// adjustment (`handleConfigChanges.inc:76-98`) -- the values actually
  /// tested for validity and (if `accepted`) applied.
  std::vector<double> deltas;
  /// Whether `cfg_is_valid` accepted these deltas (and they were applied) --
  /// mirrors `checkAndApplyConfigChanges`'s `return 1`/`return 0`.
  bool accepted = false;
};

/**
 * Resolver-wide mutable state, threaded by reference through every resolver
 * function -- replaces the vendored `vrna_plot_options_puzzler_t`'s two
 * mutable fields (`numberOfChangesAppliedToConfig`,
 * `maximumNumberOfConfigChangesAllowed`) that live OUTSIDE `PuzzlerOptions`
 * here (`PuzzlerOptions` stays a pure, unmutated set of levers -- see
 * `types.hpp`).
 */
struct ResolverState {
  /// Running count of `check_and_apply_config_changes` calls (accepted or
  /// not) -- mirrors `numberOfChangesAppliedToConfig`.
  int changes_applied = 0;
  /// The effective config-change budget for this resolve pass: the caller
  /// (`puzzler.cpp`) applies the vendored `<= 0 -> 25000` fallback
  /// (`RNApuzzler.c:460-461`) once, up front, and stores the result here --
  /// mirrors `maximumNumberOfConfigChangesAllowed`.
  int max_config_changes = 25000;
  /// The ordered change trace (see `ChangeTraceEntry`); empty unless a
  /// caller wants it (`puzzler.cpp`'s production call site does not read
  /// it -- only the parity-instrumentation entry points in `bindings.cpp`
  /// do).
  std::vector<ChangeTraceEntry> trace;
};

/**
 * Resolve intersections in @p node's subtree in place, restarting its own
 * traversal as needed. Ported from `checkAndFixIntersections`
 * (`resolveIntersections.inc:24`); see this file's header for scope.
 *
 * @param node The (sub)tree to resolve; mutated in place.
 * @param opts Resolver options; `check_ancestor` and `optimize` MUST both be
 *     `false` (see this file's header).
 * @param state Resolver-wide mutable state (change counter + trace),
 *     threaded through and updated in place.
 * @return `nullptr` always, in this step's scope (the non-null "propagate
 *     to an ancestor" return is Milestone A step 8's -- see this file's
 *     header).
 * @throws std::logic_error If @p opts requests the not-yet-ported ancestor
 *     or optimize passes.
 */
TreeNode* check_and_fix_intersections(TreeNode* node, const PuzzlerOptions& opts,
                                      ResolverState& state);

}  // namespace rna_layout
