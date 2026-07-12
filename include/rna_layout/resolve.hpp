#pragma once

/**
 * @file resolve.hpp
 * @brief The intersection-resolver driver, ported from
 *        `resolveIntersections.inc:24 checkAndFixIntersections` (Milestone A
 *        step 7 landed the SIBLING path; step 8 the ANCESTOR path; step 9
 *        the OPTIMIZE pass).
 *
 * SCOPE (Milestone A step 9, COMPLETE): all three resolver stages are real
 * -- `check_and_fix_intersections` accepts any combination of
 * `opts.check_sibling`/`opts.check_ancestor`/`opts.optimize`. With every
 * flag `true` (`PuzzlerOptions`'s own defaults), this driver reproduces the
 * vendored `checkAndFixIntersections`'s full default (PRODUCTION) behavior
 * -- see `puzzler.cpp`'s call site.
 *
 * MUTATION MODEL: unlike every earlier Milestone A step (pure construction/
 * detection over an already-built tree), this driver MUTATES the tree in
 * place -- a resolved sibling or ancestor intersection changes a node's
 * `Config` (radius, arc angles) and recomputes that node's subtree's boxes
 * (`apply_changes_to_config_and_bounding_boxes`, `config_tree.hpp`). No node
 * is ever reparented, added, or deleted (`TreeNode`'s `unique_ptr` ownership,
 * `tree.hpp`, needs no change for this) -- only `cfg`/`lbox`/`sbox`/`aabb`
 * fields already declared mutable-by-value are written through the SAME
 * `TreeNode*` addresses `puzzler.cpp`'s `build_config_tree` returned, which
 * is why this driver takes raw `TreeNode*` (never re-allocating), matching
 * the reference's own in-place-mutation model over `treeNode*`.
 *
 * RETURN-VALUE PROPAGATION (new in step 8): a resolved ancestor intersection
 * can make this function return non-`nullptr` from a RECURSIVE call (not
 * just the top-level one) -- see `resolve.cpp`'s "PROPAGATION" comment at
 * the recursive-call site for the id-comparison argument that makes
 * re-dispatching that return (propagate further up vs. restart this level)
 * sound.
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
 * @param opts Resolver options; `check_ancestor`/`check_sibling`/`optimize`
 *     may each independently be `true` or `false`.
 * @param state Resolver-wide mutable state (change counter + trace),
 *     threaded through and updated in place.
 * @return Non-`nullptr` ONLY from an internal RECURSIVE call, meaning "an
 *     ancestor intersection was fixed by rotating a node ABOVE @p node's own
 *     tree position -- the caller (an ancestor frame) must handle it, not
 *     this one" (see `resolve.cpp`'s "PROPAGATION" comment). Calling this at
 *     the tree's true root (as `puzzler.cpp` does) always returns `nullptr`:
 *     a rotation candidate can never be the root itself (`is_interior_loop`/
 *     `is_multi_loop` are both `false` for the root by definition), so every
 *     non-null return is fully absorbed at some node strictly below the
 *     root before ever reaching this function's outermost call.
 */
TreeNode* check_and_fix_intersections(TreeNode* node, const PuzzlerOptions& opts,
                                      ResolverState& state);

}  // namespace rna_layout
