#pragma once

/**
 * @file resolve_internal.hpp
 * @brief Private declarations shared across the resolver's translation
 *        units (`deltas.cpp`, `config_changes.cpp`, `resolve_siblings.cpp`,
 *        `resolve.cpp`) -- the `*_impl.hpp` convention
 *        `.claude/plans/current-plan.md`'s directory layout note describes,
 *        never installed under `include/rna_layout`.
 *
 * Every vendored counterpart declared here was `PRIVATE` (a `static`
 * function local to the `.inc` amalgam, `resolveIntersections.inc`'s
 * `#include` chain); splitting the port across several `.cpp` files (one
 * per plan's module-size guidance) needs SOME shared declaration point for
 * those same cross-file calls, which this file is.
 */

#include <cstdint>
#include <vector>

#include "rna_layout/config_tree.hpp"
#include "rna_layout/intersect_tree.hpp"
#include "rna_layout/resolve.hpp"
#include "rna_layout/tree.hpp"
#include "rna_layout/types.hpp"

namespace rna_layout {

/**
 * The area between @p node's children @p index_left and @p index_right
 * (traversing clockwise starting at @p index_left) is enlarged by
 * @p delta_angle, radians; every other arc absorbs the increase by
 * shrinking. Ported from `calcDeltas` (`calcDeltas.inc:37`).
 *
 * NOTE: the vendored doc comment says its `deltas` out-param is "in degree
 * format" -- that is a vendored inaccuracy, not a unit conversion this port
 * is missing: `deltaAngle`'s only caller (`fixIntersectionOfSiblings`)
 * passes a radian `targetAngle`, and every downstream arithmetic op
 * (`asin`, `MATH_PI` comparisons) is radian-native throughout.
 *
 * @param node The loop whose arcs are being redistributed.
 * @param recursive_end Do not walk the "can this bend go to a higher tree
 *     level" search (see `calc_deltas.cpp`) past this ancestor; siblings
 *     always pass `node.parent` here, which makes that search provably a
 *     no-op for the sibling path (see `deltas.cpp`'s doc comment).
 * @param index_left,index_right Child indices bounding the area to enlarge.
 * @param delta_angle The desired increase, radians; `< 0` is rejected
 *     (returns `0.0` immediately).
 * @param paired Distance between the two bases of a base pair.
 * @param clearance `PuzzlerOptions::clearance`.
 * @param deltas Out-param: per-arc angle deltas, resized to
 *     `node.children.size() + 1` and fully overwritten.
 * @return The amount of @p delta_angle actually achieved by @p deltas
 *     (`<= delta_angle`).
 */
[[nodiscard]] double calc_deltas(const TreeNode& node, const TreeNode* recursive_end,
                                 int index_left, int index_right, double delta_angle, double paired,
                                 double clearance, std::vector<double>& deltas);

/**
 * Apply @p delta_cfg to @p tree's `Config` and recompute its subtree's boxes,
 * if @p delta_cfg passes `cfg_is_valid` (after a "fix too small changes"
 * adjustment done in place). Ported from `checkAndApplyConfigChanges`
 * (`handleConfigChanges.inc:55`); always increments
 * `state.changes_applied` and appends one `ChangeTraceEntry`, whether or not
 * the change was accepted.
 *
 * @param tree The node whose `Config` may change; must not be the root.
 * @param delta_cfg Per-arc angle deltas (radians), one per `tree.cfg->arcs`;
 *     mutated in place by the "fix too small changes" adjustment.
 * @param it The intersection context, recorded in the change trace (and,
 *     for the sibling path, always `IntersectionType::siblings`).
 * @param unpaired Default backbone-step distance.
 * @param paired Distance between the two bases of a base pair.
 * @param state Resolver-wide mutable state; updated in place.
 * @return Whether the change was accepted and applied.
 */
bool check_and_apply_config_changes(TreeNode& tree, std::vector<double>& delta_cfg,
                                    IntersectionType it, double unpaired, double paired,
                                    ResolverState& state);

/// `checkSiblings`'s outcome (`_false` / `_intersect[|_changed]` / the
/// `numberOfChangesAppliedToConfig` budget exceeded), modeled as an enum
/// here instead of the vendored bit-flag `short` -- `resolveIntersections
/// .inc:107-113`'s caller only ever asks "abort?" / "retry this node?" /
/// "nothing to do".
///
/// DOCUMENTED DEVIATION from the vendored bit-flag semantics
/// (`ret = _intersect | (retFix ? _changed : 0)`, truthy -- "restart" --
/// whenever ANY intersection was found, fixed or not): this port returns
/// `none` (not `restart`) when intersections were found but NOTHING was
/// actually fixed, because `layout_puzzler` GUARANTEES
/// `check_ancestor == optimize == false` for this whole port (Milestone A
/// step 7's scope), which makes an unproductive restart PROVABLY a no-op
/// (the tree is byte-for-byte unchanged, so re-checking is deterministic
/// and returns the identical answer) -- restarting anyway would hang
/// forever on any structure where a detected pair's bounding wedges never
/// overlap in angle-space (verified on `hard_set`: the VENDORED oracle
/// hangs identically, on the identical structures, under the equivalent
/// options -- see `resolve_siblings.cpp`'s `check_siblings` doc comment for
/// the full analysis). MUST BE REVISITED at Milestone A step 8, once
/// `check_ancestor` can mutate the tree between sibling-check attempts.
enum class SiblingCheckResult : std::uint8_t {
  /// No sibling of @p node's children intersects (or intersections were
  /// found but PROVABLY could not be fixed by this call -- see the enum's
  /// doc comment above); the driver does nothing further at this node.
  none,
  /// At least one sibling intersection was found AND actually fixed; the
  /// driver restarts its `while (check_tree)` loop.
  restart,
  /// `state.changes_applied` exceeded `state.max_config_changes` mid-fix;
  /// the driver aborts (`checkAndFixIntersections` returns `NULL`
  /// immediately, matching `resolveIntersections.inc:108-109`).
  aborted,
};

/**
 * Check and fix intersections between all of @p node's DIRECT children.
 * Ported from `checkSiblings` (`handleSiblingIntersections.inc:32`).
 *
 * @param node The node whose children's subtrees are checked pairwise.
 * @param opts Resolver options (`clearance`, `paired`/`unpaired`).
 * @param state Resolver-wide mutable state; updated in place.
 */
[[nodiscard]] SiblingCheckResult check_siblings(TreeNode& node, const PuzzlerOptions& opts,
                                                ResolverState& state);

}  // namespace rna_layout
