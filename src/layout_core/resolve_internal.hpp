#pragma once

/**
 * @file resolve_internal.hpp
 * @brief Private declarations shared across the resolver's translation
 *        units (`deltas.cpp`, `config_changes.cpp`, `resolve_siblings.cpp`,
 *        `resolve_ancestors.cpp`, `exterior_boxes.cpp`, `rotation_angle.cpp`,
 *        `rotation_angle2.cpp`, `resolve.cpp`) -- the `*_impl.hpp`
 *        convention `.claude/plans/current-plan.md`'s directory layout note
 *        describes, never installed under `include/rna_layout`.
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
/// actually fixed.
///
/// RE-DERIVED AT MILESTONE A STEP 8 (ancestor resolution now real -- see
/// `resolve.cpp`'s driver): this shortcut remains behavior-preserving even
/// with `check_ancestor == true`. The reasoning that "found but not fixed"
/// implies "the tree is byte-for-byte unchanged" is INDEPENDENT of whether
/// ancestor-checking is enabled -- `handle_intersection_of_siblings`
/// returning `false` means NO `check_and_apply_config_changes` call for
/// THIS node's subtree ever accepted a delta (see that function; a `false`
/// outcome is only reached when every candidate pair's `target_angle >= 0`,
/// so `check_and_apply_config_changes` is never even attempted), so no
/// `Config`/box anywhere in the tree actually changed. A restart of
/// `check_and_fix_intersections`'s `while` loop at THIS node therefore
/// re-runs (a) the ancestor-check on this exact, unchanged node against its
/// unchanged ancestors, (b) the children recursion over unchanged children,
/// and (c) this exact sibling check again -- all three are pure functions
/// of tree state, so all three reproduce their identical prior answer.
/// `restart` here would still hang forever on the same structures the
/// sibling-only case did; `none` is provably equivalent to `restart` for
/// this outcome, ancestor on or off. (What DOES change with ancestor on: a
/// restart triggered by the ancestor branch's OWN "found and fixed" path,
/// or by a child's ancestor-check reaching up and mutating this node's own
/// `Config` -- both flow through `resolve.cpp`'s separate `checkTree = true`
/// paths, never through this enum, and are real restarts with real new
/// state to re-check.) Verified empirically too: re-ran the vendored
/// oracle's sibling-only hang list under `checkSibling=1, checkAncestor=1,
/// optimize=0` -- see `resolve.cpp`'s driver doc comment and this port's
/// parity harness for the (possibly narrower) re-derived exclusion set.
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

/*=============================================================================
 *  Ancestor-intersection resolution (Milestone A step 8), ported from
 *  `handleAncestorIntersections.inc` (989 LOC) -- `resolve_ancestors.cpp`
 *  implements everything below except `get_rotation_angle` (`rotation_angle
 *  2.cpp`, dispatching to the per-intersection-type solvers in
 *  `rotation_angle.cpp`/`rotation_angle2.cpp`) and
 *  `setup_exterior_bounding_boxes` (`exterior_boxes.cpp`).
 *
 *  SCOPE NOTE (dead code NOT ported, verified by `grep` over the whole
 *  vendored call graph): `handleAncestorIntersections.inc` defines SIX
 *  functions never called by anything -- `getRotationSign` (the un-prefixed,
 *  3-line-shorter variant; only `TENTATIVE2_getRotationSign` below is ever
 *  called), `TENTATIVE3_getRotationSign`, `TENTATIVE_updateExteriorBounding
 *  Boxes`, `TENTATIVE2_setupExteriorBoundingBoxes`,
 *  `TENTATIVE_setupExteriorBoundingBoxes`, and the un-prefixed
 *  `setupExteriorBoundingBoxes` (only `TENTATIVE3_setupExteriorBoundingBoxes`
 *  below is ever called, despite its "TENTATIVE" name and the plan's
 *  `:754` line anchor pointing at the dead un-prefixed one) -- none of these
 *  six are ported; this port's function names below drop the vendored
 *  "TENTATIVE[23]_" prefixes since exactly one variant of each survives.
 *============================================================================*/

/// Whether @p node is an interior loop whose single child sits exactly
/// opposite its stem (`get_child_angle_by_index(node, 0) == geom::kPi`,
/// exact float equality preserved for fidelity). Mirrors
/// `isStraightInteriorLoop` (`handleAncestorIntersections.inc:297`); such a
/// loop can never usefully rotate (there is no "more room" on either side),
/// so `construct_reduced_intersection_path` skips it.
[[nodiscard]] bool is_straight_interior_loop(const TreeNode& node);

/// The path of tree nodes from @p intersector up to (and, for most
/// intersection types, including) @p ancestor, in root-to-leaf order,
/// SKIPPING any `is_straight_interior_loop` node in between (they cannot
/// usefully rotate -- see that function). Mirrors
/// `constructReducedIntersectionPath` (`handleAncestorIntersections.inc:311`).
///
/// @param ancestor The ancestor-intersection's ancestor endpoint.
/// @param intersector The ancestor-intersection's descendant endpoint;
///     `path.back() == &intersector` always.
/// @param it The intersection type (`Lx?` types exclude @p ancestor itself
///     from the returned path, unless it is itself a straight interior
///     loop that was already going to be excluded -- see the vendored
///     `switch` this ports verbatim).
[[nodiscard]] std::vector<TreeNode*> construct_reduced_intersection_path(TreeNode& ancestor,
                                                                          TreeNode& intersector,
                                                                          IntersectionType it);

/// The clockwise (`1`) / counter-clockwise (`-1`) / degenerate (`0`)
/// rotation sense of @p path (as `construct_reduced_intersection_path`
/// returns it): the signed sum of each consecutive pair's `get_child_angle`
/// minus `pi`. Mirrors the LIVE variant, `TENTATIVE2_getRotationSign`
/// (`handleAncestorIntersections.inc:101`) -- see this file's scope note for
/// why the un-prefixed `getRotationSign` and `TENTATIVE3_getRotationSign`
/// (both dead code) are not ported.
[[nodiscard]] short get_rotation_sign(const std::vector<TreeNode*>& path);

/**
 * Try to fix the ancestor-intersection between @p ancestor and
 * @p intersector by rotating @p rotation_node's children (a node somewhere
 * on the path between them). Ported from `fixIntersectionWithAncestor`
 * (`handleAncestorIntersections.inc:214`).
 *
 * @param ancestor The ancestor-intersection's ancestor endpoint.
 * @param rotation_node The candidate node to rotate (an interior or multi
 *     loop strictly between @p intersector and @p ancestor).
 * @param intersector The ancestor-intersection's descendant endpoint.
 * @param rotation_index The child index of @p rotation_node's PATH
 *     child (`get_child_index`), i.e. which side of @p rotation_node's loop
 *     the rotation should widen away from.
 * @param rotation_sign `get_rotation_sign`'s result for the whole path.
 * @param it The intersection type.
 * @param opts Resolver options.
 * @param state Resolver-wide mutable state; updated in place.
 * @return @p rotation_node if a config change was accepted and applied;
 *     `nullptr` otherwise.
 */
[[nodiscard]] TreeNode* fix_intersection_with_ancestor(TreeNode& ancestor, TreeNode& rotation_node,
                                                        TreeNode& intersector, int rotation_index,
                                                        short rotation_sign, IntersectionType it,
                                                        const PuzzlerOptions& opts,
                                                        ResolverState& state);

/**
 * Try to fix the ancestor-intersection between @p ancestor and
 * @p intersector by trying every candidate rotation node on the reduced
 * path between them (interior loops first, then multi loops, each nearest-
 * intersector-first), stopping at the first one that actually changes
 * something. Ported from `handleIntersectionWithAncestor`
 * (`handleAncestorIntersections.inc:361`).
 *
 * @return The node whose `Config` changed, or `nullptr` if @p ancestor and
 *     @p intersector do not actually intersect, or if no candidate rotation
 *     fixed it.
 */
[[nodiscard]] TreeNode* handle_intersection_with_ancestor(TreeNode& ancestor, TreeNode& intersector,
                                                          const PuzzlerOptions& opts,
                                                          ResolverState& state);

/**
 * Check @p node against every one of its ancestors (root-to-leaves order,
 * nearest ancestor first... actually nearest-to-furthest, matching the
 * vendored `while` loop), and, if none of those intersect, against the
 * exterior baseline. Ported from `checkNodeAgainstAncestors`
 * (`handleAncestorIntersections.inc:947`); this is `resolve.cpp`'s ancestor
 * branch entry point.
 *
 * @return The node whose `Config` changed to fix an ancestor intersection,
 *     or `nullptr` if @p node has no ancestor (or exterior) intersection,
 *     or none could be fixed.
 */
[[nodiscard]] TreeNode* check_node_against_ancestors(TreeNode& node, const PuzzlerOptions& opts,
                                                      ResolverState& state);

/**
 * Rebuild @p exterior's synthetic `LoopBox`/`StemBox`/`Aabb` so it can be
 * tested (via `intersect_node_node`) as if it were an ordinary box-bearing
 * node, spanning from @p top_level_ancestor's loop center out to
 * @p intersector's actual extent. Ported from the LIVE variant,
 * `TENTATIVE3_setupExteriorBoundingBoxes`
 * (`handleAncestorIntersections.inc:477`) -- see this file's scope note.
 *
 * DEVIATION (established at Milestone A step 4, not new here): the native
 * `LoopBox`/`StemBox` (`types.hpp`) have no `parent` field -- the vendored
 * `loop->parent`/`stem->parent` writes this function's vendored counterpart
 * makes are write-only dead fields (verified: `grep` finds no read of
 * `lBox->parent`/`sBox->parent` anywhere in the vendored sources), so they
 * are not ported.
 */
void setup_exterior_bounding_boxes(TreeNode& exterior, const TreeNode& top_level_ancestor,
                                   const TreeNode& intersector, const PuzzlerOptions& opts);

/**
 * The rotation angle needed to fix the intersection between @p ancestor and
 * @p intersector, rotating about @p rotation_node, for intersection type
 * @p it. Ported from `getRotationAngle` (`rotationAngle.inc:715`), which
 * dispatches to one of nine `getRotationAngleXxY` solvers
 * (`rotation_angle2.cpp`).
 *
 * @param rotation_sign `+1`/`-1` (see `get_rotation_sign`); `0` is never
 *     passed (callers only reach here after `get_rotation_sign != 0`).
 * @param clearance `PuzzlerOptions::clearance`; threaded down to
 *     `fix_intersection_of_rectangle_and_circle`/`fix_intersection_of_circles`,
 *     which mirror the vendored `epsilonFix` macro's clearance-scaling
 *     (`definitions.inc:60`) -- read there from a process-global in the
 *     reference, an explicit parameter throughout this port.
 */
[[nodiscard]] double get_rotation_angle(const TreeNode& ancestor, const TreeNode& rotation_node,
                                        const TreeNode& intersector, IntersectionType it,
                                        short rotation_sign, double clearance);

/**
 * The rotation angle to resolve a circle-vs-oriented-rectangle intersection:
 * rotate the circle at @p mobile_circ_center (radius @p mobile_circ_radius)
 * about @p rotation_center by the returned angle so it clears the rectangle
 * centered at @p static_rect_center (axes @p static_rect_vec_a/@p
 * static_rect_vec_b, half-extents @p static_rect_length_a/@p
 * static_rect_length_b). Ported from `fixIntersectionOfRectangleAndCircle`
 * (`rotationAngle.inc:64`); shared by `rotation_angle2.cpp`'s LxS/SxL/SxB/BxS
 * solvers.
 */
[[nodiscard]] double fix_intersection_of_rectangle_and_circle(
    Vec2 static_rect_center, Vec2 static_rect_vec_a, Vec2 static_rect_vec_b,
    double static_rect_length_a, double static_rect_length_b, Vec2 mobile_circ_center,
    double mobile_circ_radius, Vec2 rotation_center, short rotation_sign, double clearance);

/**
 * The rotation angle to resolve a circle-vs-circle intersection: rotate the
 * circle at @p mobile_circle_center (radius @p mobile_circle_radius) about
 * @p rotation_center by the returned angle so it clears the circle at
 * @p static_circle_center (radius @p static_circle_radius). Ported from
 * `fixIntersectionOfCircles` (`rotationAngle.inc:196`); shared by
 * `rotation_angle2.cpp`'s LxL/LxB/BxL/BxB solvers.
 *
 * HARDENED (not a vendored behavior change, a UB guard -- see
 * `rotation_angle.cpp`): the vendored function only early-returns `0.0` when
 * `getCutPointsOfCircles` finds EXACTLY zero cut points; when the two
 * circles COINCIDE (`get_cut_points_of_circles`'s `count == -1`), the
 * vendored code reads its cut-point buffers UNINITIALIZED. This port treats
 * `count <= 0` as the early-return case (the `count == -1` branch is
 * otherwise unreachable in practice -- it requires two RNA-structure loop
 * circles within 1.0 unit of both center and radius).
 */
[[nodiscard]] double fix_intersection_of_circles(Vec2 static_circle_center,
                                                 double static_circle_radius,
                                                 Vec2 mobile_circle_center,
                                                 double mobile_circle_radius, Vec2 rotation_center,
                                                 short rotation_sign, double clearance);

}  // namespace rna_layout
