#pragma once

/**
 * @file config_tree.hpp
 * @brief Config-tree construction + navigation, ported from `configtree.inc`.
 *
 * `build_config_tree` replaces `buildConfigtree` + its `treeHandle{Stem,Loop}`
 * recursion (`RAII` via `TreeNode::add_child` stands in for the vendored
 * `vrna_alloc`'d `children` array + `freeTree`), and `update_bounding_boxes`
 * replaces `updateBoundingBoxes`: `RNApuzzler.c`'s entry point runs BOTH, in
 * that order, before the (not-yet-ported) resolver ever sees the tree --
 * see `RNApuzzler.c:448-477`.
 *
 * SCOPE (Milestone A step 4): only construction + the navigation getters the
 * config-tree/bounding-box/wedge code itself needs (`is_exterior`,
 * `get_loop_center`, `get_stem_center`, `get_child_angle`) are ported this
 * slice. `configtree.inc`'s resolver-only helpers
 * (`applyChangesToConfigAndBoundingBoxes`, `getChildIndex`/`getChildNode`,
 * `countSubtreeNodes`/`countAncestorNodes`, `isInteriorLoop`/`isMultiLoop`,
 * `getPairedAngle`) are deferred to the resolver steps (Milestone A step 7+)
 * that actually call them.
 *
 * `translate_bounding_boxes` (Milestone A step 6) is the one exception:
 * `resolveExteriorChildrenIntersectionXY` (`resolve_exterior.cpp`) needs it
 * to shift a whole exterior child's subtree, and that pass runs regardless
 * of the resolver's `check_*`/`optimize` options (`RNApuzzler.c:498-507`
 * calls it unconditionally) -- see `resolve_exterior.hpp`'s file header.
 */

#include <memory>
#include <vector>

#include "rna_layout/tree.hpp"
#include "rna_layout/turtle.hpp"
#include "rna_layout/types.hpp"

namespace rna_layout {

/**
 * The extra clearance distance bulge points get in stem-box construction
 * (vendored `distBulge`, `RNApuzzler.c:442-446`): `sqrt(unpaired^2 -
 * 0.25*unpaired^2)`.
 */
[[nodiscard]] double stem_bulge_distance(double unpaired);

/**
 * Build the config tree for @p pair_table: one node per loop (plus the
 * exterior loop as the root), each carrying its `Config` (from @p configs)
 * and an initial `LoopBox`/`StemBox` built from the turtle-base pass's
 * Cartesian coordinates (@p coords). Ported from `buildConfigtree`
 * (`configtree.inc:867`) plus its `treeHandleStem`/`treeHandleLoop`/
 * `buildBoundingBoxes` helpers.
 *
 * The initial boxes this attaches are NOT the canonical, config-geometry-
 * derived boxes the resolver operates on -- call `update_bounding_boxes` on
 * the result next, matching `RNApuzzler.c`'s own two-stage sequence.
 *
 * @param pair_table 1-indexed pair table.
 * @param base_info Per-base state from the turtle pass (`run_turtle_layout`);
 *     only `loop_id` is read.
 * @param configs The loop `Config`s `generate_config` produced.
 * @param coords The turtle-base pass's Cartesian coordinates.
 * @param bulge_dist Extra bulge clearance (`stem_bulge_distance`).
 * @return The tree's root (the exterior loop; `parent == nullptr`,
 *     `cfg == std::nullopt`, `lbox`/`sbox == std::nullopt`).
 */
[[nodiscard]] std::unique_ptr<TreeNode> build_config_tree(const std::vector<int>& pair_table,
                                                          const std::vector<BaseInfo>& base_info,
                                                          const std::vector<Config>& configs,
                                                          const Coords& coords, double bulge_dist);

/**
 * Recompute every node's canonical `StemBox`/`LoopBox`/`Aabb` from the
 * config-tree geometry (loop radii, per-arc angles) top-down, replacing the
 * turtle-coordinate-derived initial boxes `build_config_tree` attached.
 * Ported from `updateBoundingBoxes` (`configtree.inc:328`).
 *
 * @param node The tree (or subtree) root to update, in place; recurses into
 *     every descendant.
 * @param paired Distance between the two bases of a base pair.
 * @param unpaired Default backbone-step distance.
 */
void update_bounding_boxes(TreeNode& node, double paired, double unpaired);

/// Whether @p node is the tree's root (the exterior loop). Mirrors
/// `isExterior` (`configtree.inc:561`); the vendored check is `id == 0`,
/// equivalent here to "has no parent" since only the root has no parent.
[[nodiscard]] bool is_exterior(const TreeNode& node);

/// @p node's loop center, i.e. its `LoopBox::center`. Mirrors
/// `getLoopCenter` (`configtree.inc:1021`). `node.lbox` must be set (not
/// the root).
[[nodiscard]] Vec2 get_loop_center(const TreeNode& node);

/// @p node's stem center, i.e. its `StemBox::c`. Mirrors `getStemCenter`
/// (`configtree.inc:1037`). `node.sbox` must be set (not the root).
[[nodiscard]] Vec2 get_stem_center(const TreeNode& node);

/**
 * The clockwise angle, at @p parent's loop center, from @p parent's stem
 * center to @p child's loop center (possibly outside `[0, 2*pi)`). Ported
 * from `getChildAngle` (`configtree.inc:965`); `child` must be a direct
 * child of `parent` with both boxes already built (`parent` must not be the
 * root -- it has no `sbox`).
 */
[[nodiscard]] double get_child_angle(const TreeNode& parent, const TreeNode& child);

/**
 * Translate @p node's `StemBox`/`LoopBox`/`Aabb` by @p vector, and recurse
 * into every descendant (translating its boxes by the SAME vector). Ported
 * from `translateBoundingBoxes` (`configtree.inc:906`).
 *
 * @param node The (sub)tree root to translate, in place; must not be the
 *     true exterior root (it has no boxes to translate).
 * @param vector The translation to apply.
 */
void translate_bounding_boxes(TreeNode& node, Vec2 vector);

}  // namespace rna_layout
