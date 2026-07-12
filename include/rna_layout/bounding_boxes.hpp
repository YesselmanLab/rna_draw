#pragma once

/**
 * @file bounding_boxes.hpp
 * @brief `LoopBox`/`StemBox` construction and bulge geometry, ported from
 *        `boundingBoxes.inc`.
 *
 * These are the pure box-BUILDING primitives `config_tree.hpp`'s
 * `build_config_tree`/`update_bounding_boxes` compose into a whole tree;
 * see that header for the two-stage build ("initial boxes from the
 * turtle's raw x/y" then "canonical boxes from the config-tree geometry")
 * `RNApuzzler.c` itself performs.
 *
 * Preserves the reference's arithmetic and operation order throughout (per
 * `.claude/plans/current-plan.md`'s style notes and Task 1's fidelity rule):
 * this is a deterministic pass, parity-gated against the vendored oracle at
 * a tight numeric tolerance.
 */

#include <vector>

#include "rna_layout/turtle.hpp"
#include "rna_layout/types.hpp"

namespace rna_layout {

/// The three points `getBulgeCoordinatesExtraDistance`/`getBulgeCoordinates`
/// (`boundingBoxes.inc:126/177`) return: the bulge base's position just
/// before, at, and just after its peak, in world coordinates.
struct BulgePoints {
  Vec2 prev;
  Vec2 at;
  Vec2 next;
};

/**
 * Build the `LoopBox` for the loop starting at @p start: its center (offset
 * from the closing base pair by the loop's fitted radius, in the direction
 * away from the pair, oriented by the loop's clockwise/counter-clockwise
 * sense) and radius. Ported from `buildLoopBox`/`getLoopData`
 * (`boundingBoxes.inc:331/265`).
 *
 * @param start The loop-opening base (1-indexed, matching `pair_table`).
 * @param pair_table 1-indexed pair table.
 * @param base_info Per-base state; `base_info[start].loop_id` must be set
 *     (i.e. @p start must not be a bulge-folded loop-opening base).
 * @param configs The loop `Config`s `generate_config` produced.
 * @param coords The turtle-base pass's Cartesian coordinates (0-indexed).
 */
[[nodiscard]] LoopBox build_loop_box(int start, const std::vector<int>& pair_table,
                                     const std::vector<BaseInfo>& base_info,
                                     const std::vector<Config>& configs, const Coords& coords);

/**
 * Build the `StemBox` for the stem spanning @p start (its first base) to
 * @p end (the loop it opens), including every bulge notch along it. Ported
 * from `buildStemBox`/`createStemBox`/`countBulges`/`setBulges`
 * (`boundingBoxes.inc:511/350/389/468`).
 *
 * @param start The stem's first base (`TreeNode::stem_start`).
 * @param end The stem's enclosed loop's opening base (`TreeNode::loop_start`).
 * @param pair_table 1-indexed pair table.
 * @param coords The turtle-base pass's Cartesian coordinates (0-indexed).
 * @param bulge_dist Extra clearance distance for bulge points (vendored
 *     `distBulge`, `RNApuzzler.c:442`); see `config_tree.hpp`'s
 *     `stem_bulge_distance`.
 */
[[nodiscard]] StemBox build_stem_box(int start, int end, const std::vector<int>& pair_table,
                                     const Coords& coords, double bulge_dist);

/**
 * The axis-aligned bounding box enclosing @p stem_box's rectangle corners,
 * @p loop_box's circle extrema, and every bulge peak. Ported from
 * `updateAABB` (`configtree.inc:255`).
 */
[[nodiscard]] Aabb compute_aabb(const StemBox& stem_box, const LoopBox& loop_box);

/**
 * The bulge at @p index's before/at/after positions, offset an extra
 * @p extra_distance further out along the stem's `b` axis. Ported from
 * `getBulgeCoordinatesExtraDistance` (`boundingBoxes.inc:126`).
 */
[[nodiscard]] BulgePoints bulge_coordinates_extra_distance(const StemBox& stem, int index,
                                                           double extra_distance);

/**
 * `bulge_coordinates_extra_distance(stem, index, 0.0)`. Ported from
 * `getBulgeCoordinates` (`boundingBoxes.inc:177`).
 */
[[nodiscard]] BulgePoints bulge_coordinates(const StemBox& stem, int index);

/// Translate @p box's center by @p vector. Ported from `translateLoopBox`
/// (`boundingBoxes.inc:187`).
void translate_loop_box(LoopBox& box, Vec2 vector);

/// Translate @p box's center by @p vector. Ported from `translateStemBox`
/// (`boundingBoxes.inc:219`).
void translate_stem_box(StemBox& box, Vec2 vector);

/**
 * A `LoopBox` with the given @p center/@p radius, no other bookkeeping.
 * Ported from `createLoopBox` (`boundingBoxes.inc:60`); used by the
 * ancestor-resolution exterior-box setup (`exterior_boxes.cpp`, Milestone A
 * step 8) to build a SYNTHETIC loop box directly from computed geometry,
 * unlike `build_loop_box`'s pair-table-driven construction.
 */
[[nodiscard]] LoopBox create_loop_box(Vec2 center, double radius);

/**
 * A `StemBox` from three raw corner points (`s`: south/start corner,
 * `e`: north/end corner, `sp`: south-prime, the opposite long edge's start
 * corner). Ported from `createStemBox` (`boundingBoxes.inc:350`); the same
 * primitive `build_stem_box` uses internally, exposed here for
 * `exterior_boxes.cpp` (Milestone A step 8), which builds synthetic stem
 * boxes directly from computed geometry rather than pair-table positions.
 */
[[nodiscard]] StemBox create_stem_box(Vec2 s, Vec2 e, Vec2 sp);

}  // namespace rna_layout
