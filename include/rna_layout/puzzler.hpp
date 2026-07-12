#pragma once

/**
 * @file puzzler.hpp
 * @brief The RNApuzzler orchestration entry point, ported from
 *        `vrna_plot_coords_puzzler_pt` (`RNApuzzler.c:391`).
 *
 * SCOPE (Milestone A step 8): `layout_puzzler` now runs both the SIBLING and
 * ANCESTOR intersection resolvers (`check_and_fix_intersections`,
 * `resolve.hpp`) when `PuzzlerOptions::check_sibling`/`check_ancestor` are
 * `true`. `optimize` MUST still be `false` (it throws otherwise) --
 * `optimizeTree` is Milestone A step 9, not yet ported.
 * `determine_nucleotide_coordinates` and
 * `resolve_exterior_children_intersection` (Milestone A step 6) are
 * unaffected -- always-on finalization downstream of whatever the resolver
 * (or its absence) left behind.
 */

#include <string>
#include <vector>

#include "rna_layout/turtle.hpp"
#include "rna_layout/types.hpp"

namespace rna_layout {

/**
 * Lay out @p pair_table with the RNApuzzler pipeline, resolver disabled:
 * `generate_config` -> turtle base layout -> `build_config_tree` ->
 * (if any of `check_exterior`/`check_sibling`/`check_ancestor` is set)
 * `update_bounding_boxes` -> `determine_nucleotide_coordinates` ->
 * `resolve_exterior_children_intersection`.
 *
 * @param pair_table 1-indexed pair table (`make_pair_table`'s output).
 * @param opts Layout options; `optimize` MUST be `false` (see this file's
 *     header); `check_sibling`/`check_ancestor` may each be `true`.
 * @return `Coords` of size `pair_table[0]`.
 * @throws std::logic_error If `optimize` is `true`.
 * @throws std::invalid_argument If @p pair_table is empty (length `<= 0`).
 */
[[nodiscard]] Coords layout_puzzler(const std::vector<int>& pair_table, const PuzzlerOptions& opts);

/**
 * Convenience overload: validate @p dot_bracket the same way
 * `layout_turtle(const std::string&)` does (non-empty, well-nested, no bare
 * `"()"` empty loop), build its pair table, and delegate to the pair-table
 * overload.
 *
 * @throws std::invalid_argument If @p dot_bracket is malformed, or (see the
 *     pair-table overload) if @p pair_table ends up empty.
 * @throws std::logic_error See the pair-table overload.
 */
[[nodiscard]] Coords layout_puzzler(const std::string& dot_bracket, const PuzzlerOptions& opts);

}  // namespace rna_layout
