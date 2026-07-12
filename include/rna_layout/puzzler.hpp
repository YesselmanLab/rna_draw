#pragma once

/**
 * @file puzzler.hpp
 * @brief The RNApuzzler orchestration entry point, ported from
 *        `vrna_plot_coords_puzzler_pt` (`RNApuzzler.c:391`).
 *
 * SCOPE (Milestone A step 9, COMPLETE): `layout_puzzler` runs all three
 * resolver stages -- SIBLING, ANCESTOR, and OPTIMIZE
 * (`check_and_fix_intersections`, `resolve.hpp`) -- gated independently by
 * `PuzzlerOptions::check_sibling`/`check_ancestor`/`optimize`.
 * `PuzzlerOptions{}`'s own defaults (all three `true`, plus
 * `check_exterior = true`, `allow_flipping = false`,
 * `max_config_changes = 25000`, `clearance = 1.0`) reproduce the vendored
 * `vrna_plot_options_puzzler()`'s defaults field-for-field -- i.e. calling
 * `layout_puzzler` with a default-constructed `PuzzlerOptions` is the
 * native counterpart of the vendored PRODUCTION pipeline
 * (`_vienna_layout.plot_coords_puzzler`/`plot_coords_puzzler_opts`, the
 * path `rna_draw.layout.production` actually ships).
 * `determine_nucleotide_coordinates` and
 * `resolve_exterior_children_intersection` (Milestone A step 6) are
 * unaffected -- always-on finalization downstream of whatever the resolver
 * left behind.
 */

#include <string>
#include <vector>

#include "rna_layout/turtle.hpp"
#include "rna_layout/types.hpp"

namespace rna_layout {

/**
 * Lay out @p pair_table with the RNApuzzler pipeline: `generate_config` ->
 * turtle base layout -> `build_config_tree` -> (if any of
 * `check_exterior`/`check_sibling`/`check_ancestor`/`optimize` is set)
 * `update_bounding_boxes` -> `check_and_fix_intersections` ->
 * `determine_nucleotide_coordinates` ->
 * `resolve_exterior_children_intersection`.
 *
 * @param pair_table 1-indexed pair table (`make_pair_table`'s output).
 * @param opts Layout options; every field may be set independently (see
 *     this file's header for the default-`PuzzlerOptions` == vendored-
 *     production-config correspondence).
 * @return `Coords` of size `pair_table[0]`.
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
 */
[[nodiscard]] Coords layout_puzzler(const std::string& dot_bracket, const PuzzlerOptions& opts);

}  // namespace rna_layout
