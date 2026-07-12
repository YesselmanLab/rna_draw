#pragma once

/**
 * @file puzzler.hpp
 * @brief The RNApuzzler orchestration entry point, ported from
 *        `vrna_plot_coords_puzzler_pt` (`RNApuzzler.c:391`) -- Milestone A
 *        step 6 ("finalization with resolver OFF").
 *
 * THIS SLICE ONLY WIRES THE RESOLVER-OFF PATH: `layout_puzzler` requires
 * `PuzzlerOptions::check_sibling == false`, `check_ancestor == false`, and
 * `optimize == false` (it throws otherwise) -- the intersection RESOLVER
 * itself (`checkNodeAgainstAncestors`, `checkSiblings`, `optimizeTree`) is
 * Milestone A steps 7-9, not yet ported. With those three options false,
 * the vendored driver `checkAndFixIntersections` is PROVABLY a no-op (see
 * `puzzler.cpp`'s doc comment on `run_config_tree_pipeline`), so this port
 * skips implementing/calling it rather than porting dead code early --
 * everything downstream of it (`determine_nucleotide_coordinates`,
 * `resolve_exterior_children_intersection`) is real, always-on
 * finalization, ported in full this step.
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
 * @param opts Layout options; `check_sibling`, `check_ancestor`, and
 *     `optimize` MUST all be `false` (see this file's header).
 * @return `Coords` of size `pair_table[0]`.
 * @throws std::logic_error If any of `check_sibling`/`check_ancestor`/
 *     `optimize` is `true`.
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
