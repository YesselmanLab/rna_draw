/**
 * @file puzzler.cpp
 * @brief Implements `puzzler.hpp`.
 */

#include "rna_layout/puzzler.hpp"

#include <memory>
#include <stdexcept>

#include "rna_layout/config_tree.hpp"
#include "rna_layout/nucleotide_coords.hpp"
#include "rna_layout/pair_table.hpp"
#include "rna_layout/resolve_exterior.hpp"
#include "rna_layout/tree.hpp"

namespace rna_layout {

namespace {

/// Throws if @p opts requests the not-yet-ported resolver. See
/// `puzzler.hpp`'s file header.
void require_resolver_off(const PuzzlerOptions& opts) {
  if (opts.check_sibling || opts.check_ancestor || opts.optimize) {
    throw std::logic_error(
        "rna_layout::layout_puzzler: the intersection resolver "
        "(check_sibling/check_ancestor/optimize) is not implemented yet "
        "(Milestone A steps 7-9); pass a PuzzlerOptions with all three "
        "false for the resolver-off finalization path.");
  }
}

/**
 * Build + (if requested) canonicalize the config tree. Ported from
 * `vrna_plot_coords_puzzler_pt`'s setup through its `checkAndFixIntersections`
 * call (`RNApuzzler.c:441-478`) -- MINUS that call itself.
 *
 * WHY SKIPPING `checkAndFixIntersections` HERE IS EXACT, NOT AN
 * APPROXIMATION: `require_resolver_off` (called by every public entry
 * point before this runs) guarantees `check_sibling == check_ancestor ==
 * optimize == false`. Tracing `checkAndFixIntersections`'s own driver
 * (`resolveIntersections.inc:24`) under that exact condition: its
 * `while (checkTree)` body only calls `checkNodeAgainstAncestors` when
 * `puzzler->checkAncestorIntersections` is set (false here) and only calls
 * `checkSiblings` when `puzzler->checkSiblingIntersections` is set (false
 * here) -- so the loop body never mutates anything and `checkTree` is never
 * set back to `1`, meaning the `while` runs its guarded checks once (both
 * skipped) then falls through; the "OPTIMIZATIONS" block is gated by
 * `puzzler->optimize` (false here) and is skipped too. The function returns
 * `NULL` having touched nothing. So, for this exact option combination,
 * calling `checkAndFixIntersections` is PROVABLY equivalent to not calling
 * it at all -- the tree `update_bounding_boxes` leaves behind is already
 * the function's fixed point.
 */
std::unique_ptr<TreeNode> run_config_tree_pipeline(const std::vector<int>& pair_table,
                                                   const TurtleLayout& turtle,
                                                   const PuzzlerOptions& opts) {
  const double bulge_dist = stem_bulge_distance(opts.unpaired);
  std::unique_ptr<TreeNode> tree =
      build_config_tree(pair_table, turtle.base_info, turtle.configs, turtle.coords, bulge_dist);

  // Mirrors `RNApuzzler.c:472-478`'s outer gate: `updateBoundingBoxes` runs
  // whenever ANY of the three checks is requested (even though, per
  // `require_resolver_off`, `check_sibling`/`check_ancestor` are always
  // false here -- `check_exterior` alone can still trigger this, matching
  // the vendored default of `checkExteriorIntersections = 1`).
  if (opts.check_exterior || opts.check_sibling || opts.check_ancestor) {
    update_bounding_boxes(*tree, opts.paired, opts.unpaired);
  }
  return tree;
}

void validate_pair_table_nonempty(const std::vector<int>& pair_table) {
  if (pair_table.empty() || pair_table[0] <= 0) {
    throw std::invalid_argument("layout_puzzler requires a non-empty pair table");
  }
}

void validate_nonempty(const std::string& structure) {
  if (structure.empty()) {
    throw std::invalid_argument("layout_puzzler requires a non-empty structure");
  }
}

void validate_no_empty_loop(const std::string& structure) {
  if (structure.find("()") != std::string::npos) {
    throw std::invalid_argument(
        "structure contains an empty loop \"()\": not supported by layout_puzzler");
  }
}

}  // namespace

Coords layout_puzzler(const std::vector<int>& pair_table, const PuzzlerOptions& opts) {
  require_resolver_off(opts);
  validate_pair_table_nonempty(pair_table);

  const int length = pair_table[0];
  const TurtleLayout turtle = run_turtle_layout(pair_table, opts.paired, opts.unpaired);
  const std::unique_ptr<TreeNode> tree = run_config_tree_pipeline(pair_table, turtle, opts);

  Coords coords{std::vector<double>(static_cast<std::size_t>(length)),
                std::vector<double>(static_cast<std::size_t>(length))};
  determine_nucleotide_coordinates(*tree, pair_table, opts.unpaired, opts.paired, coords);

  // Runs unconditionally (matching `RNApuzzler.c:498-507`'s hardcoded
  // `checkIntersectionsOfExteriorBranches = 1`) -- see `resolve_exterior.hpp`.
  resolve_exterior_children_intersection(*tree, pair_table, opts.unpaired, opts.allow_flipping,
                                         opts.clearance, coords);

  return coords;
}

Coords layout_puzzler(const std::string& dot_bracket, const PuzzlerOptions& opts) {
  validate_nonempty(dot_bracket);
  validate_no_empty_loop(dot_bracket);
  const std::vector<int> pair_table = make_pair_table(dot_bracket);  // validates well-nestedness
  return layout_puzzler(pair_table, opts);
}

}  // namespace rna_layout
