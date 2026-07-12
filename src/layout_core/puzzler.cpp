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
#include "rna_layout/resolve.hpp"
#include "rna_layout/resolve_exterior.hpp"
#include "rna_layout/tree.hpp"

namespace rna_layout {

namespace {

/// The vendored `<= 0 -> 25000` config-change-budget fallback
/// (`RNApuzzler.c:460-461`), applied once per resolve pass.
constexpr int kDefaultMaxConfigChanges = 25000;

/// Throws if @p opts requests a not-yet-ported resolver stage. See
/// `puzzler.hpp`'s file header: `check_sibling` (Milestone A step 7) and
/// `check_ancestor` (step 8) are both real; `optimize` (step 9) is not.
void require_optimize_off(const PuzzlerOptions& opts) {
  if (opts.optimize) {
    throw std::logic_error(
        "rna_layout::layout_puzzler: the optimize pass is not implemented "
        "yet (Milestone A step 9); pass a PuzzlerOptions with optimize = "
        "false (check_sibling/check_ancestor may each be true or false).");
  }
}

/**
 * Build the config tree and, if any of `check_exterior`/`check_sibling`/
 * `check_ancestor` is requested, canonicalize its boxes and run the
 * SIBLING/ANCESTOR resolver. Ported from `vrna_plot_coords_puzzler_pt`'s
 * setup through its `checkAndFixIntersections` call (`RNApuzzler.c:441-478`).
 *
 * `require_optimize_off` (called by every public entry point before this
 * runs) guarantees `optimize == false`; if BOTH `check_sibling` and
 * `check_ancestor` are also false, `check_and_fix_intersections` is a no-op
 * (see `resolve.cpp`'s driver: with every gate false, its `while` loop runs
 * its guarded checks once, all skipped, and falls through untouched).
 */
std::unique_ptr<TreeNode> run_config_tree_pipeline(const std::vector<int>& pair_table,
                                                   const TurtleLayout& turtle,
                                                   const PuzzlerOptions& opts) {
  const double bulge_dist = stem_bulge_distance(opts.unpaired);
  std::unique_ptr<TreeNode> tree =
      build_config_tree(pair_table, turtle.base_info, turtle.configs, turtle.coords, bulge_dist);

  // Mirrors `RNApuzzler.c:472-478`'s outer gate: both `updateBoundingBoxes`
  // and `checkAndFixIntersections` run whenever ANY of the three checks is
  // requested (matching the vendored default of `checkExteriorIntersections
  // = 1`, even when `check_sibling`/`check_ancestor` are both off).
  if (opts.check_exterior || opts.check_sibling || opts.check_ancestor) {
    update_bounding_boxes(*tree, opts.paired, opts.unpaired);

    ResolverState state;
    state.max_config_changes =
        opts.max_config_changes <= 0 ? kDefaultMaxConfigChanges : opts.max_config_changes;
    check_and_fix_intersections(tree.get(), opts, state);
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
  require_optimize_off(opts);
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
