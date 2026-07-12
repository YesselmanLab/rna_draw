/**
 * @file optimize2.cpp
 * @brief `optimize_node`/`optimize_tree_recursive`/`optimize_tree`, ported
 *        from `optimize.inc` -- see `resolve_internal.hpp`'s "Optimization
 *        pass" section for the module split. `optimize.cpp` holds this
 *        file's standalone helper dependencies.
 */

#include <algorithm>
#include <cmath>

#include "resolve_internal.hpp"
#include "rna_layout/config_tree.hpp"
#include "rna_layout/geometry.hpp"

namespace rna_layout {

namespace {

/// `optimize_node`'s per-arc search state, threaded through its three
/// extracted helpers below (a plain aggregate, not exposed outside this
/// file) -- splits the vendored `optimizeNode`'s single ~190-line function
/// (`optimize.inc:686`) along its own internal blank-line-delimited seams
/// without altering any expression, matching this port's fidelity rule.
struct OptimizeSearchState {
  int config_size = 0;
  int min_sorted_index = 0;
  double unpaired_angle = 0.0;
  std::vector<double> alphas;
  std::vector<double> spaces;
  std::vector<int> sorted;
  Config best_config;
};

/// Sort `indices` (initialized `0..values_level1.size()-1`) so
/// `values_level1[indices[i]]` is DECREASING, breaking near-ties
/// (`geom::kEpsilon7`) by `values_level2`. Ported from `bubblesort`
/// (`definitions.inc:101`); used only by `optimize_node` below (the
/// vendored function's one and only caller, `optimize.inc:780`).
void bubble_sort_descending(const std::vector<double>& values_level1,
                            const std::vector<double>& values_level2, std::vector<int>& indices) {
  const auto num_values = static_cast<int>(values_level1.size());
  indices.resize(static_cast<std::size_t>(num_values));
  for (int i = 0; i < num_values; ++i) {
    indices[static_cast<std::size_t>(i)] = i;
  }

  for (int i = 0; i < num_values - 1; ++i) {
    for (int j = 0; j < num_values - i - 1; ++j) {
      double this_value =
          values_level1[static_cast<std::size_t>(indices[static_cast<std::size_t>(j)])];
      double next_value =
          values_level1[static_cast<std::size_t>(indices[static_cast<std::size_t>(j + 1)])];
      bool swap = false;
      if (next_value - this_value > geom::kEpsilon7) {
        swap = true;
      } else if (std::fabs(next_value - this_value) < geom::kEpsilon7) {
        this_value = values_level2[static_cast<std::size_t>(indices[static_cast<std::size_t>(j)])];
        next_value =
            values_level2[static_cast<std::size_t>(indices[static_cast<std::size_t>(j + 1)])];
        if (next_value - this_value > geom::kEpsilon7) {
          swap = true;
        }
      }
      if (swap) {
        std::swap(indices[static_cast<std::size_t>(j)], indices[static_cast<std::size_t>(j + 1)]);
      }
    }
  }
}

/// The `if (configChanged) {...} else {...}` block of `optimizeNode`
/// (`optimize.inc:744-781`): refresh `s.alphas`/`s.unpaired_angle`, try
/// shrinking @p node's radius (keeping the shrunk config as the new
/// `s.best_config` iff it actually helped), and, only when
/// `s.min_sorted_index` is back at `0`, recompute `s.spaces`/`s.sorted`.
void refresh_alphas_and_maybe_shrink(TreeNode& node, OptimizeSearchState& s,
                                     const std::vector<const TreeNode*>& subtree,
                                     const std::vector<const TreeNode*>& ancestor_list,
                                     const PuzzlerOptions& opts) {
  Config& cfg = *node.cfg;  // NOLINT(bugprone-unchecked-optional-access)
  s.unpaired_angle = geom::distance_to_angle(cfg.radius, opts.unpaired);
  compute_alphas(s.alphas, cfg, opts.paired);

  if (can_shrink(s.alphas, s.unpaired_angle) &&
      shrink_loop_radius(node, subtree, ancestor_list, opts) < 1.0) {
    s.best_config = cfg;
    s.min_sorted_index = 0;
    s.unpaired_angle = geom::distance_to_angle(cfg.radius, opts.unpaired);
    compute_alphas(s.alphas, cfg, opts.paired);
  } else {
    apply_config(node, s.best_config, opts);
  }

  if (s.min_sorted_index == 0) {
    const double paired_angle = geom::distance_to_angle(cfg.radius, opts.paired);
    get_spaces(node, s.config_size, paired_angle, opts.clearance, s.spaces);
    bubble_sort_descending(s.alphas, s.spaces, s.sorted);
  }
}

/// The "Find first arc with sufficient space" loop of `optimizeNode`
/// (`optimize.inc:787-811`): the first arc (scanning `s.sorted` from
/// `s.min_sorted_index`) whose available space exceeds `2 * s.unpaired_
/// angle`, or `-1` if none does. Advances `s.min_sorted_index` past
/// whatever it scanned.
int find_decrease_index(OptimizeSearchState& s, double min_multiple) {
  for (int index = s.min_sorted_index; index < s.config_size; ++index) {
    const int current_arc = s.sorted[static_cast<std::size_t>(index)];
    double space = s.spaces[static_cast<std::size_t>(current_arc)];
    if (space > geom::kPi) {
      space = geom::kPi;
    }

    const double min_space = min_multiple * s.unpaired_angle;
    if (space > min_space) {
      s.min_sorted_index = index + 1;
      return current_arc;
    }
  }
  return -1;
}

/// The "Check if current arc has sufficient space" block of `optimizeNode`
/// (`optimize.inc:817-825`): how much to narrow @p decrease_index by.
double compute_decrease_angle(const Config& cfg, int decrease_index, const OptimizeSearchState& s) {
  const double space = s.spaces[static_cast<std::size_t>(decrease_index)];
  const int segments = cfg.arcs[static_cast<std::size_t>(decrease_index)].segments;
  const double min_necessary_space = segments * s.unpaired_angle;
  const double current_necessary_space =
      segments * s.alphas[static_cast<std::size_t>(decrease_index)];

  constexpr double kFactor = 0.5;
  return kFactor * std::fmin(current_necessary_space - min_necessary_space, space);
}

}  // namespace

double optimize_node(TreeNode& node, const std::vector<const TreeNode*>& subtree,
                     const std::vector<const TreeNode*>& ancestor_list, const PuzzlerOptions& opts,
                     ResolverState& state) {
  if (node.children.empty()) {
    // Nothing to do for hairpin loops.
    return 1.0;
  }

  Config& cfg = *node.cfg;  // NOLINT(bugprone-unchecked-optional-access)
  constexpr double kMinRadiusIncrease = 5.0;
  if (cfg.radius - cfg.default_radius < kMinRadiusIncrease) {
    // Nothing to do if the radius increase is small.
    return 1.0;
  }

  constexpr double kMinMultiple = 2.0;
  const auto config_size = static_cast<int>(cfg.arcs.size());

  const Config initial_config = cfg;
  const double initial_radius = cfg.radius;

  OptimizeSearchState s;
  s.config_size = config_size;
  s.best_config = initial_config;
  s.alphas.assign(static_cast<std::size_t>(config_size), 0.0);
  s.spaces.assign(static_cast<std::size_t>(config_size), 0.0);
  s.sorted.assign(static_cast<std::size_t>(config_size), 0);

  std::vector<double> deltas(static_cast<std::size_t>(config_size), 0.0);
  std::vector<int> increase;

  s.min_sorted_index = 0;
  int run_nr = 0;
  const int run_nr_max = 100 * config_size;  // just in case ...
  bool config_changed = true;

  while (s.min_sorted_index < config_size && run_nr < run_nr_max) {
    ++run_nr;

    if (config_changed) {
      refresh_alphas_and_maybe_shrink(node, s, subtree, ancestor_list, opts);
    } else {
      apply_config(node, s.best_config, opts);
    }

    const int decrease_index = find_decrease_index(s, kMinMultiple);
    if (decrease_index < 0) {
      // No arc found: leave.
      break;
    }

    const double decrease_angle = compute_decrease_angle(cfg, decrease_index, s);
    if (decrease_angle < geom::kEpsilon3) {
      // decreaseAngle too small: leave.
      continue;
    }

    compute_increases(increase, decrease_index, config_size);
    compute_deltas(deltas, decrease_index, decrease_angle, cfg.arcs, s.alphas, increase);

    config_changed = search_best_config(node, deltas, subtree, ancestor_list, opts);
  }

  // Apply the best configuration found so far.
  apply_config(node, s.best_config, opts);

  if (s.best_config.radius < initial_config.radius) {
    // Best radius smaller than the initial radius: log the change. (The
    // vendored `optimizeNode` also computes a per-arc `deltas` diff here
    // that it never reads again before freeing -- omitted; see
    // `optimize.inc:857-859`.)
    ++state.changes_applied;
  } else {
    // Otherwise revert to the initial state.
    apply_config(node, initial_config, opts);
  }

  return cfg.radius / initial_radius;
}

namespace {

/// @p node itself, then every descendant, DFS pre-order. Ported from
/// `countSubtreeNodes`/`collectSubtreeNodes` (`configtree.inc:475,505`),
/// combined into one pass since this port has no separate "count first,
/// allocate, then fill" step (`std::vector` grows itself).
void collect_subtree_nodes(const TreeNode& node, std::vector<const TreeNode*>& out) {
  out.push_back(&node);
  for (const auto& child : node.children) {
    collect_subtree_nodes(*child, out);
  }
}

/// @p node's ancestor chain, nearest first, up to and including the tree's
/// root. Ported from `countAncestorNodes`/`collectAncestorNodes`
/// (`configtree.inc:489,521`).
std::vector<const TreeNode*> collect_ancestor_nodes(const TreeNode& node) {
  std::vector<const TreeNode*> ancestors;
  const TreeNode* ancestor = node.parent;
  while (ancestor != nullptr) {
    ancestors.push_back(ancestor);
    ancestor = ancestor->parent;
  }
  return ancestors;
}

}  // namespace

double optimize_tree_recursive(TreeNode& node, const std::vector<const TreeNode*>& subtree,
                               const std::vector<const TreeNode*>& ancestor_list,
                               const PuzzlerOptions& opts, ResolverState& state) {
  double shrinking_ratio = 1.0;
  double ratio = 1.0;
  double min_ratio = 1.0;

  // Fidelity, not style (`optimize.inc:895-933`): the vendored
  // `optimizeTreeRecursive` is itself a `do { ... } while (minRatio < 1.0)`
  // -- it must run its body (recurse into children, then this node) at
  // least once before ever consulting `minRatio`, so a `while` loop
  // wouldn't be an equivalent rewrite here.
  do {  // NOLINT(cppcoreguidelines-avoid-do-while)
    if (state.changes_applied > state.max_config_changes) {
      // Also a vendored dead store (`optimize.inc:897`): `break` exits the
      // loop without ever consulting `min_ratio` again. Kept for fidelity
      // (same reasoning as the `do`/`while` shape above) rather than
      // dropped as "obviously" redundant.
      min_ratio = 1.0;  // NOLINT(clang-analyzer-deadcode.DeadStores)
      break;
    }

    min_ratio = 1.0;
    // Do-loop until nothing improves further: recursive optimization of
    // all children first.
    for (const auto& child : node.children) {
      ratio = optimize_tree_recursive(*child, subtree, ancestor_list, opts, state);
      min_ratio = std::fmin(ratio, min_ratio);
      shrinking_ratio *= ratio;
    }

    if (min_ratio < 1.0) {
      continue;
    }

    if (!is_exterior(node)) {
      // Shrink the current node.
      ratio = optimize_node(node, subtree, ancestor_list, opts, state);
      min_ratio = std::fmin(ratio, min_ratio);
      shrinking_ratio *= ratio;
    }
  } while (min_ratio < 1.0);

  return shrinking_ratio;
}

double optimize_tree(TreeNode& node, const PuzzlerOptions& opts, ResolverState& state) {
  if (!opts.optimize) {
    return 1.0;
  }

  double shrinking_ratio = 1.0;

  std::vector<const TreeNode*> subtree;
  collect_subtree_nodes(node, subtree);
  const std::vector<const TreeNode*> ancestor_list = collect_ancestor_nodes(node);

  if (!check_optimize_intersections(subtree, ancestor_list, opts)) {
    // Only start if the subtree does not intersect with the ancestor tree.
    shrinking_ratio = optimize_tree_recursive(node, subtree, ancestor_list, opts, state);
  }

  return shrinking_ratio;
}

}  // namespace rna_layout
