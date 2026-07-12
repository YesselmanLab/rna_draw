/**
 * @file optimize.cpp
 * @brief Standalone geometry/config helpers for the optimize pass, ported
 *        from `optimize.inc` -- see `resolve_internal.hpp`'s "Optimization
 *        pass" section for the module split and its dead-code SCOPE NOTE.
 *        `optimize_node`/`optimize_tree_recursive`/`optimize_tree` (the
 *        driver-facing entry point) are `optimize2.cpp`.
 */

#include <cmath>

#include "resolve_internal.hpp"
#include "rna_layout/bounding_wedge.hpp"
#include "rna_layout/config.hpp"
#include "rna_layout/config_tree.hpp"
#include "rna_layout/geometry.hpp"
#include "rna_layout/intersect_tree.hpp"

namespace rna_layout {

namespace {

/// Whether applying @p new_radius (with @p deltas, possibly empty) would
/// actually change anything about @p cfg. Ported from `somethingChanged`
/// (`optimize.inc:308`).
bool something_changed(const Config& cfg, double old_radius, double new_radius,
                       const std::vector<double>& deltas) {
  bool changed = (new_radius - old_radius != 0.0);
  if (!changed && !deltas.empty()) {
    for (std::size_t current_arc = 0; current_arc < cfg.arcs.size(); ++current_arc) {
      if (deltas[current_arc] != 0.0) {
        changed = true;
        break;
      }
    }
  }
  return changed;
}

/// Apply @p deltas to @p node's `Config` at @p target_radius, but only if
/// `something_changed` says the result would differ from the current
/// state. Ported from `applyDeltas` (`optimize.inc:340`).
void apply_deltas(TreeNode& node, std::vector<double>& deltas, double target_radius,
                  const PuzzlerOptions& opts) {
  const Config& cfg = *node.cfg;  // NOLINT(bugprone-unchecked-optional-access)
  if (something_changed(cfg, cfg.radius, target_radius, deltas)) {
    apply_changes_to_config_and_bounding_boxes(node, deltas, target_radius, opts.paired,
                                               opts.unpaired);
  }
}

}  // namespace

bool check_optimize_intersections(const std::vector<const TreeNode*>& subtree,
                                  const std::vector<const TreeNode*>& ancestor_list,
                                  const PuzzlerOptions& opts) {
  return intersect_node_lists(subtree, subtree, opts.check_exterior, opts.clearance) ||
         intersect_node_lists(subtree, ancestor_list, opts.check_exterior, opts.clearance);
}

double shrink_loop_radius(TreeNode& node, const std::vector<const TreeNode*>& subtree,
                          const std::vector<const TreeNode*>& ancestor_list,
                          const PuzzlerOptions& opts) {
  Config& cfg = *node.cfg;  // NOLINT(bugprone-unchecked-optional-access)
  const double max_radius = cfg.radius;
  const double min_radius = cfg.min_radius;
  const double min_valid_radius = max_radius;

  constexpr double kMinAbsoluteDelta = 1.0;
  if (max_radius - min_radius < kMinAbsoluteDelta) {
    // No change -> shrinkingRatio == 1.0.
    return 1.0;
  }

  double radius = min_radius;
  const double delta = 0.1 * (max_radius - min_radius);

  int current_step = 0;
  constexpr int kMaxSteps = 10;

  while (current_step < kMaxSteps) {
    apply_changes_to_config_and_bounding_boxes(node, {}, radius, opts.paired, opts.unpaired);

    const bool intersecting = check_optimize_intersections(subtree, ancestor_list, opts);
    if (intersecting) {
      radius += delta;
    } else {
      break;
    }
    ++current_step;
  }

  if (current_step >= kMaxSteps || cfg.radius > max_radius) {
    apply_changes_to_config_and_bounding_boxes(node, {}, min_valid_radius, opts.paired,
                                               opts.unpaired);
  }

  return cfg.radius / max_radius;
}

void get_spaces(const TreeNode& node, int config_size, double paired_angle, double clearance,
                std::vector<double>& space) {
  std::vector<double> bounds_left(static_cast<std::size_t>(config_size));
  std::vector<double> bounds_right(static_cast<std::size_t>(config_size));

  bounds_left[0] = 0.0 + 0.5 * paired_angle;
  for (int i = 0; i < config_size - 1; ++i) {
    const AngleRange wedge = bounding_wedge(node, i, clearance);
    bounds_right[static_cast<std::size_t>(i)] = wedge.min_angle;
    bounds_left[static_cast<std::size_t>(i + 1)] = wedge.max_angle;
  }
  bounds_right[static_cast<std::size_t>(config_size - 1)] = geom::kTwoPi - 0.5 * paired_angle;

  space.assign(static_cast<std::size_t>(config_size), 0.0);
  for (int i = 0; i < config_size; ++i) {
    space[static_cast<std::size_t>(i)] =
        bounds_right[static_cast<std::size_t>(i)] - bounds_left[static_cast<std::size_t>(i)];
  }
}

void apply_config(TreeNode& node, const Config& target_config, const PuzzlerOptions& opts) {
  const Config& cfg = *node.cfg;  // NOLINT(bugprone-unchecked-optional-access)
  std::vector<double> deltas(cfg.arcs.size());
  for (std::size_t current_arc = 0; current_arc < cfg.arcs.size(); ++current_arc) {
    deltas[current_arc] = target_config.arcs[current_arc].angle - cfg.arcs[current_arc].angle;
  }
  apply_deltas(node, deltas, target_config.radius, opts);
}

void compute_alphas(std::vector<double>& alphas, const Config& cfg, double paired_distance) {
  // VENDORED QUIRK preserved for fidelity -- see this function's
  // declaration (`resolve_internal.hpp`).
  const auto truncated_paired_distance = static_cast<double>(static_cast<int>(paired_distance));
  const double paired_angle = geom::distance_to_angle(cfg.radius, truncated_paired_distance);

  alphas.assign(cfg.arcs.size(), 0.0);
  for (std::size_t current_arc = 0; current_arc < cfg.arcs.size(); ++current_arc) {
    alphas[current_arc] =
        (cfg.arcs[current_arc].angle - paired_angle) / cfg.arcs[current_arc].segments;
  }
}

void compute_increases(std::vector<int>& increase, int decrease_index, int config_size) {
  increase.clear();
  for (int i = 0; i < config_size; ++i) {
    if (i != decrease_index) {
      increase.push_back(i);
    }
  }
}

void compute_deltas(std::vector<double>& deltas, int decrease_index, double decrease_angle,
                    const std::vector<ConfigArc>& cfg_arcs, const std::vector<double>& alphas,
                    const std::vector<int>& increase) {
  double sum_increase_alphas = 0.0;
  for (int index : increase) {
    sum_increase_alphas += cfg_arcs[static_cast<std::size_t>(index)].segments *
                           alphas[static_cast<std::size_t>(index)];
  }

  for (int index : increase) {
    deltas[static_cast<std::size_t>(index)] =
        (cfg_arcs[static_cast<std::size_t>(index)].segments *
         alphas[static_cast<std::size_t>(index)] / sum_increase_alphas) *
        decrease_angle;
  }
  deltas[static_cast<std::size_t>(decrease_index)] = -decrease_angle;
}

bool search_best_config(TreeNode& node, std::vector<double>& deltas,
                        const std::vector<const TreeNode*>& subtree,
                        const std::vector<const TreeNode*>& ancestor_list,
                        const PuzzlerOptions& opts) {
  // A live reference (not a copy) into `node`'s own `Config`: every
  // `apply_deltas` call below mutates `node.cfg` in place, so `cfg.radius`
  // always reads back the FRESH radius each call actually applied --
  // exactly the vendored `cfg->radius`-re-read-at-each-call-site pattern
  // (`optimize.inc:626,644`), not a stale snapshot.
  const Config& cfg = *node.cfg;  // NOLINT(bugprone-unchecked-optional-access)
  apply_deltas(node, deltas, cfg.radius, opts);

  constexpr int kNumSteps = 10;
  constexpr double kFactor = 1.0 / kNumSteps;
  for (double& delta : deltas) {
    delta *= -kFactor;
  }

  bool intersecting = check_optimize_intersections(subtree, ancestor_list, opts);

  if (intersecting) {
    for (int current_step = 0; current_step < kNumSteps - 1; ++current_step) {
      apply_deltas(node, deltas, cfg.radius, opts);
      intersecting = check_optimize_intersections(subtree, ancestor_list, opts);
      if (!intersecting) {
        break;
      }
    }
  }

  return !intersecting;
}

bool can_shrink(const std::vector<double>& alphas, double unpaired_angle) {
  for (double alpha : alphas) {
    if (alpha <= unpaired_angle) {
      return false;
    }
  }
  return true;
}

}  // namespace rna_layout
