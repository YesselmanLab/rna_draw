/**
 * @file deltas.cpp
 * @brief `calc_deltas` + its three private helpers, ported from
 *        `calcDeltas.inc`.
 *
 * Pure computation: given a loop and two of its children (indices
 * `index_left`/`index_right`) that need more angular space between them,
 * decides how to redistribute the loop's `Config::arcs` angles to make
 * room, WITHOUT mutating anything -- `check_and_apply_config_changes`
 * (`config_changes.cpp`) is the only thing that applies the result.
 *
 * SIBLING-PATH DEAD CODE (documented, not simplified away -- fidelity
 * rule): `calc_deltas`'s "can this bend go to a higher tree level instead"
 * search is UNREACHABLE for every call this port makes (Milestone A step 7,
 * sibling-only): `fix_intersection_of_siblings` (`resolve_siblings.cpp`)
 * always passes `recursive_end = node.parent`, and the search's own first
 * step reads `TreeNode* parent = node.parent` -- so `parent == recursive_end`
 * is true before the loop's first condition check ever runs, and the `while`
 * body never executes (`can_go_higher` stays `false`). This is ported
 * verbatim anyway (not simplified to an `if (true)`) because Milestone A
 * step 8's ancestor path calls `calc_deltas` with a DIFFERENT,
 * more-distant `recursive_end` where this search is real and load-bearing.
 */

#include <algorithm>
#include <cmath>

#include "resolve_internal.hpp"
#include "rna_layout/bounding_wedge.hpp"
#include "rna_layout/config.hpp"
#include "rna_layout/geometry.hpp"

namespace rna_layout {

namespace {

/// Distribute @p target_angle equally across every arc @p increase marks,
/// adding to @p delta_cfg in place. Ported from
/// `calcDeltasEquidistantIncrease` (`calcDeltas.inc:47`).
void calc_deltas_equidistant_increase(double target_angle, int config_size,
                                      const std::vector<bool>& increase,
                                      std::vector<double>& delta_cfg) {
  int increase_count = 0;
  for (int i = 0; i < config_size; ++i) {
    if (increase[static_cast<std::size_t>(i)]) {
      ++increase_count;
    }
  }
  const double delta_per_increase = target_angle / increase_count;

  for (int i = 0; i < config_size; ++i) {
    if (increase[static_cast<std::size_t>(i)]) {
      delta_cfg[static_cast<std::size_t>(i)] += delta_per_increase;
    }
  }
}

/// Greedily shrink the single largest available "space" (an arc's current
/// angle minus `2 * min_angle_half`) between @p index_left and
/// @p index_right, repeating until @p target_angle_in is exhausted or no
/// space remains. Ported from `calcDeltasMaximumFirstDecrease`
/// (`calcDeltas.inc:73`).
double calc_deltas_maximum_first_decrease(double target_angle_in, int index_left, int index_right,
                                          int config_size, std::vector<double>& delta_cfg,
                                          const std::vector<double>& current_angles,
                                          double min_angle_half) {
  double target_angle = target_angle_in;

  bool do_loop = true;
  while (do_loop) {
    double max_space = 0.0;
    int max_space_index = -1;
    int i = 0;

    if (index_left == -1) {
      double sum_angles = 0.0;
      i = -1;
      while (i != index_right) {
        ++i;
        const double cfg = current_angles[static_cast<std::size_t>(i)] +
                           delta_cfg[static_cast<std::size_t>(i)] - 2 * min_angle_half;
        sum_angles += cfg;
      }
      while (i != config_size - 1) {
        ++i;
        const double cfg = current_angles[static_cast<std::size_t>(i)] +
                           delta_cfg[static_cast<std::size_t>(i)] - 2 * min_angle_half;
        if (sum_angles < geom::kPi) {
          if (cfg > max_space) {
            max_space = cfg;
            max_space_index = i;
          }
        } else {
          break;
        }
        sum_angles += cfg;
      }
    } else if (index_right == -1) {
      double sum_angles = 0.0;
      i = config_size - 1;
      while (i != index_left) {
        const double cfg = current_angles[static_cast<std::size_t>(i)] +
                           delta_cfg[static_cast<std::size_t>(i)] - 2 * min_angle_half;
        sum_angles += cfg;
        --i;
      }
      while (i != -1) {
        const double cfg = current_angles[static_cast<std::size_t>(i)] +
                           delta_cfg[static_cast<std::size_t>(i)] - 2 * min_angle_half;
        if (sum_angles < geom::kPi) {
          if (cfg > max_space) {
            max_space = cfg;
            max_space_index = i;
          }
        } else {
          break;
        }
        sum_angles += cfg;
        --i;
      }
    } else {
      i = index_right;
      if (i == config_size - 1) {
        i = -1;
      }
      while (i != index_left) {
        const double cfg = current_angles[static_cast<std::size_t>(i) + 1] +
                           delta_cfg[static_cast<std::size_t>(i) + 1] - 2 * min_angle_half;
        if (cfg > max_space) {
          max_space = cfg;
          max_space_index = i + 1;
        }
        ++i;
        if (i == config_size - 1) {
          i = -1;
        }
      }
    }

    double diff = 0.0;
    if (max_space_index != -1) {
      const double factor = (target_angle < 0.1 * target_angle_in) ? 1.0 : 0.5;
      diff = -std::fmin(factor * max_space, target_angle);
      delta_cfg[static_cast<std::size_t>(max_space_index)] += diff;
      target_angle += diff;
    }

    do_loop = target_angle > 0.0 && std::fabs(diff) > geom::kEpsilon3;
  }

  return target_angle;
}

/// Alternately shrink the arcs nearest @p index_left/@p index_right
/// (walking outward from both ends toward the midpoint between them),
/// splitting @p target_angle_in evenly across every arc @p decrease marks
/// in each pass. Ported from `calcDeltasNearestNeighborsFirstDecrease`
/// (`calcDeltas.inc:183`).
double calc_deltas_nearest_neighbors_first_decrease(double target_angle_in, int index_left,
                                                    int index_right, int config_size,
                                                    const std::vector<bool>& decrease,
                                                    const std::vector<double>& space,
                                                    std::vector<double>& delta_cfg) {
  double target_angle = target_angle_in;

  int steps = 0;
  int stem_it = index_right;
  while (stem_it != index_left) {
    ++stem_it;
    if (stem_it == config_size) {
      stem_it = -1;
    }
    ++steps;
  }
  const int num_it = steps / 2;  // implicit floor()

  std::vector<int> index(static_cast<std::size_t>(steps));
  bool changed = true;
  while (changed) {
    changed = false;
    int count = 0;

    int i_l = index_left;
    if (i_l == -1) {
      i_l = config_size - 1;
    }
    int i_r = index_right + 1;
    if (i_r == config_size) {
      i_r = 0;
    }

    for (int i = 0; i < num_it; ++i) {
      if (decrease[static_cast<std::size_t>(i_l)]) {
        index[static_cast<std::size_t>(count)] = i_l;
        ++count;
      }
      if (decrease[static_cast<std::size_t>(i_r)]) {
        index[static_cast<std::size_t>(count)] = i_r;
        ++count;
      }

      --i_l;
      if (i_l == -1) {
        i_l = config_size - 1;
      }
      ++i_r;
      if (i_r == config_size) {
        i_r = 0;
      }
    }

    if (num_it < 0.5 * steps) {
      index[static_cast<std::size_t>(count)] = i_l;
      ++count;
      --i_l;
      if (i_l == -1) {
        // `i_l`'s new value is never read again this call -- ported
        // verbatim from `calcDeltasNearestNeighborsFirstDecrease`
        // (`calcDeltas.inc:262-263`), which has the identical dead store;
        // preserved for expression-tree fidelity, not a bug.
        i_l = config_size - 1;  // NOLINT(clang-analyzer-deadcode.DeadStores)
      }
    }

    if (count > 0) {
      const double part_angle = target_angle / count;
      for (int k = 0; k < count; ++k) {
        const int j = index[static_cast<std::size_t>(k)];
        if (decrease[static_cast<std::size_t>(j)]) {
          const double diff = -std::fmin(
              space[static_cast<std::size_t>(j)] + delta_cfg[static_cast<std::size_t>(j)],
              part_angle);
          delta_cfg[static_cast<std::size_t>(j)] += diff;
          target_angle += diff;
          changed = changed || (diff != 0.0);
        }
      }
    }
  }

  return target_angle;
}

}  // namespace

double calc_deltas(const TreeNode& node, const TreeNode* recursive_end, int index_left,
                   int index_right, double delta_angle, double paired, double clearance,
                   std::vector<double>& deltas) {
  if (delta_angle < 0.0) {
    return 0.0;
  }

  const int child_count = static_cast<int>(node.children.size());
  const int config_size = child_count + 1;

  const Config& cfg =
      *node.cfg;  // NOLINT(bugprone-unchecked-optional-access) -- see resolve_siblings.cpp
  const double min_outer_angle = std::asin(paired / (2 * cfg.radius));

  std::vector<double> angles_min(static_cast<std::size_t>(child_count));
  std::vector<double> angles_max(static_cast<std::size_t>(child_count));
  std::vector<double> space(static_cast<std::size_t>(config_size));
  std::vector<double> delta_cfg(static_cast<std::size_t>(config_size), 0.0);
  std::vector<bool> increase(static_cast<std::size_t>(config_size));
  std::vector<bool> decrease(static_cast<std::size_t>(config_size));
  std::vector<double> current_angles(static_cast<std::size_t>(config_size));

  for (std::size_t current_arc = 0; current_arc < cfg.arcs.size(); ++current_arc) {
    current_angles[current_arc] = cfg.arcs[current_arc].angle;  // getArcAngle(cfg, currentArc)
  }

  for (int current_child = 0; current_child < child_count; ++current_child) {
    const AngleRange wedge = bounding_wedge(node, current_child, clearance);
    angles_min[static_cast<std::size_t>(current_child)] = wedge.min_angle;
    angles_max[static_cast<std::size_t>(current_child)] = wedge.max_angle;
  }

  // Convert bounding wedges to "free" areas usable for compensation.
  space[0] = angles_min[0] - (0 + min_outer_angle);
  for (int i = 1; i < config_size - 1; ++i) {
    space[static_cast<std::size_t>(i)] =
        angles_min[static_cast<std::size_t>(i)] - angles_max[static_cast<std::size_t>(i - 1)];
  }
  space[static_cast<std::size_t>(config_size - 1)] =
      (geom::kTwoPi - min_outer_angle) - angles_max[static_cast<std::size_t>(config_size - 2)];

  // Fix too-big spaces (may exceed the config for very large loops).
  for (int i = 0; i < config_size; ++i) {
    space[static_cast<std::size_t>(i)] =
        std::fmin(space[static_cast<std::size_t>(i)],
                  cfg.arcs[static_cast<std::size_t>(i)].angle - 2 * min_outer_angle);
  }

  // Mark increase and decrease areas (every index is visited exactly once
  // by the two `while` loops below, so no entry is ever read "unset").
  int current_index = index_left;
  while (current_index != index_right) {
    increase[static_cast<std::size_t>(current_index) + 1] = true;
    decrease[static_cast<std::size_t>(current_index) + 1] = false;
    ++current_index;
    if (current_index == config_size - 1) {
      current_index = -1;
    }
  }
  while (current_index != index_left) {
    increase[static_cast<std::size_t>(current_index) + 1] = false;
    decrease[static_cast<std::size_t>(current_index) + 1] =
        space[static_cast<std::size_t>(current_index) + 1] > 0.0;
    ++current_index;
    if (current_index == config_size - 1) {
      current_index = -1;
    }
  }

  double target_angle = delta_angle;

  // Step 1: equidistant increase.
  calc_deltas_equidistant_increase(target_angle, config_size, increase, delta_cfg);

  // Step 2: nearest-neighbors-first decrease.
  target_angle = calc_deltas_nearest_neighbors_first_decrease(
      target_angle, index_left, index_right, config_size, decrease, space, delta_cfg);

  // Step 3: check if intersections are fixed.
  const bool not_fixed_yet = (target_angle != 0.0);
  if (not_fixed_yet) {
    // See this file's header: `can_go_higher` is provably always `false`
    // for the sibling-only call path this step exercises.
    const TreeNode* parent = node.parent;
    bool can_go_higher = false;

    while (parent != recursive_end && !is_exterior(*parent)) {
      if (is_multi_loop(*parent)) {
        can_go_higher = true;
        break;
      }
      // `parent` is a real loop, not the exterior root (guarded by
      // `!is_exterior(*parent)` above), so it always has a `cfg` -- see
      // `config_tree.cpp`'s invariant note. getArcAngle(parent->cfg, 0).
      const double child_angle =
          parent->cfg->arcs[0].angle;  // NOLINT(bugprone-unchecked-optional-access)
      if (std::fabs(child_angle - geom::kPi) < geom::kEpsilon3) {
        // no-op
      } else if (child_angle > geom::kPi) {
        if (index_left == 0) {
          can_go_higher = true;
          break;
        }
      } else if (child_angle < geom::kPi) {
        if (index_left == -1) {
          can_go_higher = true;
          break;
        }
      }
      parent = parent->parent;
    }

    if (!can_go_higher) {
      target_angle =
          calc_deltas_maximum_first_decrease(target_angle, index_left, index_right, config_size,
                                             delta_cfg, current_angles, min_outer_angle);
    }
  }

  // Step 4: equidistant increase with the negative remaining target angle.
  calc_deltas_equidistant_increase(-target_angle, config_size, increase, delta_cfg);

  // The vendored "fix too small changes" re-run here is disabled at this
  // call site (`fixTooSmallChanges = 0`, `calcDeltas.inc:471`) -- unlike
  // `check_and_apply_config_changes`'s own (enabled) copy of the same fix,
  // which is what actually runs; omitted here rather than ported dead.

  deltas.resize(static_cast<std::size_t>(config_size));
  for (int current_arc = 0; current_arc < config_size; ++current_arc) {
    deltas[static_cast<std::size_t>(current_arc)] =
        delta_cfg[static_cast<std::size_t>(current_arc)];
  }

  // Check that all deltas sum to zero.
  double check_sum = 0.0;
  for (int current_arc = 0; current_arc < config_size; ++current_arc) {
    check_sum += deltas[static_cast<std::size_t>(current_arc)];
  }
  if (std::fabs(check_sum) > geom::kEpsilon3) {
    std::fill(deltas.begin(), deltas.end(), 0.0);
    target_angle = delta_angle;
  }

  if (!cfg_is_valid(cfg, deltas)) {
    std::fill(deltas.begin(), deltas.end(), 0.0);
    target_angle = delta_angle;
  }

  return delta_angle - target_angle;
}

}  // namespace rna_layout
