/**
 * @file config.cpp
 * @brief Implements per-loop config generation (`config.hpp`), ported from
 *        `drawingconfig.inc`.
 */

#include "rna_layout/config.hpp"

#include <cmath>

#include "rna_layout/geometry.hpp"

namespace rna_layout {

namespace {

/// Newton-iteration cap, matching the vendored `MAX_ITERATIONS`
/// (`drawingconfig.inc:276`).
constexpr int kMaxRadiusIterations = 1000;

/**
 * One Newton step of `approximate_config_arc_radius`'s iteration: the
 * closed-form derivative of `stems * asin(paired / (2r)) + backbones *
 * asin(unpaired / (2r)) - angle / 2` with respect to `r`, evaluated at
 * `radius`. Extracted from the loop body purely to keep the iteration
 * itself readable; the arithmetic and evaluation order are preserved
 * exactly from `drawingconfig.inc:314` for parity.
 */
double newton_step(double paired, double unpaired, int stems, int backbones, double angle,
                   double radius) {
  const double numerator = 2 * (stems * std::asin(paired / (2 * radius)) +
                                backbones * std::asin(unpaired / (2 * radius)) - (angle / 2));
  const double denominator =
      -(paired * stems / (radius * std::sqrt(radius * radius - paired * paired / 4)) +
        unpaired * backbones / (radius * std::sqrt(radius * radius - unpaired * unpaired / 4)));
  return numerator / denominator;
}

}  // namespace

double approximate_config_arc_radius(double paired, double unpaired, int stems, int backbones,
                                     double angle) {
  const double lower_bound = (unpaired / 2) / std::sin((angle / (stems + backbones)) / 2);
  const double upper_bound = (paired / 2) / std::sin((angle / (stems + backbones)) / 2);

  double radius = 0.5 * (lower_bound + upper_bound);
  radius = std::fmax(radius, 0.5 * paired);
  radius = std::fmax(radius, 0.5 * unpaired);

  for (int iteration = 0; iteration < kMaxRadiusIterations; ++iteration) {
    const double dx = newton_step(paired, unpaired, stems, backbones, angle, radius);
    radius -= dx;
    if (std::fabs(dx) < geom::kEpsilon3) {
      break;
    }
  }

  if (radius < lower_bound) {
    radius = lower_bound;
  } else if (radius > upper_bound) {
    radius = upper_bound;
  }
  return radius;
}

double approximate_config_radius(const Config& cfg, double unpaired, double paired) {
  double radius = 0;
  for (const ConfigArc& arc : cfg.arcs) {
    const double candidate =
        approximate_config_arc_radius(paired, unpaired, /*stems=*/1, arc.segments, arc.angle);
    if (candidate > radius) {
      radius = candidate;
    }
  }
  return radius;
}

namespace {

/// Count how many stems open directly off the loop starting at
/// `pair_table[start]`'s partner -- i.e. how many `ConfigArc`s the loop
/// needs -- without allocating anything yet. First pass of
/// `cfgGenerateDefaultConfig` (`drawingconfig.inc:400`).
int count_loop_arcs(const std::vector<int>& pair_table, int start) {
  int num_arcs = 0;
  int i = start + 1;
  while (i <= pair_table[start]) {
    if (pair_table[i] == 0) {
      ++i;
    } else {
      ++num_arcs;
      if (i != pair_table[start]) {
        i = pair_table[i] + 1;
      } else {
        break;
      }
    }
  }
  return num_arcs;
}

/// Fill in each `ConfigArc`'s angle/segment-count for the loop starting at
/// `start`, given the paired/unpaired per-base angles already computed.
/// Second pass of `cfgGenerateDefaultConfig` (`drawingconfig.inc:421`).
std::vector<ConfigArc> fill_loop_arcs(const std::vector<int>& pair_table, int start,
                                      double angle_paired, double angle_unpaired, int num_arcs) {
  std::vector<ConfigArc> arcs(num_arcs);
  int arc_unpaired = 0;
  int current_arc = 0;
  int i = start + 1;
  while (i <= pair_table[start]) {
    if (pair_table[i] == 0) {
      ++arc_unpaired;
      ++i;
    } else {
      arcs[current_arc] =
          ConfigArc{arc_unpaired + 1, angle_paired + (arc_unpaired + 1) * angle_unpaired};
      ++current_arc;
      if (i != pair_table[start]) {
        arc_unpaired = 0;
        i = pair_table[i] + 1;
      } else {
        break;
      }
    }
  }
  return arcs;
}

/// Build the default `Config` for the loop starting at `start`, given a
/// pre-fitted `radius`. Ported from `cfgGenerateDefaultConfig`
/// (`drawingconfig.inc:379`).
Config generate_default_config(const std::vector<int>& pair_table, int start, double unpaired,
                               double paired, double radius) {
  Config cfg;
  cfg.radius = radius;
  cfg.min_radius = radius;
  cfg.default_radius = radius;

  const double angle_paired = 2 * std::asin(paired / (2 * radius));
  const double angle_unpaired = 2 * std::asin(unpaired / (2 * radius));
  const int num_arcs = count_loop_arcs(pair_table, start);
  cfg.arcs = fill_loop_arcs(pair_table, start, angle_paired, angle_unpaired, num_arcs);
  return cfg;
}

void gen_handle_stem(int base_nr, const std::vector<int>& pair_table,
                     std::vector<BaseInfo>& base_info, std::vector<Config>& configs,
                     double unpaired, double paired);

/// Count the stems and unpaired bases directly enclosed by the loop
/// starting at `start`, used to detect the bulge special case. Ported from
/// the counting loop at the top of `cfgGenHandleLoop`
/// (`drawingconfig.inc:487`).
void count_loop_members(const std::vector<int>& pair_table, int start, int end, int& unpaired_count,
                        int& stem_count) {
  unpaired_count = 0;
  stem_count = 1;
  int i = start + 1;
  while (i < end) {
    if (pair_table[i] == 0) {
      ++unpaired_count;
      ++i;
    } else if (pair_table[i] > i) {
      ++stem_count;
      i = pair_table[i];
    } else {
      ++i;
    }
  }
}

/// Recurse into every child stem of the loop starting at `start`, giving
/// each its own `Config` via `gen_handle_stem`. Second half of the
/// non-bulge branch of `cfgGenHandleLoop` (`drawingconfig.inc:524`).
void gen_handle_loop_children(int start, int end, const std::vector<int>& pair_table,
                              std::vector<BaseInfo>& base_info, std::vector<Config>& configs,
                              double unpaired, double paired) {
  int i = start + 1;
  while (i < end) {
    // The "unpaired base" and "returned from a stem" arms both just
    // advance `i` -- identical bodies for two conceptually different
    // reasons, matching the reference's own three-way branch
    // (`drawingconfig.inc:524`) exactly; not a copy-paste mistake.
    if (pair_table[i] == 0) {  // NOLINT(bugprone-branch-clone)
      ++i;
    } else if (pair_table[i] > i) {
      gen_handle_stem(i, pair_table, base_info, configs, unpaired, paired);
      i = pair_table[i];
    } else {
      ++i;
    }
  }
}

/// Ported from `cfgGenHandleLoop` (`drawingconfig.inc:475`): recursively
/// generate `Config`s for the loop starting at `base_nr` and every loop
/// nested inside it, skipping the single-unpaired-base bulge case (which
/// has no `Config` of its own -- see `config.hpp`).
void gen_handle_loop(int base_nr, const std::vector<int>& pair_table,
                     std::vector<BaseInfo>& base_info, std::vector<Config>& configs,
                     double unpaired, double paired) {
  const int start = base_nr;
  const int end = pair_table[base_nr];

  int unpaired_count = 0;
  int stem_count = 0;
  count_loop_members(pair_table, start, end, unpaired_count, stem_count);

  // This bulge test (`stem_count == 2 && unpaired_count == 1`) MUST stay in
  // lockstep with `turtle.cpp`'s independent bulge detector
  // (`detect_bulge` + `n - m == 1`, `handle_loop`): if they ever disagree,
  // a loop gets skipped here (no `Config`, no `loop_id`) but is walked as
  // the non-bulge case there, and `configs[base_info[start].loop_id.value()]`
  // (`turtle.cpp:277`) throws on the unset `optional`.
  const bool is_bulge = (stem_count == 2 && unpaired_count == 1);
  if (is_bulge) {
    const int stem_start = (pair_table[start + 1] == 0) ? start + 2 : start + 1;
    gen_handle_stem(stem_start, pair_table, base_info, configs, unpaired, paired);
    return;
  }

  const double default_radius = approximate_config_arc_radius(
      paired, unpaired, stem_count, unpaired_count + stem_count, geom::kTwoPi);
  configs.push_back(generate_default_config(pair_table, start, unpaired, paired, default_radius));
  base_info[start].loop_id = static_cast<int>(configs.size()) - 1;

  gen_handle_loop_children(start, end, pair_table, base_info, configs, unpaired, paired);
}

/// Ported from `cfgGenHandleStem` (`drawingconfig.inc:553`): walk the
/// stem starting at `base_nr` to its first loop (a stem is not a bulge
/// candidate itself -- `pair_table[i+1] == pair_table[i] - 1` just says
/// "the next base pair is nested one step in, i.e. still a helix").
void gen_handle_stem(int base_nr, const std::vector<int>& pair_table,
                     std::vector<BaseInfo>& base_info, std::vector<Config>& configs,
                     double unpaired, double paired) {
  int i = base_nr;
  while (pair_table[i + 1] == pair_table[i] - 1) {
    ++i;
  }
  gen_handle_loop(i, pair_table, base_info, configs, unpaired, paired);
}

}  // namespace

bool cfg_is_valid(const Config& cfg, const std::vector<double>& delta_cfg) {
  double sum_angles = 0.0;
  bool valid_single_angles = true;

  for (std::size_t current_arc = 0; current_arc < cfg.arcs.size(); ++current_arc) {
    const double angle = cfg.arcs[current_arc].angle + delta_cfg[current_arc];
    sum_angles += angle;

    const bool valid_angle = 0.0 < angle && angle < geom::kTwoPi;
    valid_single_angles = valid_single_angles && valid_angle;
  }

  const bool valid_sum_angles = std::fabs(sum_angles - geom::kTwoPi) < geom::kEpsilon3;
  return valid_single_angles && valid_sum_angles;
}

namespace {

/// Ported from `cfgUpdateMinRadius` (`drawingconfig.inc:636`): exactly
/// `approximate_config_radius`, given a new name here since
/// `cfg_apply_changes` reads it as "the smallest radius this config's
/// current arcs can be drawn at", not as a general-purpose radius fit.
double cfg_min_radius(const Config& cfg, double unpaired, double paired) {
  return approximate_config_radius(cfg, unpaired, paired);
}

}  // namespace

double cfg_apply_changes(Config& cfg, const std::vector<double>& delta_cfg, double radius_new,
                         double unpaired, double paired) {
  if (!delta_cfg.empty()) {
    for (std::size_t current_arc = 0; current_arc < cfg.arcs.size(); ++current_arc) {
      cfg.arcs[current_arc].angle += delta_cfg[current_arc];
    }
  }

  const double old_radius = cfg.radius;
  double new_radius = -1.0;
  if (radius_new > 0.0) {
    // The minimum of valid and input as the new radius.
    cfg.min_radius = cfg_min_radius(cfg, unpaired, paired);
    new_radius = std::fmax(radius_new, cfg.min_radius);
    cfg.radius = new_radius;
  } else if (radius_new == 0.0) {
    // The minRadius as new value (allows shrinking a loop).
    cfg.min_radius = cfg_min_radius(cfg, unpaired, paired);
    new_radius = cfg.min_radius;
    cfg.radius = new_radius;
  } else if (radius_new == -1.0) {
    // The minRadius as new value (forbids shrinking a loop).
    cfg.min_radius = cfg_min_radius(cfg, unpaired, paired);
    if (cfg.min_radius - geom::kEpsilon0 > old_radius) {
      new_radius = cfg.min_radius;
    } else {
      constexpr double kDefaultIncrease = 1.05;
      new_radius = old_radius * kDefaultIncrease;
    }
    cfg.radius = new_radius;
  }
  // All other inputs are unhandled (vendored default: `newRadius = -1.0`,
  // returned WITHOUT touching `cfg.radius`) -- no caller in this port passes
  // any other value, so that branch is omitted rather than modeled.

  return new_radius;
}

std::vector<Config> generate_config(const std::vector<int>& pair_table, double unpaired,
                                    double paired, std::vector<BaseInfo>& base_info) {
  std::vector<Config> configs;
  const int length = pair_table[0];
  int i = 1;
  while (i < length) {
    // See `gen_handle_loop_children`'s comment above: the "unpaired" and
    // "returned from stem" arms intentionally share a body, matching
    // `cfgGenerateConfig` (`drawingconfig.inc:578`).
    if (pair_table[i] == 0) {  // NOLINT(bugprone-branch-clone)
      ++i;
    } else if (pair_table[i] > i) {
      gen_handle_stem(i, pair_table, base_info, configs, unpaired, paired);
      i = pair_table[i];
    } else {
      ++i;
    }
  }
  return configs;
}

}  // namespace rna_layout
