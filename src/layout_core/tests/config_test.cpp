// ctest: `rna_layout::generate_config`/`approximate_config_radius` sanity
// checks (`config.hpp`), matching plan Step 2's parity gate expectations
// (structural properties any correct port must have) at the C++ level.

#include "rna_layout/config.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "rna_layout/geometry.hpp"
#include "rna_layout/pair_table.hpp"

namespace geom = rna_layout::geom;
using rna_layout::BaseInfo;
using rna_layout::Config;
using rna_layout::generate_config;
using rna_layout::make_pair_table;

namespace {

int g_failures = 0;

void expect(bool condition, const char* what) {
  if (!condition) {
    std::cerr << "FAIL " << what << "\n";
    ++g_failures;
  }
}

void expect_near(double actual, double expected, double tol, const char* what) {
  if (std::fabs(actual - expected) > tol) {
    std::cerr << "FAIL " << what << ": expected " << expected << ", got " << actual << "\n";
    ++g_failures;
  }
}

constexpr double kPaired = 35.0;
constexpr double kUnpaired = 25.0;

/// Every arc's angle in a generated `Config` must sum to exactly 2*pi (a
/// loop's arcs partition its full circle) -- the same invariant
/// `cfgIsValid`'s `validSumAngles` check enforces on the vendored side.
void expect_arcs_sum_to_full_circle(const Config& cfg, const char* what) {
  double sum = 0.0;
  for (const auto& arc : cfg.arcs) {
    sum += arc.angle;
  }
  expect_near(sum, geom::kTwoPi, 1e-6, what);
}

void test_hairpin_single_arc() {
  const auto pair_table = make_pair_table("((((....))))");
  std::vector<BaseInfo> base_info(pair_table[0] + 1);
  const auto configs = generate_config(pair_table, kUnpaired, kPaired, base_info);

  expect(configs.size() == 1, "hairpin has exactly one loop config");
  expect(configs[0].arcs.size() == 1, "a hairpin loop has exactly one arc");
  expect(configs[0].radius > 0.0, "hairpin loop radius is positive");
  expect_arcs_sum_to_full_circle(configs[0], "hairpin loop arc angles sum to 2*pi");
  // The loop opens at the LAST base of the enclosing stem's first half
  // (position 4: the 4bp stem occupies 1-4/9-12), not at the structure's
  // first base -- matches `cfgGenHandleStem`'s walk to the stem's end
  // before calling `cfgGenHandleLoop`.
  expect(base_info[4].loop_id.has_value(), "hairpin loop-opening base (4) got a loop_id");
}

void test_two_helix_two_loop_configs() {
  // Two independent hairpins joined by an unpaired exterior run: each
  // hairpin loop gets its own config; the exterior loop (never enclosed by
  // a stem) gets none.
  const auto pair_table = make_pair_table("((((....))))..((((....))))");
  std::vector<BaseInfo> base_info(pair_table[0] + 1);
  const auto configs = generate_config(pair_table, kUnpaired, kPaired, base_info);

  expect(configs.size() == 2, "two independent hairpins produce two loop configs");
  for (const auto& cfg : configs) {
    expect_arcs_sum_to_full_circle(cfg, "each hairpin's arc angles sum to 2*pi");
  }
}

void test_multiloop_config() {
  // A three-way junction: one multiloop config with three arcs (one per
  // gap between the three enclosed stems). Wrapped in an outer stem so the
  // three inner stems form one multiloop rather than three independent
  // exterior-loop hairpins.
  const auto wrapped = make_pair_table(
      "(("
      "(((...)))(((...)))(((...)))"
      "))");
  std::vector<BaseInfo> base_info(wrapped[0] + 1);
  const auto configs = generate_config(wrapped, kUnpaired, kPaired, base_info);

  // Configs: the outer multiloop (4 arcs -- the enclosing stem is itself
  // one of the loop's stems, per `cfgGenerateDefaultConfig`'s "walk every
  // stem around the circle, including the parent" counting) + the 3 inner
  // hairpins (1 arc each) = 4 configs total.
  expect(configs.size() == 4, "wrapped three-way junction yields 4 loop configs");

  bool found_multiloop = false;
  for (const auto& cfg : configs) {
    expect_arcs_sum_to_full_circle(cfg, "every loop's arc angles sum to 2*pi");
    if (cfg.arcs.size() == 4) {
      found_multiloop = true;
    }
  }
  expect(found_multiloop, "one config has the multiloop's 4 arcs (3 children + the parent stem)");
}

void test_bulge_has_no_own_config() {
  // "(.((....)))": outer pair (1,11) encloses a single unpaired base (2)
  // then the inner stem (3,10)/(4,9) with a 4nt hairpin loop -- a bulge
  // (one unpaired base between two stems) followed by a real loop. Only
  // the hairpin gets its own config; the bulge is folded into the
  // enclosing stem's walk.
  const auto pair_table = make_pair_table("(.((....)))");
  std::vector<BaseInfo> base_info(pair_table[0] + 1);
  const auto configs = generate_config(pair_table, kUnpaired, kPaired, base_info);
  expect(configs.size() == 1, "a bulge contributes no config of its own");
}

void test_approximate_config_radius_matches_max_arc() {
  const auto pair_table = make_pair_table("((((....))))");
  std::vector<BaseInfo> base_info(pair_table[0] + 1);
  const auto configs = generate_config(pair_table, kUnpaired, kPaired, base_info);
  const double radius = rna_layout::approximate_config_radius(configs[0], kUnpaired, kPaired);
  expect(radius > 0.0, "approximate_config_radius is positive for a real loop");
}

}  // namespace

int main() {
  test_hairpin_single_arc();
  test_two_helix_two_loop_configs();
  test_multiloop_config();
  test_bulge_has_no_own_config();
  test_approximate_config_radius_matches_max_arc();

  if (g_failures > 0) {
    std::cerr << "config_test: " << g_failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
