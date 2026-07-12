// ctest: `rna_layout::build_loop_box`/`build_stem_box`/`compute_aabb`/bulge
// geometry sanity checks (`bounding_boxes.hpp`). Parity against the
// vendored oracle is the Python-level gate
// (`tests/test_native_parity.py::TestTreeParity`); this file only proves
// the port's own internal shape/invariants and a couple of hand-computed
// analytic cases.

#include "rna_layout/bounding_boxes.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "rna_layout/config_tree.hpp"
#include "test_tree_helper.hpp"

using rna_layout::Aabb;
using rna_layout::LoopBox;
using rna_layout::StemBox;
using rna_layout::TreeNode;
using rna_layout::Vec2;

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

// NOLINTBEGIN(bugprone-unchecked-optional-access) -- see the matching note
// in `config_tree_test.cpp`: the `expect(...has_value()...)` calls just
// above each `->` access are the (clang-tidy-invisible) gate.
void test_initial_boxes_from_hairpin() {
  // `build_config_tree` itself calls `build_loop_box`/`build_stem_box` --
  // exercise them indirectly via a freshly-built (pre-`update_bounding_boxes`)
  // tree, then check the initial boxes are geometrically sane.
  const std::vector<int> pair_table = rna_layout::make_pair_table("((((....))))");
  const rna_layout::TurtleLayout turtle =
      rna_layout::run_turtle_layout(pair_table, kPaired, kUnpaired);
  const double bulge_dist = rna_layout::stem_bulge_distance(kUnpaired);
  const auto tree = rna_layout::build_config_tree(pair_table, turtle.base_info, turtle.configs,
                                                  turtle.coords, bulge_dist);

  const TreeNode& loop = *tree->children[0];
  expect(loop.lbox.has_value(), "initial lbox is set");
  expect(loop.sbox.has_value(), "initial sbox is set");
  expect(loop.lbox->radius > 0.0, "initial loop radius is positive");
  expect(loop.sbox->bulges.empty(), "hairpin stem has no bulges");
  expect_near(loop.sbox->bulge_dist, bulge_dist, 1e-12, "sbox.bulge_dist matches the input");
}
// NOLINTEND(bugprone-unchecked-optional-access)

void test_compute_aabb_bounds_stem_and_loop() {
  const StemBox stem{
      /*a=*/Vec2{1.0, 0.0},
      /*b=*/Vec2{0.0, 1.0},
      /*c=*/Vec2{0.0, 0.0},
      /*e=*/Vec2{2.0, 3.0},
      /*bulges=*/{},
      /*bulge_dist=*/0.0,
  };
  const LoopBox loop{Vec2{10.0, 0.0}, 4.0};

  const Aabb aabb = rna_layout::compute_aabb(stem, loop);

  // Stem corners span x in [-2, 2], y in [-3, 3]; loop spans x in [6, 14],
  // y in [-4, 4] -- the union's AABB must cover both.
  expect(aabb.min.x <= -2.0, "aabb.min.x covers the stem's left edge");
  expect(aabb.max.x >= 14.0, "aabb.max.x covers the loop's right edge");
  expect(aabb.min.y <= -4.0, "aabb.min.y covers the loop's bottom edge");
  expect(aabb.max.y >= 4.0, "aabb.max.y covers the loop's top edge");
}

void test_bulge_coordinates_extra_distance_offsets_along_b() {
  StemBox stem{
      /*a=*/Vec2{1.0, 0.0},
      /*b=*/Vec2{0.0, 1.0},
      /*c=*/Vec2{0.0, 0.0},
      /*e=*/Vec2{2.0, 1.0},
      /*bulges=*/{},
      /*bulge_dist=*/0.5,
  };
  stem.bulges.push_back(rna_layout::Bulge{/*sign=*/1.0, /*a_prev=*/-1.0, /*a_this=*/0.0,
                                          /*a_next=*/1.0});

  const rna_layout::BulgePoints at_zero = rna_layout::bulge_coordinates(stem, 0);
  const rna_layout::BulgePoints at_extra =
      rna_layout::bulge_coordinates_extra_distance(stem, 0, 2.0);

  // Only `at` (the peak) reads `extra_distance`; `prev`/`next` do not
  // (`boundingBoxes.inc:165-172`).
  expect_near(at_zero.prev.y, at_extra.prev.y, 1e-12, "extra_distance does not move 'prev'");
  expect_near(at_zero.next.y, at_extra.next.y, 1e-12, "extra_distance does not move 'next'");
  expect(at_extra.at.y > at_zero.at.y, "extra_distance moves 'at' further along +b");
  expect_near(at_extra.at.y - at_zero.at.y, 2.0, 1e-12,
              "extra_distance=2 moves 'at' by exactly 2 along +b (sign=+1, b=(0,1))");
}

void test_translate_moves_center_only() {
  LoopBox loop{Vec2{1.0, 2.0}, 5.0};
  rna_layout::translate_loop_box(loop, Vec2{10.0, -1.0});
  expect_near(loop.center.x, 11.0, 1e-12, "translate_loop_box moves center.x");
  expect_near(loop.center.y, 1.0, 1e-12, "translate_loop_box moves center.y");
  expect_near(loop.radius, 5.0, 1e-12, "translate_loop_box leaves radius untouched");

  StemBox stem;
  stem.c = Vec2{-3.0, 4.0};
  rna_layout::translate_stem_box(stem, Vec2{1.0, 1.0});
  expect_near(stem.c.x, -2.0, 1e-12, "translate_stem_box moves c.x");
  expect_near(stem.c.y, 5.0, 1e-12, "translate_stem_box moves c.y");
}

}  // namespace

int main() {
  test_initial_boxes_from_hairpin();
  test_compute_aabb_bounds_stem_and_loop();
  test_bulge_coordinates_extra_distance_offsets_along_b();
  test_translate_moves_center_only();

  if (g_failures > 0) {
    std::cerr << "bounding_boxes_test: " << g_failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
