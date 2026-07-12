// ctest: `rna_layout`'s box-level and node-level intersection DETECTION
// (`intersect_tree.hpp`), Milestone A step 5.
//
// Uses `clearance = 0.0` throughout (`geom::epsilon_recognize(0.0) == 0`)
// unless a test specifically exercises the clearance growth, so expected
// true/false results match plain Euclidean geometry directly -- the
// vendored stock `epsilonRecognize` (14 units at `clearance = 1.0`) is
// large relative to this file's small hand-built fixtures, which would
// otherwise force every "far apart" negative case to also place objects
// implausibly far apart just to clear that margin.

#include "rna_layout/intersect_tree.hpp"

#include <cstdlib>
#include <iostream>
#include <memory>

#include "rna_layout/bounding_boxes.hpp"
#include "rna_layout/tree.hpp"
#include "rna_layout/types.hpp"
#include "test_tree_helper.hpp"

using rna_layout::Aabb;
using rna_layout::Bulge;
using rna_layout::IntersectionType;
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

void expect_type(IntersectionType actual, IntersectionType expected, const char* what) {
  if (actual != expected) {
    std::cerr << "FAIL " << what << ": expected type " << static_cast<int>(expected) << ", got "
              << static_cast<int>(actual) << "\n";
    ++g_failures;
  }
}

/*=============================================================================
 *  Box-level predicates
 *============================================================================*/

StemBox make_stem(Vec2 c, Vec2 e) {
  StemBox stem;
  stem.a = Vec2{1.0, 0.0};
  stem.b = Vec2{0.0, 1.0};
  stem.c = c;
  stem.e = e;
  return stem;
}

void test_intersect_stem_stem() {
  const StemBox overlapping1 = make_stem(Vec2{0.0, 0.0}, Vec2{5.0, 2.0});
  const StemBox overlapping2 = make_stem(Vec2{3.0, 0.0}, Vec2{5.0, 2.0});
  expect(rna_layout::intersect_stem_stem(overlapping1, overlapping2),
         "intersect_stem_stem: overlapping rectangles");

  const StemBox far = make_stem(Vec2{1000.0, 1000.0}, Vec2{5.0, 2.0});
  expect(!rna_layout::intersect_stem_stem(overlapping1, far),
         "intersect_stem_stem: far-apart rectangles do not intersect");
}

void test_intersect_loop_loop() {
  const LoopBox loop1{Vec2{0.0, 0.0}, 5.0};
  const LoopBox loop2{Vec2{8.0, 0.0}, 5.0};
  expect(rna_layout::intersect_loop_loop(loop1, loop2, /*clearance=*/0.0),
         "intersect_loop_loop: overlapping circles (distance 8 < r1+r2=10)");

  const LoopBox loop3{Vec2{1000.0, 0.0}, 5.0};
  expect(!rna_layout::intersect_loop_loop(loop1, loop3, /*clearance=*/0.0),
         "intersect_loop_loop: far-apart circles do not overlap");
}

void test_intersect_stem_loop() {
  const StemBox stem = make_stem(Vec2{0.0, 0.0}, Vec2{2.0, 1.0});  // spans x[-2,2], y[-1,1]
  const LoopBox near_loop{Vec2{2.5, 0.0}, 1.0};  // circle touches/overlaps the stem's right edge
  expect(rna_layout::intersect_stem_loop(stem, near_loop, /*clearance=*/0.0),
         "intersect_stem_loop: loop overlapping the stem rectangle");

  const LoopBox far_loop{Vec2{1000.0, 1000.0}, 1.0};
  expect(!rna_layout::intersect_stem_loop(stem, far_loop, /*clearance=*/0.0),
         "intersect_stem_loop: far-away loop does not intersect");
}

/// A stem with one bulge whose peak is at `(0, 1.5)`: `bulge_dist = 0.5`
/// (so `at` pokes 0.5 further than `prev`/`next`, which sit at `e.y = 1.0`)
/// -- matching the vendored `distBulge` invariant (always > 0 in real
/// usage; `bounding_boxes_test.cpp` covers `bulge_coordinates` itself in
/// depth, this fixture just needs A real, non-degenerate triangle).
StemBox make_stem_with_one_bulge() {
  StemBox stem = make_stem(Vec2{0.0, 0.0}, Vec2{5.0, 1.0});
  stem.bulge_dist = 0.5;
  Bulge bulge;
  bulge.sign = 1.0;
  bulge.a_prev = -1.0;
  bulge.a_this = 0.0;
  bulge.a_next = 1.0;
  stem.bulges.push_back(bulge);
  return stem;
}

void test_intersect_loop_bulges() {
  const StemBox stem = make_stem_with_one_bulge();
  const rna_layout::BulgePoints points = rna_layout::bulge_coordinates(stem, 0);
  expect(points.at.x == 0.0 && points.at.y == 1.5, "setup: bulge peak is at (0, 1.5)");

  const LoopBox at_peak{points.at, 0.1};
  const rna_layout::BulgeHit hit =
      rna_layout::intersect_loop_bulges(at_peak, stem, /*clearance=*/0.0);
  expect(hit.intersects, "intersect_loop_bulges: loop centered on the bulge peak intersects it");
  expect(hit.bulge == 0, "intersect_loop_bulges: reports the intersecting bulge's index");

  const LoopBox far{Vec2{1000.0, 1000.0}, 0.1};
  expect(!rna_layout::intersect_loop_bulges(far, stem, /*clearance=*/0.0).intersects,
         "intersect_loop_bulges: far-away loop does not intersect");
}

void test_intersect_bulges_bulges() {
  const StemBox stem1 = make_stem_with_one_bulge();
  StemBox stem2 = make_stem_with_one_bulge();
  stem2.c = Vec2{0.0, 0.0};  // identical bulge peak -> triangles coincide

  const rna_layout::BulgeBulgeHit hit =
      rna_layout::intersect_bulges_bulges(stem1, stem2, /*clearance=*/0.0);
  expect(hit.intersects, "intersect_bulges_bulges: identical bulge triangles intersect");

  StemBox stem3 = make_stem_with_one_bulge();
  stem3.c = Vec2{1000.0, 1000.0};
  expect(!rna_layout::intersect_bulges_bulges(stem1, stem3, /*clearance=*/0.0).intersects,
         "intersect_bulges_bulges: far-apart bulge triangles do not intersect");
}

void test_intersect_stem_bulges() {
  const StemBox stem1 = make_stem(Vec2{0.0, 0.0}, Vec2{5.0, 1.0});  // spans x[-5,5], y[-1,1]
  // stem2's bulge peak lands just inside stem1's left edge (x=-5).
  StemBox stem2 = make_stem_with_one_bulge();
  stem2.c = Vec2{-5.0, 0.0};

  const rna_layout::BulgeHit hit =
      rna_layout::intersect_stem_bulges(stem1, stem2, /*clearance=*/0.0);
  expect(hit.intersects, "intersect_stem_bulges: a bulge crossing the stem's side intersects");

  StemBox stem3 = make_stem_with_one_bulge();
  stem3.c = Vec2{1000.0, 1000.0};
  expect(!rna_layout::intersect_stem_bulges(stem1, stem3, /*clearance=*/0.0).intersects,
         "intersect_stem_bulges: far-away bulge does not intersect");
}

/*=============================================================================
 *  Node-level dispatch (`intersect_node_node`)
 *============================================================================*/

/// A standalone, non-root node (via a throwaway parent) with the given
/// stem/loop boxes -- enough state for `intersect_node_node` without
/// building a whole real config tree.
struct StandaloneNode {
  std::unique_ptr<TreeNode> parent = std::make_unique<TreeNode>();
  TreeNode* node = nullptr;
};

StandaloneNode make_standalone_node(const StemBox& sbox, const LoopBox& lbox) {
  StandaloneNode result;
  result.node = result.parent->add_child();
  result.node->sbox = sbox;
  result.node->lbox = lbox;
  result.node->aabb = rna_layout::compute_aabb(sbox, lbox);
  return result;
}

void test_intersect_node_node_no_intersection() {
  const StandaloneNode a = make_standalone_node(make_stem(Vec2{0.0, 0.0}, Vec2{2.0, 1.0}),
                                                LoopBox{Vec2{0.0, 20.0}, 2.0});
  const StandaloneNode b = make_standalone_node(make_stem(Vec2{1000.0, 1000.0}, Vec2{2.0, 1.0}),
                                                LoopBox{Vec2{1000.0, 1020.0}, 2.0});
  const rna_layout::NodeIntersection result =
      rna_layout::intersect_node_node(*a.node, *b.node, /*clearance=*/0.0);
  expect_type(result.type, IntersectionType::none, "intersect_node_node: far-apart nodes");
}

void test_intersect_node_node_stem_stem() {
  // Different (unrelated) parents, overlapping stems, loops placed well
  // away from either stem so the AABB reject/SxS check is what fires.
  const StandaloneNode a = make_standalone_node(make_stem(Vec2{0.0, 0.0}, Vec2{5.0, 2.0}),
                                                LoopBox{Vec2{0.0, 100.0}, 2.0});
  const StandaloneNode b = make_standalone_node(make_stem(Vec2{3.0, 0.0}, Vec2{5.0, 2.0}),
                                                LoopBox{Vec2{3.0, 100.0}, 2.0});
  const rna_layout::NodeIntersection result =
      rna_layout::intersect_node_node(*a.node, *b.node, /*clearance=*/0.0);
  expect_type(result.type, IntersectionType::stem_stem,
              "intersect_node_node: overlapping stems, unrelated parents -> SxS");
}

void test_intersect_node_node_loop_loop() {
  // Stems far apart (no SxS); loops overlapping.
  const StandaloneNode a = make_standalone_node(make_stem(Vec2{0.0, -100.0}, Vec2{2.0, 1.0}),
                                                LoopBox{Vec2{0.0, 0.0}, 5.0});
  const StandaloneNode b = make_standalone_node(make_stem(Vec2{50.0, -100.0}, Vec2{2.0, 1.0}),
                                                LoopBox{Vec2{8.0, 0.0}, 5.0});
  const rna_layout::NodeIntersection result =
      rna_layout::intersect_node_node(*a.node, *b.node, /*clearance=*/0.0);
  expect_type(result.type, IntersectionType::loop_loop,
              "intersect_node_node: overlapping loops, non-overlapping stems -> LxL");
}

void test_intersect_node_node_parent_child_excludes_stem_stem_and_loop_loop() {
  auto parent = std::make_unique<TreeNode>();
  TreeNode* p = parent->add_child();
  p->sbox = make_stem(Vec2{0.0, 0.0}, Vec2{5.0, 2.0});
  p->lbox = LoopBox{Vec2{0.0, 0.0}, 5.0};
  p->aabb = rna_layout::compute_aabb(*p->sbox, *p->lbox);

  TreeNode* c = p->add_child();          // c.parent == p
  c->sbox = *p->sbox;                    // NOLINT(bugprone-unchecked-optional-access) --
  c->lbox = *p->lbox;                    // NOLINT(bugprone-unchecked-optional-access) -- set two
  c->aabb = p->aabb;                     // lines above; identical (overlapping) boxes on purpose

  const rna_layout::NodeIntersection result =
      rna_layout::intersect_node_node(*p, *c, /*clearance=*/0.0);
  expect(result.type != IntersectionType::stem_stem,
         "intersect_node_node: parent-child excludes SxS");
  expect(result.type != IntersectionType::loop_loop,
         "intersect_node_node: parent-child excludes LxL");
}

/*=============================================================================
 *  Whole-tree detection set (`detect_intersections`)
 *============================================================================*/

void test_detect_intersections_is_empty_on_a_clean_structure() {
  // A well-formed multi-branch structure: the turtle/config-tree pipeline
  // fits radii so siblings do not overlap by construction -- a real
  // (non-degenerate) regression check that `detect_intersections` does not
  // false-positive on ordinary geometry.
  const auto tree =
      rna_layout::test::build_updated_tree("((((...)))(((...)))(((...))))", 35.0, 25.0);
  const std::vector<rna_layout::Detection> detections =
      rna_layout::detect_intersections(*tree, /*clearance=*/1.0);
  expect(detections.empty(),
         "detect_intersections: a clean multi-branch structure has zero detections");
}

}  // namespace

int main() {
  test_intersect_stem_stem();
  test_intersect_loop_loop();
  test_intersect_stem_loop();
  test_intersect_loop_bulges();
  test_intersect_bulges_bulges();
  test_intersect_stem_bulges();

  test_intersect_node_node_no_intersection();
  test_intersect_node_node_stem_stem();
  test_intersect_node_node_loop_loop();
  test_intersect_node_node_parent_child_excludes_stem_stem_and_loop_loop();

  test_detect_intersections_is_empty_on_a_clean_structure();

  if (g_failures > 0) {
    std::cerr << "intersect_tree_test: " << g_failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
