/**
 * @file intersect_boxes.cpp
 * @brief Implements `intersect_tree.hpp`'s box-level (`StemBox`/`LoopBox`)
 *        predicates, ported from `intersectLevelBoundingBoxes.inc`.
 */

#include "rna_layout/bounding_boxes.hpp"
#include "rna_layout/geometry.hpp"
#include "rna_layout/intersect_tree.hpp"
#include "rna_layout/intersections.hpp"

namespace rna_layout {

bool intersect_stem_stem(const StemBox& stem1, const StemBox& stem2) {
  // Ported verbatim from `intersectStemStem` (`intersectLevelBoundingBoxes.inc:375-442`):
  // brute-force rectangle-vs-rectangle via each stem's two long sides.
  const Vec2 stem1_ea{stem1.e.x * stem1.a.x, stem1.e.x * stem1.a.y};
  const Vec2 stem1_eb{stem1.e.y * stem1.b.x, stem1.e.y * stem1.b.y};
  const Vec2 b1{stem1.c.x + stem1_ea.x + stem1_eb.x, stem1.c.y + stem1_ea.y + stem1_eb.y};
  const Vec2 c1{stem1.c.x + stem1_ea.x - stem1_eb.x, stem1.c.y + stem1_ea.y - stem1_eb.y};
  const Vec2 d1{stem1.c.x - stem1_ea.x - stem1_eb.x, stem1.c.y - stem1_ea.y - stem1_eb.y};
  const Vec2 a1{stem1.c.x - stem1_ea.x + stem1_eb.x, stem1.c.y - stem1_ea.y + stem1_eb.y};

  const Vec2 stem2_ea{stem2.e.x * stem2.a.x, stem2.e.x * stem2.a.y};
  const Vec2 stem2_eb{stem2.e.y * stem2.b.x, stem2.e.y * stem2.b.y};
  const Vec2 b2{stem2.c.x + stem2_ea.x + stem2_eb.x, stem2.c.y + stem2_ea.y + stem2_eb.y};
  const Vec2 c2{stem2.c.x + stem2_ea.x - stem2_eb.x, stem2.c.y + stem2_ea.y - stem2_eb.y};
  const Vec2 d2{stem2.c.x - stem2_ea.x - stem2_eb.x, stem2.c.y - stem2_ea.y - stem2_eb.y};
  const Vec2 a2{stem2.c.x - stem2_ea.x + stem2_eb.x, stem2.c.y - stem2_ea.y + stem2_eb.y};

  return geom::intersect_line_segments(a1, b1, a2, b2) ||
         geom::intersect_line_segments(a1, b1, c2, d2) ||
         geom::intersect_line_segments(c1, d1, a2, b2) ||
         geom::intersect_line_segments(c1, d1, c2, d2);
}

bool intersect_loop_loop(const LoopBox& loop1, const LoopBox& loop2, double clearance) {
  const double r1 = loop1.radius + 0.5 * geom::epsilon_recognize(clearance);
  const double r2 = loop2.radius + 0.5 * geom::epsilon_recognize(clearance);
  return geom::intersect_circle_circle(loop1.center, r1, loop2.center, r2);
}

bool intersect_stem_loop(const StemBox& stem, const LoopBox& loop, double clearance) {
  const Vec2 p = geom::closest_pt_point_obb(stem, loop.center);
  const Vec2 v_c_to_p = geom::vector_from_to(loop.center, p);
  const double distance_squared = geom::dot(v_c_to_p, v_c_to_p);
  const double grown_radius = loop.radius + geom::epsilon_recognize(clearance);
  return distance_squared < (grown_radius * grown_radius);
}

BulgeHit intersect_loop_bulges(const LoopBox& loop, const StemBox& stem, double clearance) {
  const Vec2 c = loop.center;
  const double r = loop.radius + geom::epsilon_recognize(clearance);

  for (int i = 0; i < static_cast<int>(stem.bulges.size()); ++i) {
    const BulgePoints points = bulge_coordinates(stem, i);
    if (geom::test_circle_triangle(c, r, points.prev, points.at, points.next)) {
      return BulgeHit{true, i};
    }
  }
  return BulgeHit{};
}

BulgeBulgeHit intersect_bulges_bulges(const StemBox& stem1, const StemBox& stem2,
                                      double clearance) {
  const double distance = 0.5 * geom::epsilon_recognize(clearance);

  for (int i = 0; i < static_cast<int>(stem1.bulges.size()); ++i) {
    const BulgePoints pi = bulge_coordinates_extra_distance(stem1, i, distance);
    for (int j = 0; j < static_cast<int>(stem2.bulges.size()); ++j) {
      const BulgePoints pj = bulge_coordinates_extra_distance(stem2, j, distance);
      if (geom::intersect_line_segments(pi.prev, pi.at, pj.prev, pj.at) ||
          geom::intersect_line_segments(pi.prev, pi.at, pj.at, pj.next) ||
          geom::intersect_line_segments(pi.at, pi.next, pj.prev, pj.at) ||
          geom::intersect_line_segments(pi.at, pi.next, pj.at, pj.next)) {
        return BulgeBulgeHit{true, i, j};
      }
    }
  }
  return BulgeBulgeHit{};
}

BulgeHit intersect_stem_bulges(const StemBox& stem1, const StemBox& stem2, double clearance) {
  if (stem2.bulges.empty()) {
    return BulgeHit{};
  }

  // N/E/S/W corners of `stem1`'s rectangle -- only the LEFT (`NW`/`SW`) and
  // RIGHT (`NE`/`SE`) sides are checked (the top/bottom sides border the
  // adjacent loops, not siblings -- see the vendored comment this preserves).
  const Vec2 p_nw{stem1.c.x + stem1.e.x * stem1.a.x - stem1.e.y * stem1.b.x,
                  stem1.c.y + stem1.e.x * stem1.a.y - stem1.e.y * stem1.b.y};
  const Vec2 p_sw{stem1.c.x - stem1.e.x * stem1.a.x - stem1.e.y * stem1.b.x,
                  stem1.c.y - stem1.e.x * stem1.a.y - stem1.e.y * stem1.b.y};
  const Vec2 p_ne{stem1.c.x + stem1.e.x * stem1.a.x + stem1.e.y * stem1.b.x,
                  stem1.c.y + stem1.e.x * stem1.a.y + stem1.e.y * stem1.b.y};
  const Vec2 p_se{stem1.c.x - stem1.e.x * stem1.a.x + stem1.e.y * stem1.b.x,
                  stem1.c.y - stem1.e.x * stem1.a.y + stem1.e.y * stem1.b.y};

  const double distance = geom::epsilon_recognize(clearance);

  for (int j = 0; j < static_cast<int>(stem2.bulges.size()); ++j) {
    const BulgePoints p = bulge_coordinates_extra_distance(stem2, j, distance);
    if (geom::intersect_line_segments(p_nw, p_sw, p.prev, p.at) ||
        geom::intersect_line_segments(p_nw, p_sw, p.at, p.next) ||
        geom::intersect_line_segments(p_ne, p_se, p.prev, p.at) ||
        geom::intersect_line_segments(p_ne, p_se, p.at, p.next)) {
      return BulgeHit{true, j};
    }
  }
  return BulgeHit{};
}

}  // namespace rna_layout
