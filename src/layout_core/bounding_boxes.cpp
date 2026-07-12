/**
 * @file bounding_boxes.cpp
 * @brief Implements `bounding_boxes.hpp`, ported from `boundingBoxes.inc`.
 */

#include "rna_layout/bounding_boxes.hpp"

#include <cmath>

#include "rna_layout/geometry.hpp"

namespace rna_layout {

namespace {

/*===========================================================================
 *  LoopBox construction
 *==========================================================================*/

/// The clockwise/counter-clockwise sense and center of the loop opened at
/// `start`, given the closing base pair's coordinates. Ported from
/// `getLoopData` (`boundingBoxes.inc:265`).
struct LoopData {
  Vec2 center;
  double radius = 0.0;
};

LoopData get_loop_data(int start, const std::vector<int>& pair_table,
                       const std::vector<BaseInfo>& base_info, const std::vector<Config>& configs,
                       const Coords& coords) {
  const int end = pair_table[start];
  // Same invariant/tradeoff as `turtle.cpp`'s `configs[base_info[start]
  // .loop_id.value()]`: callers only ever pass a `start` that is a real
  // (non-bulge-folded) loop-opening base, which always has `loop_id` set.
  const int loop_id =
      base_info[start].loop_id.value();  // NOLINT(bugprone-unchecked-optional-access)
  const double r = configs[loop_id].radius;

  // Reminder to offset -1 (`pair_table`/`base_info` are 1-indexed, `coords`
  // is 0-indexed) -- comment kept from `boundingBoxes.inc:280` since the
  // same off-by-one convention threads through the whole config-tree port.
  const Vec2 current{coords.x[start - 1], coords.y[start - 1]};
  const Vec2 next{coords.x[(start + 1) - 1], coords.y[(start + 1) - 1]};
  const Vec2 last{coords.x[end - 1], coords.y[end - 1]};
  const bool go_clockwise = geom::is_to_the_right_point_point(current, next, last);

  const Vec2 v_pair = geom::vector_from_to(last, current);
  const Vec2 v_normal = geom::normal(v_pair);
  const double pair_length = geom::length(v_pair);  // == paired
  const double center_dist = std::sqrt(r * r - 0.25 * pair_length * pair_length);

  const double dir = go_clockwise ? 1.0 : -1.0;
  LoopData data;
  data.center = Vec2{
      last.x + 0.5 * v_pair.x + dir * center_dist * v_normal.x,
      last.y + 0.5 * v_pair.y + dir * center_dist * v_normal.y,
  };
  data.radius = r;
  return data;
}

/*===========================================================================
 *  StemBox construction
 *==========================================================================*/

/// Ported from `createStemBox` (`boundingBoxes.inc:350`).
StemBox create_stem_box(Vec2 s, Vec2 e, Vec2 sp) {
  Vec2 a{0.5 * (e.x - s.x), 0.5 * (e.y - s.y)};
  const Vec2 b{0.5 * (s.x - sp.x), 0.5 * (s.y - sp.y)};

  double length_a = geom::length(a);
  const double length_b = geom::length(b);

  if (length_a == 0) {
    // Degenerate stem (zero backbone steps): solve using b's normal vector.
    // The reference scales the normal to length 0.1 and then immediately
    // divides by that same 0.1 below -- preserved as written (Task 1's
    // fidelity rule), not simplified to "a := normal(b)".
    a = geom::normal(b);
    length_a = 0.1;
    a.x = a.x * length_a;
    a.y = a.y * length_a;
  }

  StemBox box;
  box.a = Vec2{a.x / length_a, a.y / length_a};
  box.b = Vec2{b.x / length_b, b.y / length_b};
  box.c = Vec2{s.x + a.x - b.x, s.y + a.y - b.y};
  box.e = Vec2{length_a, length_b};
  return box;
}

/// Ported from `countBulges` (`boundingBoxes.inc:389`).
int count_bulges(const std::vector<int>& pair_table, int start, int end) {
  int bulge_count = 0;
  for (int i = start; i < end; ++i) {
    if (pair_table[i] == 0) {
      ++bulge_count;
    }
  }
  for (int i = pair_table[end]; i < pair_table[start]; ++i) {
    if (pair_table[i] == 0) {
      ++bulge_count;
    }
  }
  return bulge_count;
}

/// The stem-local coordinate of world point (x, y) along the stem's `a`
/// axis, used to place a bulge along that axis. Ported from `getA`
/// (`boundingBoxes.inc:412`).
double stem_box_a_coordinate(const StemBox& box, double x, double y) {
  const Vec2 p{x - box.c.x, y - box.c.y};
  if (box.b.x == 0.0) {
    return p.x / box.a.x;
  }
  if (box.b.y == 0.0) {
    return p.y / box.a.y;
  }
  return (p.x * box.b.y - p.y * box.b.x) / (box.a.x * box.b.y - box.a.y * box.b.x);
}

/// Ported from `createBulge` (`boundingBoxes.inc:445`); `i` is the
/// unpaired base's 1-indexed `pair_table` position.
Bulge create_bulge(const StemBox& box, const Coords& coords, int i, double sign) {
  // Remember -1 offset between pair_table and coords (`boundingBoxes.inc:454`).
  Bulge bulge;
  bulge.sign = sign;
  bulge.a_prev = stem_box_a_coordinate(box, coords.x[(i - 1) - 1], coords.y[(i - 1) - 1]);
  bulge.a_this = stem_box_a_coordinate(box, coords.x[(i - 1) + 0], coords.y[(i - 1) + 0]);
  bulge.a_next = stem_box_a_coordinate(box, coords.x[(i - 1) + 1], coords.y[(i - 1) + 1]);
  return bulge;
}

/// Ported from `setBulges` (`boundingBoxes.inc:468`).
void set_bulges(StemBox& box, const std::vector<int>& pair_table, int start, int end,
                const Coords& coords, int bulge_count, double bulge_dist) {
  box.bulge_dist = bulge_dist;
  if (bulge_count <= 0) {
    box.bulges.clear();
    return;
  }

  box.bulges.reserve(static_cast<std::size_t>(bulge_count));
  for (int i = start; i < end; ++i) {
    if (pair_table[i] == 0) {
      box.bulges.push_back(create_bulge(box, coords, i, /*sign=*/1.0));
    }
  }
  for (int i = pair_table[end]; i < pair_table[start]; ++i) {
    if (pair_table[i] == 0) {
      box.bulges.push_back(create_bulge(box, coords, i, /*sign=*/-1.0));
    }
  }
}

}  // namespace

LoopBox build_loop_box(int start, const std::vector<int>& pair_table,
                       const std::vector<BaseInfo>& base_info, const std::vector<Config>& configs,
                       const Coords& coords) {
  const LoopData data = get_loop_data(start, pair_table, base_info, configs, coords);
  return LoopBox{data.center, data.radius};
}

StemBox build_stem_box(int start, int end, const std::vector<int>& pair_table, const Coords& coords,
                       double bulge_dist) {
  const int i_s = start;
  const int i_e = end;
  const int i_sp = pair_table[start];

  // Coordinates for the rectangle corners; -1 for the pair_table/coords
  // offset (`boundingBoxes.inc:523-535`).
  const Vec2 s{coords.x[i_s - 1], coords.y[i_s - 1]};
  const Vec2 e{coords.x[i_e - 1], coords.y[i_e - 1]};
  const Vec2 sp{coords.x[i_sp - 1], coords.y[i_sp - 1]};

  StemBox box = create_stem_box(s, e, sp);

  const int bulge_count = count_bulges(pair_table, i_s, i_e);
  set_bulges(box, pair_table, i_s, i_e, coords, bulge_count, bulge_dist);
  return box;
}

Aabb compute_aabb(const StemBox& stem_box, const LoopBox& loop_box) {
  const Vec2 stem_ea{stem_box.e.x * stem_box.a.x, stem_box.e.x * stem_box.a.y};
  const Vec2 stem_eb{stem_box.e.y * stem_box.b.x, stem_box.e.y * stem_box.b.y};

  std::vector<Vec2> points;
  points.reserve(6 + stem_box.bulges.size());
  // Corners of the stem.
  points.push_back(
      Vec2{stem_box.c.x - stem_ea.x + stem_eb.x, stem_box.c.y - stem_ea.y + stem_eb.y});
  points.push_back(
      Vec2{stem_box.c.x + stem_ea.x + stem_eb.x, stem_box.c.y + stem_ea.y + stem_eb.y});
  points.push_back(
      Vec2{stem_box.c.x + stem_ea.x - stem_eb.x, stem_box.c.y + stem_ea.y - stem_eb.y});
  points.push_back(
      Vec2{stem_box.c.x - stem_ea.x - stem_eb.x, stem_box.c.y - stem_ea.y - stem_eb.y});
  // Lower-left / upper-right of the loop AABB.
  points.push_back(Vec2{loop_box.center.x - loop_box.radius, loop_box.center.y - loop_box.radius});
  points.push_back(Vec2{loop_box.center.x + loop_box.radius, loop_box.center.y + loop_box.radius});
  // Bulge peaks.
  for (std::size_t i = 0; i < stem_box.bulges.size(); ++i) {
    points.push_back(bulge_coordinates(stem_box, static_cast<int>(i)).at);
  }

  Aabb aabb{points[0], points[0]};
  for (std::size_t i = 1; i < points.size(); ++i) {
    aabb.min.x = std::fmin(aabb.min.x, points[i].x);
    aabb.min.y = std::fmin(aabb.min.y, points[i].y);
    aabb.max.x = std::fmax(aabb.max.x, points[i].x);
    aabb.max.y = std::fmax(aabb.max.y, points[i].y);
  }
  return aabb;
}

BulgePoints bulge_coordinates_extra_distance(const StemBox& stem, int index,
                                             double extra_distance) {
  const Bulge& bulge = stem.bulges[static_cast<std::size_t>(index)];
  BulgePoints points;
  points.prev = Vec2{
      stem.c.x + bulge.a_prev * stem.a.x + bulge.sign * stem.b.x * stem.e.y,
      stem.c.y + bulge.a_prev * stem.a.y + bulge.sign * stem.b.y * stem.e.y,
  };
  points.at = Vec2{
      stem.c.x + bulge.a_this * stem.a.x +
          bulge.sign * stem.b.x * (stem.e.y + extra_distance + stem.bulge_dist),
      stem.c.y + bulge.a_this * stem.a.y +
          bulge.sign * stem.b.y * (stem.e.y + extra_distance + stem.bulge_dist),
  };
  points.next = Vec2{
      stem.c.x + bulge.a_next * stem.a.x + bulge.sign * stem.b.x * stem.e.y,
      stem.c.y + bulge.a_next * stem.a.y + bulge.sign * stem.b.y * stem.e.y,
  };
  return points;
}

BulgePoints bulge_coordinates(const StemBox& stem, int index) {
  return bulge_coordinates_extra_distance(stem, index, 0.0);
}

void translate_loop_box(LoopBox& box, Vec2 vector) {
  box.center = Vec2{box.center.x + vector.x, box.center.y + vector.y};
}

void translate_stem_box(StemBox& box, Vec2 vector) {
  box.c = Vec2{box.c.x + vector.x, box.c.y + vector.y};
}

}  // namespace rna_layout
