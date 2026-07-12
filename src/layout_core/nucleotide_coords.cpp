/**
 * @file nucleotide_coords.cpp
 * @brief Implements `nucleotide_coords.hpp`.
 */

#include "rna_layout/nucleotide_coords.hpp"

#include <cmath>

#include "rna_layout/bounding_boxes.hpp"
#include "rna_layout/geometry.hpp"

namespace rna_layout {

namespace {

/**
 * Place every base of @p node's stem along its two long rectangle edges
 * (bulge bases at their own triangle peak instead). Ported verbatim from
 * `determineNucleotideCoordinates`'s "Handle stem of current node" block
 * (`RNApuzzler.c:200-265`), including its INDEX ARITHMETIC quirk: `coords`
 * is written at `[nt - 1]` throughout (the ordinary 1-indexed-`pair_table`-
 * to-0-indexed-`coords` offset) -- unlike `handle_loop`'s arcs, this block
 * uses that offset consistently.
 */
void handle_stem(const TreeNode& node, const std::vector<int>& pair_table, Coords& coords) {
  const StemBox& s_box = *node.sbox;  // NOLINT(bugprone-unchecked-optional-access) -- guarded by
                                      // `node.stem_start >= 1` at the call site (only a real stem
                                      // node reaches here, which always has an `sbox`).

  int left_bulges = 0;
  int right_bulges = 0;
  int current_bulge = 0;
  for (const Bulge& bulge : s_box.bulges) {
    if (bulge.sign < 0.0) {
      ++right_bulges;
    } else {
      ++left_bulges;
    }
  }

  // Left side: `node.stem_start` .. `node.loop_start`.
  int nt_start = node.stem_start;
  int nt_end = node.loop_start;
  const int left_segments = nt_end - nt_start - left_bulges;
  Vec2 p_start{s_box.c.x - s_box.e.x * s_box.a.x + s_box.e.y * s_box.b.x,
               s_box.c.y - s_box.e.x * s_box.a.y + s_box.e.y * s_box.b.y};
  Vec2 p_end{s_box.c.x + s_box.e.x * s_box.a.x + s_box.e.y * s_box.b.x,
             s_box.c.y + s_box.e.x * s_box.a.y + s_box.e.y * s_box.b.y};

  for (int nt = nt_start; nt < nt_end; ++nt) {
    if (pair_table[nt] == 0) {
      const Vec2 xy = bulge_coordinates(s_box, current_bulge).at;
      coords.x[nt - 1] = xy.x;
      coords.y[nt - 1] = xy.y;
      ++current_bulge;
    } else {
      coords.x[nt - 1] =
          p_start.x + (nt - nt_start - current_bulge) * (p_end.x - p_start.x) / left_segments;
      coords.y[nt - 1] =
          p_start.y + (nt - nt_start - current_bulge) * (p_end.y - p_start.y) / left_segments;
    }
  }
  coords.x[nt_end - 1] = p_end.x;
  coords.y[nt_end - 1] = p_end.y;

  // Right side: `pair_table[node.loop_start]` .. `pair_table[node.stem_start]`.
  nt_start = pair_table[node.loop_start];
  nt_end = pair_table[node.stem_start];
  const int right_segments = nt_end - nt_start - right_bulges;
  p_start = Vec2{s_box.c.x + s_box.e.x * s_box.a.x - s_box.e.y * s_box.b.x,
                 s_box.c.y + s_box.e.x * s_box.a.y - s_box.e.y * s_box.b.y};
  p_end = Vec2{s_box.c.x - s_box.e.x * s_box.a.x - s_box.e.y * s_box.b.x,
               s_box.c.y - s_box.e.x * s_box.a.y - s_box.e.y * s_box.b.y};

  for (int nt = nt_start; nt < nt_end; ++nt) {
    if (pair_table[nt] == 0) {
      const Vec2 xy = bulge_coordinates(s_box, current_bulge).at;
      coords.x[nt - 1] = xy.x;
      coords.y[nt - 1] = xy.y;
      ++current_bulge;
    } else {
      coords.x[nt - 1] = p_start.x + (nt - nt_start - current_bulge + left_bulges) *
                                         (p_end.x - p_start.x) / right_segments;
      coords.y[nt - 1] = p_start.y + (nt - nt_start - current_bulge + left_bulges) *
                                         (p_end.y - p_start.y) / right_segments;
    }
  }
  coords.x[nt_end - 1] = p_end.x;
  coords.y[nt_end - 1] = p_end.y;
}

/**
 * Place every base of @p node's loop-arc(s) on its bounding circle. Ported
 * verbatim from `determineNucleotideCoordinates`'s "loop" block
 * (`RNApuzzler.c:267-302`) -- INCLUDING its index arithmetic: `nt` starts at
 * `node.loop_start` (a 1-indexed `pair_table` position) and is used
 * DIRECTLY as a `coords` index (`coords.x[nt]`, not `coords.x[nt - 1]`) --
 * intentional, not an off-by-one bug: the loop's OWN opening base
 * (`loop_start`) was already written by `handle_stem`'s `coords.x[nt_end -
 * 1] = p_end.x` (where `nt_end == loop_start`), so this arc walk correctly
 * starts ONE PAST it, at nucleotide `loop_start + 1`.
 */
void handle_loop(const TreeNode& node, const std::vector<int>& pair_table, double paired,
                 Coords& coords) {
  if (!node.cfg.has_value()) {
    return;
  }
  const Config& cfg = *node.cfg;
  const Vec2 center = node.lbox->center;  // NOLINT(bugprone-unchecked-optional-access) -- a node
                                          // with `cfg` set always has `lbox`/`sbox` set too (see
                                          // `config_tree.cpp`'s invariant note).
  const double radius = cfg.radius;
  const double paired_angle = geom::distance_to_angle(radius, paired);

  const StemBox& s_box = *node.sbox;  // NOLINT(bugprone-unchecked-optional-access)
  double start_angle = std::atan2((s_box.c.y - center.y), (s_box.c.x - center.x));
  start_angle -= paired_angle / 2.0;

  int nt = node.loop_start;
  for (const ConfigArc& arc : cfg.arcs) {
    const int number_of_arc_segments = arc.segments;
    const double arc_angle = arc.angle;

    for (int arc_segment = 1; arc_segment < number_of_arc_segments; ++arc_segment) {
      const double angle =
          start_angle - arc_segment * ((arc_angle - paired_angle) / number_of_arc_segments);
      coords.x[nt] = center.x + radius * std::cos(angle);
      coords.y[nt] = center.y + radius * std::sin(angle);
      ++nt;
    }
    nt = pair_table[nt + 1];
    start_angle -= arc_angle;
  }
}

void determine_nucleotide_coordinates_recursive(const TreeNode& node,
                                                const std::vector<int>& pair_table, double paired,
                                                Coords& coords) {
  if (node.stem_start >= 1) {
    handle_stem(node, pair_table, coords);
  }
  handle_loop(node, pair_table, paired, coords);

  for (const auto& child : node.children) {
    determine_nucleotide_coordinates_recursive(*child, pair_table, paired, coords);
  }
}

/// The exterior (root-level unpaired) baseline walk. Ported verbatim from
/// `determineNucleotideCoordinates`'s final block (`RNApuzzler.c:312-329`);
/// see this file's header for why it runs once here, not once per node.
void set_exterior_coordinates(const std::vector<int>& pair_table, double unpaired, Coords& coords) {
  const int length = pair_table[0];

  coords.x[0] = geom::kExteriorY;
  coords.y[0] = geom::kExteriorY;

  int start = 1;
  if (pair_table[1] != 0) {
    start = pair_table[1] + 1;
  } else {
    start = 2;
  }

  for (int nt = start; nt <= length; ++nt) {
    if (pair_table[nt] == 0) {
      coords.x[nt - 1] = coords.x[nt - 2] + unpaired;
      coords.y[nt - 1] = geom::kExteriorY;
    } else {
      nt = pair_table[nt];
    }
  }
}

}  // namespace

void determine_nucleotide_coordinates(const TreeNode& root, const std::vector<int>& pair_table,
                                      double unpaired, double paired, Coords& coords) {
  determine_nucleotide_coordinates_recursive(root, pair_table, paired, coords);
  set_exterior_coordinates(pair_table, unpaired, coords);
}

}  // namespace rna_layout
