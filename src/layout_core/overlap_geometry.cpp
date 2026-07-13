/**
 * @file overlap_geometry.cpp
 * @brief Implements `overlap_geometry.hpp` -- port of `rna_draw/geometry.py:99-286`.
 *
 * See that header's file comment for the fidelity rule this whole file
 * follows: every expression below is transcribed operator-for-operator
 * from the Python source, not merely reproduced as equivalent math.
 */

#include "rna_layout/overlap_geometry.hpp"

#include <cmath>

namespace rna_layout::overlap {

double clamp(double value, double lo, double hi) {
  // Mirrors geometry.py:99-114's if/if/else branch order exactly (not
  // std::clamp, whose tie-handling contract differs).
  if (value < lo) {
    return lo;
  }
  if (value > hi) {
    return hi;
  }
  return value;
}

double point_segment_distance(Vec2 p, Vec2 a, Vec2 b) {
  // geometry.py:117-144.
  const double abx = b.x - a.x;
  const double aby = b.y - a.y;
  const double denom = abx * abx + aby * aby;
  if (denom == 0.0) {
    // Degenerate segment (a == b): distance reduces to point-point.
    return std::hypot(p.x - a.x, p.y - a.y);
  }
  double t = ((p.x - a.x) * abx + (p.y - a.y) * aby) / denom;
  t = clamp(t, 0.0, 1.0);
  return std::hypot(p.x - (a.x + t * abx), p.y - (a.y + t * aby));
}

std::pair<double, double> closest_params(Vec2 a, Vec2 ab, Vec2 c, Vec2 cd) {
  // geometry.py:147-188 (Ericson, Real-Time Collision Detection section
  // 5.1.9). `r` is `a - c`; the guard order and branch structure below are
  // LOAD-BEARING -- see overlap_geometry.hpp's file header.
  const double r_x = a.x - c.x;
  const double r_y = a.y - c.y;
  const double a_dot_a = ab.x * ab.x + ab.y * ab.y;
  const double c_dot_c = cd.x * cd.x + cd.y * cd.y;
  const double c_dot_r = cd.x * r_x + cd.y * r_y;
  const double a_dot_r = ab.x * r_x + ab.y * r_y;
  const double a_dot_c = ab.x * cd.x + ab.y * cd.y;

  const double denom = a_dot_a * c_dot_c - a_dot_c * a_dot_c;
  double s = 0.0;
  if (denom != 0.0) {
    s = clamp((a_dot_c * c_dot_r - a_dot_r * c_dot_c) / denom, 0.0, 1.0);
  }
  double t = (c_dot_c != 0.0) ? (a_dot_c * s + c_dot_r) / c_dot_c : 0.0;
  if (t < 0.0) {
    t = 0.0;
    s = (a_dot_a != 0.0) ? clamp(-a_dot_r / a_dot_a, 0.0, 1.0) : 0.0;
  } else if (t > 1.0) {
    t = 1.0;
    s = (a_dot_a != 0.0) ? clamp((a_dot_c - a_dot_r) / a_dot_a, 0.0, 1.0) : 0.0;
  }
  return {s, t};
}

double segment_segment_distance(Vec2 p1, Vec2 q1, Vec2 p2, Vec2 q2) {
  // geometry.py:191-228. Degenerate special-cases FIRST (both, seg1,
  // seg2), then the general clamped-parametric solve.
  const double abx = q1.x - p1.x;
  const double aby = q1.y - p1.y;
  const double cdx = q2.x - p2.x;
  const double cdy = q2.y - p2.y;
  const bool seg1_degenerate = abx == 0.0 && aby == 0.0;
  const bool seg2_degenerate = cdx == 0.0 && cdy == 0.0;

  if (seg1_degenerate && seg2_degenerate) {
    return std::hypot(p1.x - p2.x, p1.y - p2.y);
  }
  if (seg1_degenerate) {
    return point_segment_distance(p1, p2, q2);
  }
  if (seg2_degenerate) {
    return point_segment_distance(p2, p1, q1);
  }

  const auto [s, t] = closest_params(p1, Vec2{abx, aby}, p2, Vec2{cdx, cdy});
  const Vec2 closest1{p1.x + s * abx, p1.y + s * aby};
  const Vec2 closest2{p2.x + t * cdx, p2.y + t * cdy};
  return std::hypot(closest1.x - closest2.x, closest1.y - closest2.y);
}

std::optional<double> disks_overlap(DiskGeom d1, DiskGeom d2, double tol) {
  // geometry.py:231-248.
  const double dist = std::hypot(d1.cx - d2.cx, d1.cy - d2.cy);
  const double required = d1.radius + d2.radius;
  if (dist < required - tol) {
    return required - dist;
  }
  return std::nullopt;
}

std::optional<double> disk_capsule_overlap(DiskGeom d, SegmentGeom c, double tol) {
  // geometry.py:251-267.
  const double dist = point_segment_distance(Vec2{d.cx, d.cy}, Vec2{c.x0, c.y0}, Vec2{c.x1, c.y1});
  const double required = d.radius + c.half_width;
  if (dist < required - tol) {
    return required - dist;
  }
  return std::nullopt;
}

std::optional<double> capsule_capsule_overlap(SegmentGeom c1, SegmentGeom c2, double tol) {
  // geometry.py:270-286.
  const double dist = segment_segment_distance(Vec2{c1.x0, c1.y0}, Vec2{c1.x1, c1.y1},
                                               Vec2{c2.x0, c2.y0}, Vec2{c2.x1, c2.y1});
  const double required = c1.half_width + c2.half_width;
  if (dist < required - tol) {
    return required - dist;
  }
  return std::nullopt;
}

}  // namespace rna_layout::overlap
