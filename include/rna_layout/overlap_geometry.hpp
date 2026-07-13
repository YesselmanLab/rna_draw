#pragma once

/**
 * @file overlap_geometry.hpp
 * @brief Pure distance/overlap math for the overlap checker, ported
 *        EXPRESSION-FOR-EXPRESSION from `rna_draw/geometry.py:99-286`.
 *
 * FIDELITY RULE (the whole point of this file, `.claude/plans/
 * current-plan-checker.md`): `rna_draw/geometry.py` is the frozen,
 * definitional arbiter -- this C++ twin's ONLY job is to agree with it
 * exactly. Every function below preserves the Python source's operator
 * order, branch order, and guard conditions verbatim (not just the
 * equivalent math): `std::hypot` where Python calls `hypot`, the same
 * `denom != 0.0` / `c_dot_c != 0.0` / `a_dot_a != 0.0` guards in the same
 * order, and `separation = required - depth` (never a `dist` shortcut) at
 * the call sites in `overlap_check.cpp`. A reordered op or a `sqrt`-for-
 * `hypot` swap is a port bug, not a stylistic choice -- see that plan's
 * Risk R1 and the differential parity harness
 * (`tests/test_overlap_native_parity.py`) that catches it.
 *
 * No primitive (`Disk`/`Capsule`/`PrimitiveId`) or spatial-hash logic
 * lives here -- only the geometry. `DiskGeom`/`SegmentGeom` are plain
 * unpacked-field views (center+radius / endpoints+half-width) so these
 * predicates take <=4 parameters without depending on `overlap_check.hpp`'s
 * primitive types (which additionally carry identity + exclusion `ends`).
 */

#include <optional>
#include <utility>

#include "rna_layout/types.hpp"

namespace rna_layout::overlap {

/// A disk's geometry only (no identity) -- mirrors `Disk.cx/cy/radius`
/// (`geometry.py:34-61`).
struct DiskGeom {
  double cx = 0.0;
  double cy = 0.0;
  double radius = 0.0;
};

/// A capsule's geometry only (no identity) -- mirrors
/// `Capsule.x0/y0/x1/y1/half_width` (`geometry.py:64-96`).
struct SegmentGeom {
  double x0 = 0.0;
  double y0 = 0.0;
  double x1 = 0.0;
  double y1 = 0.0;
  double half_width = 0.0;
};

/**
 * Clamp @p value to the closed interval `[lo, hi]`. Mirrors `clamp`
 * (`geometry.py:99-114`) branch-for-branch (not `std::clamp`, whose
 * argument order and tie-handling are not guaranteed to match).
 */
[[nodiscard]] double clamp(double value, double lo, double hi);

/**
 * Shortest distance between point @p p and segment `a`-`b`. Mirrors
 * `point_segment_distance` (`geometry.py:117-144`): the `denom == 0.0`
 * degenerate-segment special case first, then the clamped projection
 * parameter.
 */
[[nodiscard]] double point_segment_distance(Vec2 p, Vec2 a, Vec2 b);

/**
 * Clamped parametric solution for the closest points of two segments:
 * segment 1 is `a + s*ab`, segment 2 is `c + t*cd`, `s,t in [0,1]`.
 * Mirrors `_closest_params` (`geometry.py:147-188`) -- the Ericson
 * clamped-parametric method; the `denom != 0.0` / `c_dot_c != 0.0` /
 * `a_dot_a != 0.0` guards and their branch order are LOAD-BEARING (see
 * this file's header).
 *
 * @return `(s, t)`, both clamped to `[0, 1]`.
 */
[[nodiscard]] std::pair<double, double> closest_params(Vec2 a, Vec2 ab, Vec2 c, Vec2 cd);

/**
 * Shortest distance between segments `p1`-`q1` and `p2`-`q2`. Mirrors
 * `segment_segment_distance` (`geometry.py:191-228`): degenerate
 * (zero-length) segments are handled as special cases up front (both,
 * then segment 1, then segment 2) so the general `closest_params` solve
 * only ever runs on two genuine segments.
 */
[[nodiscard]] double segment_segment_distance(Vec2 p1, Vec2 q1, Vec2 p2, Vec2 q2);

/**
 * Test whether two disks overlap. Mirrors `disks_overlap`
 * (`geometry.py:231-248`): `dist = hypot(...)`, `required = r1 + r2`,
 * `dist < required - tol` (strict).
 *
 * @return Penetration depth `required - dist` if the disks overlap beyond
 *     @p tol, otherwise `std::nullopt`.
 */
[[nodiscard]] std::optional<double> disks_overlap(DiskGeom d1, DiskGeom d2, double tol);

/**
 * Test whether a disk overlaps a capsule. Mirrors `disk_capsule_overlap`
 * (`geometry.py:251-267`).
 *
 * @return Penetration depth if the disk overlaps the capsule beyond
 *     @p tol, otherwise `std::nullopt`.
 */
[[nodiscard]] std::optional<double> disk_capsule_overlap(DiskGeom d, SegmentGeom c, double tol);

/**
 * Test whether two capsules overlap. Mirrors `capsule_capsule_overlap`
 * (`geometry.py:270-286`).
 *
 * @return Penetration depth if the capsules overlap beyond @p tol,
 *     otherwise `std::nullopt`.
 */
[[nodiscard]] std::optional<double> capsule_capsule_overlap(SegmentGeom c1, SegmentGeom c2,
                                                            double tol);

}  // namespace rna_layout::overlap
