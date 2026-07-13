/**
 * @file overlap_check.cpp
 * @brief Implements `overlap_check.hpp` -- port of `rna_draw/overlap.py:36-468`.
 *
 * See that header's file comment for the R1 PrimId-ordering fix and the
 * pinned `OverlapKind` index mapping this file relies on.
 */

#include "rna_layout/overlap_check.hpp"

#include <algorithm>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <utility>

#include "rna_layout/overlap_geometry.hpp"

namespace rna_layout::overlap {

// ---------------------------------------------------------------------
// PrimId
// ---------------------------------------------------------------------

std::string kind_str(PrimKind kind) {
  switch (kind) {
    case PrimKind::bb:
      return "bb";
    case PrimKind::nt:
      return "nt";
    case PrimKind::pair:
      return "pair";
  }
  return "";  // unreachable; silences -Wreturn-type on some compilers
}

bool operator==(const PrimId& lhs, const PrimId& rhs) {
  return lhs.kind == rhs.kind && lhs.index == rhs.index;
}

bool operator<(const PrimId& lhs, const PrimId& rhs) {
  // `PrimKind`'s enum-declaration order already reproduces Python's
  // kind-STRING sort order ("bb" < "nt" < "pair") -- see overlap_check.hpp's
  // file header (R1). Comparing the enum values directly is therefore
  // equivalent to Python's `sorted((a.pid, b.pid))` string comparison.
  if (lhs.kind != rhs.kind) {
    return lhs.kind < rhs.kind;
  }
  return lhs.index < rhs.index;
}

// ---------------------------------------------------------------------
// Ends
// ---------------------------------------------------------------------

bool Ends::contains(int index) const { return a == index || b == index; }

bool Ends::intersects(const Ends& other) const {
  if (contains(other.a)) {
    return true;
  }
  return other.b != -1 && contains(other.b);
}

// ---------------------------------------------------------------------
// Primitive (Disk | Capsule) accessors
// ---------------------------------------------------------------------

PrimId pid_of(const Primitive& primitive) {
  return std::visit([](const auto& shape) { return shape.pid; }, primitive);
}

Ends ends_of(const Primitive& primitive) {
  return std::visit([](const auto& shape) { return shape.ends; }, primitive);
}

Aabb aabb_of(const Primitive& primitive) {
  // geometry.py:54-61 (Disk.aabb) / geometry.py:88-96 (Capsule.aabb).
  if (const auto* disk = std::get_if<Disk>(&primitive)) {
    return Aabb{Vec2{disk->cx - disk->radius, disk->cy - disk->radius},
                Vec2{disk->cx + disk->radius, disk->cy + disk->radius}};
  }
  const auto& capsule = std::get<Capsule>(primitive);
  const double w = capsule.half_width;
  return Aabb{Vec2{std::min(capsule.x0, capsule.x1) - w, std::min(capsule.y0, capsule.y1) - w},
              Vec2{std::max(capsule.x0, capsule.x1) + w, std::max(capsule.y0, capsule.y1) + w}};
}

// ---------------------------------------------------------------------
// OverlapKind / Witness / OverlapResult
// ---------------------------------------------------------------------

std::string kind_str(OverlapKind kind) {
  switch (kind) {
    case OverlapKind::disk_disk:
      return "disk_disk";
    case OverlapKind::disk_capsule:
      return "disk_capsule";
    case OverlapKind::capsule_capsule:
      return "capsule_capsule";
  }
  return "";  // unreachable
}

int OverlapResult::num_overlaps() const { return static_cast<int>(witnesses.size()); }

bool OverlapResult::passed() const { return witnesses.empty(); }

// ---------------------------------------------------------------------
// build_disks / build_backbone_capsules / build_pair_capsules / build_primitives
// ---------------------------------------------------------------------

std::vector<Disk> build_disks(const std::vector<double>& x, const std::vector<double>& y,
                              double node_r) {
  // overlap.py:120-134.
  std::vector<Disk> disks;
  disks.reserve(x.size());
  for (std::size_t i = 0; i < x.size(); ++i) {
    disks.push_back(Disk{PrimId{PrimKind::nt, static_cast<int>(i)}, x[i], y[i], node_r,
                         Ends{static_cast<int>(i), -1}});
  }
  return disks;
}

std::vector<Capsule> build_backbone_capsules(const std::vector<double>& x,
                                             const std::vector<double>& y, double half_width) {
  // overlap.py:137-161.
  std::vector<Capsule> capsules;
  if (x.empty()) {
    return capsules;
  }
  capsules.reserve(x.size() - 1);
  for (std::size_t i = 0; i + 1 < x.size(); ++i) {
    capsules.push_back(Capsule{PrimId{PrimKind::bb, static_cast<int>(i)}, x[i], y[i], x[i + 1],
                               y[i + 1], half_width,
                               Ends{static_cast<int>(i), static_cast<int>(i + 1)}});
  }
  return capsules;
}

std::vector<Capsule> build_pair_capsules(const std::vector<double>& x, const std::vector<double>& y,
                                         const std::vector<int>& pair_map, double half_width) {
  // overlap.py:164-191.
  std::vector<Capsule> capsules;
  for (std::size_t i = 0; i < pair_map.size(); ++i) {
    const int j = pair_map[i];
    if (j != -1 && static_cast<int>(i) < j) {
      capsules.push_back(Capsule{PrimId{PrimKind::pair, static_cast<int>(i)}, x[i], y[i], x[j],
                                 y[j], half_width, Ends{static_cast<int>(i), j}});
    }
  }
  return capsules;
}

std::vector<Primitive> build_primitives(const std::vector<double>& x, const std::vector<double>& y,
                                        const std::vector<int>& pair_map,
                                        const OverlapParams& params) {
  // overlap.py:194-216: disks, THEN backbone, THEN pairs -- index stability.
  // Disk/Capsule are small, trivially-copyable value types, so this copies
  // each element into the variant rather than moving (a move would be a
  // no-op here; see this project's clang-tidy config).
  const std::vector<Disk> disks = build_disks(x, y, params.node_r);
  const std::vector<Capsule> backbone = build_backbone_capsules(x, y, params.backbone_half_width);
  const std::vector<Capsule> pairs = build_pair_capsules(x, y, pair_map, params.pair_half_width);

  std::vector<Primitive> primitives;
  primitives.reserve(disks.size() + backbone.size() + pairs.size());
  for (const Disk& disk : disks) {
    primitives.emplace_back(disk);
  }
  for (const Capsule& capsule : backbone) {
    primitives.emplace_back(capsule);
  }
  for (const Capsule& capsule : pairs) {
    primitives.emplace_back(capsule);
  }
  return primitives;
}

// ---------------------------------------------------------------------
// is_excluded
// ---------------------------------------------------------------------

bool is_excluded(const Primitive& a, const Primitive& b, const std::vector<int>& pair_map) {
  // overlap.py:219-241.
  const auto* disk_a = std::get_if<Disk>(&a);
  const auto* disk_b = std::get_if<Disk>(&b);
  if (disk_a != nullptr && disk_b != nullptr) {
    const int i = disk_a->ends.a;
    const int j = disk_b->ends.a;
    return std::abs(i - j) == 1 || pair_map[static_cast<std::size_t>(i)] == j;
  }
  return ends_of(a).intersects(ends_of(b));
}

// ---------------------------------------------------------------------
// test_pair, split into per-kind helpers (style: <=30 lines, one job each)
// ---------------------------------------------------------------------

namespace {

/// Normalize a mixed pair so the `Disk` comes first. Mirrors `_ordered`
/// (`overlap.py:244-257`).
std::pair<const Primitive*, const Primitive*> ordered(const Primitive& a, const Primitive& b) {
  if (std::holds_alternative<Capsule>(a) && std::holds_alternative<Disk>(b)) {
    return {&b, &a};
  }
  return {&a, &b};
}

/// Shared witness-building tail of `test_pair` (`overlap.py:292-297`):
/// `None` depth -> `std::nullopt`; else sort the two ids and compute
/// `separation = required - depth` (never a `dist` shortcut -- see
/// `overlap_check.hpp`'s file header).
std::optional<Witness> build_witness(OverlapKind kind, PrimId pid_a, PrimId pid_b, double required,
                                     std::optional<double> depth) {
  if (!depth.has_value()) {
    return std::nullopt;
  }
  const PrimId id_a = pid_a < pid_b ? pid_a : pid_b;
  const PrimId id_b = pid_a < pid_b ? pid_b : pid_a;
  return Witness{kind, id_a, id_b, required - *depth, *depth};
}

std::optional<Witness> test_disk_disk(const Disk& a, const Disk& b, double tol) {
  const double required = a.radius + b.radius;
  const std::optional<double> depth =
      disks_overlap(DiskGeom{a.cx, a.cy, a.radius}, DiskGeom{b.cx, b.cy, b.radius}, tol);
  return build_witness(OverlapKind::disk_disk, a.pid, b.pid, required, depth);
}

std::optional<Witness> test_disk_capsule(const Disk& a, const Capsule& b, double tol) {
  const double required = a.radius + b.half_width;
  const std::optional<double> depth = disk_capsule_overlap(
      DiskGeom{a.cx, a.cy, a.radius}, SegmentGeom{b.x0, b.y0, b.x1, b.y1, b.half_width}, tol);
  return build_witness(OverlapKind::disk_capsule, a.pid, b.pid, required, depth);
}

std::optional<Witness> test_capsule_capsule(const Capsule& a, const Capsule& b, double tol) {
  const double required = a.half_width + b.half_width;
  const std::optional<double> depth =
      capsule_capsule_overlap(SegmentGeom{a.x0, a.y0, a.x1, a.y1, a.half_width},
                              SegmentGeom{b.x0, b.y0, b.x1, b.y1, b.half_width}, tol);
  return build_witness(OverlapKind::capsule_capsule, a.pid, b.pid, required, depth);
}

}  // namespace

std::optional<Witness> test_pair(const Primitive& a, const Primitive& b, double tol) {
  // overlap.py:260-297, dispatch-only (the three kind branches are the
  // helpers above).
  const auto [first, second] = ordered(a, b);
  if (const auto* disk_first = std::get_if<Disk>(first)) {
    if (const auto* disk_second = std::get_if<Disk>(second)) {
      return test_disk_disk(*disk_first, *disk_second, tol);
    }
    return test_disk_capsule(*disk_first, std::get<Capsule>(*second), tol);
  }
  // `ordered` guarantees a Disk-first for any mixed pair, so reaching
  // here means both are capsules (the only remaining `Primitive` case).
  return test_capsule_capsule(std::get<Capsule>(*first), std::get<Capsule>(*second), tol);
}

// ---------------------------------------------------------------------
// check_overlaps / check_overlaps_bruteforce
// ---------------------------------------------------------------------

namespace {

/// Mirrors `_validate_inputs` (`overlap.py:300-319`) -- PLAN-CRITIC R2:
/// raises rather than leaving a degenerate call UB-adjacent. Also bounds-
/// checks `pair_map` entries before indexing (Python would let an
/// out-of-range partner raise `IndexError` from the list indexing itself;
/// here that is folded into the same `std::invalid_argument` contract to
/// avoid undefined behavior in the C++ port).
void validate_inputs(const std::vector<double>& x, const std::vector<double>& y,
                     const std::vector<int>& pair_map) {
  if (!(x.size() == y.size() && y.size() == pair_map.size())) {
    throw std::invalid_argument("x, y, and pair_map must all have equal length");
  }
  if (x.empty()) {
    throw std::invalid_argument("x, y, and pair_map must be non-empty");
  }
  const auto n = static_cast<int>(pair_map.size());
  for (std::size_t i = 0; i < pair_map.size(); ++i) {
    const int j = pair_map[i];
    if (j == -1) {
      continue;
    }
    if (j < 0 || j >= n) {
      throw std::invalid_argument("pair_map entry out of range at index " + std::to_string(i));
    }
    if (pair_map[static_cast<std::size_t>(j)] != static_cast<int>(i)) {
      throw std::invalid_argument("pair_map is not symmetric at index " + std::to_string(i));
    }
  }
}

/// Insert one primitive into the spatial hash. Mirrors `_insert_primitive`
/// (`overlap.py:388-406`): capsules are tiled segments, disks a single AABB.
void insert_primitive(SpatialHash& grid, int index, const Primitive& primitive) {
  if (const auto* capsule = std::get_if<Capsule>(&primitive)) {
    grid.insert_segment(index, Vec2{capsule->x0, capsule->y0}, Vec2{capsule->x1, capsule->y1},
                        capsule->half_width);
  } else {
    grid.insert(index, aabb_of(primitive));
  }
}

/// Test every candidate index pair, skipping excluded pairs. Mirrors
/// `_collect_witnesses` (`overlap.py:322-351`).
std::vector<Witness> collect_witnesses(const std::vector<Primitive>& primitives,
                                       const std::vector<std::pair<int, int>>& candidate_pairs,
                                       const std::vector<int>& pair_map, double tol) {
  std::vector<Witness> witnesses;
  for (const auto& [i, j] : candidate_pairs) {
    const Primitive& a = primitives[static_cast<std::size_t>(i)];
    const Primitive& b = primitives[static_cast<std::size_t>(j)];
    if (is_excluded(a, b, pair_map)) {
      continue;
    }
    if (std::optional<Witness> witness = test_pair(a, b, tol)) {
      witnesses.push_back(*witness);
    }
  }
  return witnesses;
}

/// Deterministic ordering key for the report's witness list (kind, id_a,
/// id_b). NOTE: this ordering is NOT load-bearing for parity -- the
/// differential harness compares witnesses as an unordered key SET (see
/// `tests/test_overlap_native_parity.py`), matching Python's own tests'
/// `set(report.witnesses)` comparisons. Any stable order is therefore
/// correct here; `OverlapKind`'s pinned integer values are used (not
/// Python's incidental `str`-subclass sort, which orders by kind STRING
/// value and would sort differently -- irrelevant since order is unchecked).
bool witness_order(const Witness& lhs, const Witness& rhs) {
  if (lhs.kind != rhs.kind) {
    return lhs.kind < rhs.kind;
  }
  if (!(lhs.id_a == rhs.id_a)) {
    return lhs.id_a < rhs.id_a;
  }
  return lhs.id_b < rhs.id_b;
}

/// Sort witnesses deterministically and tally them by kind. Mirrors
/// `_build_report` (`overlap.py:354-367`).
OverlapResult build_report(std::vector<Witness> witnesses) {
  std::sort(witnesses.begin(), witnesses.end(), witness_order);
  OverlapResult result;
  result.counts_by_kind = {0, 0, 0};
  for (const Witness& witness : witnesses) {
    ++result.counts_by_kind.at(static_cast<std::size_t>(witness.kind));
  }
  result.witnesses = std::move(witnesses);
  return result;
}

/// Pick a spatial-hash cell size sized to a typical primitive footprint.
/// Mirrors `_hash_cell_size` (`overlap.py:370-385`).
double hash_cell_size(const OverlapParams& params) {
  const double max_half_width = std::max(params.backbone_half_width, params.pair_half_width);
  const double size = 2.0 * (params.node_r + max_half_width);
  return size > 0.0 ? size : 1.0;
}

}  // namespace

OverlapResult check_overlaps(const std::vector<double>& x, const std::vector<double>& y,
                             const std::vector<int>& pair_map, const OverlapParams& params) {
  // overlap.py:409-437.
  validate_inputs(x, y, pair_map);
  const std::vector<Primitive> primitives = build_primitives(x, y, pair_map, params);

  SpatialHash grid(hash_cell_size(params));
  for (std::size_t index = 0; index < primitives.size(); ++index) {
    insert_primitive(grid, static_cast<int>(index), primitives[index]);
  }

  const std::vector<Witness> witnesses =
      collect_witnesses(primitives, grid.candidate_pairs(), pair_map, params.tol);
  return build_report(witnesses);
}

OverlapResult check_overlaps_bruteforce(const std::vector<double>& x, const std::vector<double>& y,
                                        const std::vector<int>& pair_map,
                                        const OverlapParams& params) {
  // overlap.py:440-468: same logic as check_overlaps, all O(n^2) pairs.
  validate_inputs(x, y, pair_map);
  const std::vector<Primitive> primitives = build_primitives(x, y, pair_map, params);

  std::vector<std::pair<int, int>> all_pairs;
  const auto n = primitives.size();
  all_pairs.reserve(n * (n - 1) / 2);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j) {
      all_pairs.emplace_back(static_cast<int>(i), static_cast<int>(j));
    }
  }

  const std::vector<Witness> witnesses =
      collect_witnesses(primitives, all_pairs, pair_map, params.tol);
  return build_report(witnesses);
}

}  // namespace rna_layout::overlap
