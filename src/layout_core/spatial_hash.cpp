/**
 * @file spatial_hash.cpp
 * @brief Implements `spatial_hash.hpp` -- port of `rna_draw/spatial_hash.py:33-113`.
 */

#include "rna_layout/spatial_hash.hpp"

#include <algorithm>
#include <cmath>
#include <set>

namespace rna_layout::overlap {

std::size_t SpatialHash::CellKeyHash::operator()(const CellKey& key) const noexcept {
  // A standard 64-bit mix (boost::hash_combine-style); cell coordinates
  // are a performance detail, so any well-distributed hash is fine.
  std::size_t seed = std::hash<CellCoord>{}(key.first);
  seed ^= std::hash<CellCoord>{}(key.second) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
  return seed;
}

SpatialHash::SpatialHash(double cell_size) : cell_size_(cell_size) {}

void SpatialHash::insert(int item_index, const Aabb& aabb) {
  // spatial_hash.py:45-58.
  const auto fx0 = static_cast<CellCoord>(std::floor(aabb.min.x / cell_size_));
  const auto fx1 = static_cast<CellCoord>(std::floor(aabb.max.x / cell_size_));
  const auto fy0 = static_cast<CellCoord>(std::floor(aabb.min.y / cell_size_));
  const auto fy1 = static_cast<CellCoord>(std::floor(aabb.max.y / cell_size_));
  for (CellCoord cx = fx0; cx <= fx1; ++cx) {
    for (CellCoord cy = fy0; cy <= fy1; ++cy) {
      cells_[{cx, cy}].push_back(item_index);
    }
  }
}

void SpatialHash::insert_segment(int item_index, Vec2 p0, Vec2 p1, double pad) {
  // spatial_hash.py:60-97.
  const double length = std::hypot(p1.x - p0.x, p1.y - p0.y);
  if (length == 0.0) {
    insert(item_index, Aabb{Vec2{p0.x - pad, p0.y - pad}, Vec2{p0.x + pad, p0.y + pad}});
    return;
  }
  const long steps = std::max(1L, static_cast<long>(std::ceil(length / cell_size_)));
  const double dx = (p1.x - p0.x) / static_cast<double>(steps);
  const double dy = (p1.y - p0.y) / static_cast<double>(steps);
  for (long k = 0; k < steps; ++k) {
    const double px0 = p0.x + static_cast<double>(k) * dx;
    const double py0 = p0.y + static_cast<double>(k) * dy;
    const double px1 = px0 + dx;
    const double py1 = py0 + dy;
    const Aabb aabb{Vec2{std::min(px0, px1) - pad, std::min(py0, py1) - pad},
                    Vec2{std::max(px0, px1) + pad, std::max(py0, py1) + pad}};
    insert(item_index, aabb);
  }
}

std::vector<std::pair<int, int>> SpatialHash::candidate_pairs() const {
  // spatial_hash.py:99-113: every (i, j), i < j, sharing a cell, deduped
  // (Python uses a `set`) -- a `std::set` gives both dedup and a stable
  // ascending order for free.
  std::set<std::pair<int, int>> pairs;
  for (const auto& [key, items] : cells_) {
    (void)key;
    const std::size_t n = items.size();
    for (std::size_t a = 0; a < n; ++a) {
      for (std::size_t b = a + 1; b < n; ++b) {
        const int i = items[a];
        const int j = items[b];
        pairs.emplace(i < j ? std::make_pair(i, j) : std::make_pair(j, i));
      }
    }
  }
  return {pairs.begin(), pairs.end()};
}

}  // namespace rna_layout::overlap
