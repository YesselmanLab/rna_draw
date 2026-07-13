#pragma once

/**
 * @file spatial_hash.hpp
 * @brief Uniform-grid spatial hash for candidate overlap-pair generation.
 *        Port of `rna_draw/spatial_hash.py:33-113`.
 *
 * Correctness argument (verbatim from the Python module docstring): two
 * footprints overlap implies their (already inflated) axis-aligned
 * bounding boxes overlap, which implies some grid cell is covered by both
 * boxes, which implies both items were inserted into that cell, which
 * implies the pair is emitted by `candidate_pairs`. So `candidate_pairs`
 * is a SUPERSET of every true overlap for ANY positive `cell_size` -- the
 * cell size can only ever affect speed, never correctness (`.claude/plans/
 * current-plan-checker.md`'s Risk R2). `candidate_pairs` DEDUPES (the
 * Python side returns a `set`) -- the checker must not double-count a
 * witness.
 */

#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

#include "rna_layout/types.hpp"

namespace rna_layout::overlap {

/// Grid-cell coordinate. `std::int64_t` (not `int`): production coordinates
/// are bounded well within range, but a cell index is `floor(coord /
/// cell_size)`, which the wider type guards against overflow for even a
/// generously out-of-range input.
using CellCoord = std::int64_t;

/**
 * Buckets items by the grid cells their AABBs touch. Mirrors
 * `SpatialHash` (`spatial_hash.py:33-113`); `cell_size` is a PERFORMANCE
 * KNOB ONLY (see this file's header) -- it never affects correctness.
 */
class SpatialHash {
 public:
  /**
   * @param cell_size Side length of a square grid cell; must be positive
   *     (see `overlap_check.hpp`'s `hash_cell_size`, which guards this).
   */
  explicit SpatialHash(double cell_size);

  /**
   * Insert an item into every grid cell @p aabb touches. Mirrors `insert`
   * (`spatial_hash.py:45-58`).
   *
   * @param item_index Identifier for the item (its primitive-list index).
   * @param aabb Bounding box, `min`/`max` corners.
   */
  void insert(int item_index, const Aabb& aabb);

  /**
   * Insert an item as a segment, tiled into @p pad-padded pieces so a
   * long diagonal chord doesn't blow up to one huge AABB. Mirrors
   * `insert_segment` (`spatial_hash.py:60-97`): `length == 0.0` falls
   * back to a single padded point AABB; otherwise the segment is split
   * into `max(1, ceil(length / cell_size))` equal pieces.
   *
   * @param item_index Identifier for the item.
   * @param p0 Segment start point.
   * @param p1 Segment end point.
   * @param pad Inflation applied to each piece's bounding box.
   */
  void insert_segment(int item_index, Vec2 p0, Vec2 p1, double pad);

  /**
   * Every `(i, j)` with `i < j` sharing at least one grid cell, deduped
   * and returned in ascending order. Mirrors `candidate_pairs`
   * (`spatial_hash.py:99-113`), which returns a Python `set` -- dedup is
   * load-bearing, not an optimization (a duplicate candidate pair would
   * double-count a witness).
   */
  [[nodiscard]] std::vector<std::pair<int, int>> candidate_pairs() const;

 private:
  using CellKey = std::pair<CellCoord, CellCoord>;

  /// Hash for `CellKey`, since `std::pair` has no default `std::hash`.
  struct CellKeyHash {
    [[nodiscard]] std::size_t operator()(const CellKey& key) const noexcept;
  };

  double cell_size_;
  std::unordered_map<CellKey, std::vector<int>, CellKeyHash> cells_;
};

}  // namespace rna_layout::overlap
