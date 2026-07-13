/**
 * @file broad_phase.cpp
 * @brief Implements `intersect_tree.hpp`'s `any_intersection` -- the
 *        broad-phase spatial index behind SPEED lever A1
 *        (`.claude/plans/current-plan-speed.md`).
 *
 * CORRECTNESS ARGUMENT (this whole file's purpose -- see `any_intersection`'s
 * doc comment for the summary): the grid buckets each node's AABB, expanded
 * by `margin()`'s conservative bound, into every cell it touches; two nodes
 * are a CANDIDATE pair iff they share a cell, i.e. iff their expanded AABBs
 * overlap. Since `margin()` dominates `intersect_nodes_bounding_boxes`'s
 * (`intersect_tree.cpp:20`) real per-pair `extra_distance` for every pair
 * this file considers, "expanded AABBs overlap" is implied by "the exact
 * AABB-reject test would have accepted" -- so the candidate set is a
 * superset of every pair the ORIGINAL brute-force scan would ever have
 * tested past its own AABB reject. Every candidate is then run through the
 * exact, unmodified `intersect_node_node`; every non-candidate is therefore
 * provably `none`. The OR-over-all-pairs boolean this file returns is thus
 * bit-identical to `intersect_node_lists(subtree, subtree, ...) ||
 * intersect_node_lists(subtree, ancestor_list, ...)` -- EXACT, not
 * approximate.
 *
 * Cell-key grouping uses a sort over `(cell, node_index)` pairs rather than
 * an `unordered_map<cell, vector<node_index>>` -- avoids one small heap
 * allocation per occupied cell on every call (SPEED lever A3's "kill
 * allocation churn", folded in here since this index is exactly the
 * allocation-heavy hot loop that section calls out); a `thread_local` scratch
 * struct (`BroadPhaseScratch`) reuses its buffers' capacity across the many
 * `check_optimize_intersections` calls within one `optimize_tree` (and is
 * therefore also safe under SPEED lever B1's per-thread parallel batch: each
 * thread gets its own independent scratch, no shared mutable state).
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

#include "rna_layout/config_tree.hpp"
#include "rna_layout/geometry.hpp"
#include "rna_layout/intersect_tree.hpp"

namespace rna_layout {

namespace {

/// Which of `check_optimize_intersections`'s two lists a node came from --
/// decides which candidate pairs are even semantically valid to test (see
/// `is_relevant_pair`).
enum class ListTag : std::uint8_t { subtree, ancestor };

struct GridEntry {
  const TreeNode* node = nullptr;
  ListTag tag = ListTag::subtree;
  Aabb expanded;
};

/// Reused across calls (within one thread) to avoid rebuilding this index's
/// backing storage from scratch every `check_optimize_intersections` call --
/// SPEED lever A3, folded into A1's own hot loop (see file header).
struct BroadPhaseScratch {
  std::vector<GridEntry> entries;
  std::vector<std::pair<std::int64_t, int>> cell_hits;  // (packed cell key, entry index)
};

BroadPhaseScratch& scratch_for_this_thread() {
  static thread_local BroadPhaseScratch scratch;
  return scratch;
}

/// Below this many total (subtree + ancestor_list) nodes, building the grid
/// costs more than the brute-force O(m^2) scan it would replace -- fall
/// through to the original `intersect_node_lists` instead (still EXACT: the
/// unmodified function, unconditionally correct at any size).
constexpr std::size_t kBruteForceThreshold = 32;

/// Packs a grid cell's `(cx, cy)` into one map/sort key: a bijective
/// concatenation of each coordinate's low 32 bits. Grid cells routinely have
/// NEGATIVE coordinates (RNA layouts are not confined to the positive
/// quadrant), so this shifts the UNSIGNED reinterpretation of `cx` -- left-
/// shifting a negative signed integer is undefined behavior (caught by
/// UBSan on the pathological hard-set structures), even though it happens
/// to "work" on every mainstream two's-complement target.
[[nodiscard]] std::int64_t pack_cell(std::int64_t cx, std::int64_t cy) {
  const auto ucx = static_cast<std::uint64_t>(cx);
  const auto ucy = static_cast<std::uint64_t>(cy);
  return static_cast<std::int64_t>((ucx << 32) ^ (ucy & 0xFFFFFFFFULL));
}

/// A safe upper bound on `intersect_nodes_bounding_boxes`'s
/// (`intersect_tree.cpp:20-33`) per-pair `extra_distance`: that function
/// adds `epsilon_recognize(clearance)` plus a convex combination (weight
/// `1/count`, `count` in `{1, 2}`) of the two nodes' `StemBox::bulge_dist`,
/// which is therefore always `<= max(bulge_dist1, bulge_dist2) <=
/// max_bulge_dist` over every node this call considers -- see
/// `any_intersection`'s doc comment for how this bound is actually used.
[[nodiscard]] double compute_margin(const std::vector<const TreeNode*>& subtree,
                                    const std::vector<const TreeNode*>& ancestor_list,
                                    double clearance) {
  double max_bulge_dist = 0.0;
  for (const std::vector<const TreeNode*>* list : {&subtree, &ancestor_list}) {
    for (const TreeNode* node : *list) {
      if (is_exterior(*node)) {
        continue;
      }
      // NOLINTNEXTLINE(bugprone-unchecked-optional-access) -- non-root nodes
      // always have `sbox` set (`config_tree.cpp`'s invariant note).
      max_bulge_dist = std::fmax(max_bulge_dist, node->sbox->bulge_dist);
    }
  }
  return geom::epsilon_recognize(clearance) + max_bulge_dist;
}

/// Fills `scratch.entries` with every non-exterior node of @p subtree
/// (tagged `subtree`) and @p ancestor_list (tagged `ancestor`), each AABB
/// expanded by @p margin on every side. Exterior nodes are skipped here --
/// `has_exterior_pair` (below) handles them separately, exactly (they have
/// no `sbox`/meaningful `Aabb` to index).
void collect_geometric_entries(const std::vector<const TreeNode*>& subtree,
                               const std::vector<const TreeNode*>& ancestor_list, double margin,
                               std::vector<GridEntry>& entries) {
  entries.clear();
  for (const TreeNode* node : subtree) {
    if (!is_exterior(*node)) {
      const Aabb& box = node->aabb;
      entries.push_back(GridEntry{node, ListTag::subtree,
                                  Aabb{{box.min.x - margin, box.min.y - margin},
                                       {box.max.x + margin, box.max.y + margin}}});
    }
  }
  for (const TreeNode* node : ancestor_list) {
    if (!is_exterior(*node)) {
      const Aabb& box = node->aabb;
      entries.push_back(GridEntry{node, ListTag::ancestor,
                                  Aabb{{box.min.x - margin, box.min.y - margin},
                                       {box.max.x + margin, box.max.y + margin}}});
    }
  }
}

/// A cell size that keeps each grid cell's expected OCCUPANCY near one entry
/// (a performance knob only -- ANY positive cell size keeps the candidate
/// set a superset, see the file header): `sqrt(bounding-area / count)`, the
/// standard uniform-grid density heuristic, floored so a degenerate
/// (zero-area) node set can never divide by (or tile into) a near-zero cell
/// size. Deliberately NOT based on individual AABB extents (a loop's own
/// box can be large while its subtree is still densely packed -- that
/// heuristic under-partitions exactly the dense subtrees this index most
/// needs to prune).
[[nodiscard]] double compute_cell_size(const std::vector<GridEntry>& entries) {
  constexpr double kMinCellSize = 1.0;
  if (entries.empty()) {
    return kMinCellSize;
  }
  double min_x = entries.front().expanded.min.x;
  double max_x = entries.front().expanded.max.x;
  double min_y = entries.front().expanded.min.y;
  double max_y = entries.front().expanded.max.y;
  for (const GridEntry& entry : entries) {
    min_x = std::fmin(min_x, entry.expanded.min.x);
    max_x = std::fmax(max_x, entry.expanded.max.x);
    min_y = std::fmin(min_y, entry.expanded.min.y);
    max_y = std::fmax(max_y, entry.expanded.max.y);
  }
  const double area = std::fmax(max_x - min_x, 0.0) * std::fmax(max_y - min_y, 0.0);
  if (area <= 0.0) {
    return kMinCellSize;
  }
  return std::fmax(std::sqrt(area / static_cast<double>(entries.size())), kMinCellSize);
}

/// Fills `scratch.cell_hits` with one `(cell_key, entry_index)` per grid
/// cell each entry's expanded AABB touches, then sorts it by cell key so
/// `for_each_candidate_pair` can scan same-cell runs.
void bucket_entries_into_cells(const std::vector<GridEntry>& entries, double cell_size,
                               std::vector<std::pair<std::int64_t, int>>& cell_hits) {
  cell_hits.clear();
  for (int i = 0; i < static_cast<int>(entries.size()); ++i) {
    const Aabb& box = entries[static_cast<std::size_t>(i)].expanded;
    const auto cx0 = static_cast<std::int64_t>(std::floor(box.min.x / cell_size));
    const auto cx1 = static_cast<std::int64_t>(std::floor(box.max.x / cell_size));
    const auto cy0 = static_cast<std::int64_t>(std::floor(box.min.y / cell_size));
    const auto cy1 = static_cast<std::int64_t>(std::floor(box.max.y / cell_size));
    for (std::int64_t cx = cx0; cx <= cx1; ++cx) {
      for (std::int64_t cy = cy0; cy <= cy1; ++cy) {
        cell_hits.emplace_back(pack_cell(cx, cy), i);
      }
    }
  }
  std::sort(cell_hits.begin(), cell_hits.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });
}

/// Whether a candidate pair tagged `(tag_a, tag_b)` is one
/// `check_optimize_intersections` actually tests: subtree-subtree (both
/// `subtree`) or subtree-ancestor (exactly one of each) -- NEVER
/// ancestor-ancestor (the original two `intersect_node_lists` calls never
/// test that combination, so this file must not introduce it).
[[nodiscard]] bool is_relevant_pair(ListTag tag_a, ListTag tag_b) {
  return tag_a == ListTag::subtree || tag_b == ListTag::subtree;
}

/// Scans same-cell runs of @p cell_hits for relevant candidate pairs,
/// running the exact `intersect_node_node` on each and returning `true` on
/// the first hit (early exit, matching `intersect_node_lists`'s own
/// early-return `||` chain).
///
/// NOT deduplicated across cells: a pair whose expanded AABBs span several
/// shared cells is tested once per shared cell. `intersect_node_node` is a
/// pure, side-effect-free predicate, so retesting only costs a little time,
/// never correctness -- and, with `compute_cell_size`'s density-tuned cell
/// size, a pair sharing more than one or two cells is rare in practice, so a
/// dedup set (itself an O(k) lookup/insert per candidate, `k` = candidates
/// found so far) is not worth its own overhead here.
[[nodiscard]] bool query_candidate_pairs(const std::vector<GridEntry>& entries,
                                         const std::vector<std::pair<std::int64_t, int>>& cell_hits,
                                         double clearance) {
  std::size_t run_start = 0;
  const std::size_t total = cell_hits.size();
  while (run_start < total) {
    std::size_t run_end = run_start + 1;
    while (run_end < total && cell_hits[run_end].first == cell_hits[run_start].first) {
      ++run_end;
    }
    for (std::size_t a = run_start; a < run_end; ++a) {
      for (std::size_t b = a + 1; b < run_end; ++b) {
        const GridEntry& entry_a = entries[static_cast<std::size_t>(cell_hits[a].second)];
        const GridEntry& entry_b = entries[static_cast<std::size_t>(cell_hits[b].second)];
        if (!is_relevant_pair(entry_a.tag, entry_b.tag)) {
          continue;
        }
        if (intersect_node_node(*entry_a.node, *entry_b.node, clearance).type !=
            IntersectionType::none) {
          return true;
        }
      }
    }
    run_start = run_end;
  }
  return false;
}

/// The exterior-root half of `intersect_node_lists`'s semantics
/// (`intersect_tree.cpp:227-247`), reproduced exactly: at most one node in
/// @p list1 ∪ @p list2 is ever the true exterior root in practice (only
/// `ancestor_list`'s last entry can be), so this is O(|list1| + |list2|),
/// not the O(m^2) this file exists to avoid for the geometric case.
[[nodiscard]] bool has_exterior_pair(const std::vector<const TreeNode*>& list1,
                                     const std::vector<const TreeNode*>& list2,
                                     bool check_exterior_intersections, double clearance) {
  for (const TreeNode* node1 : list1) {
    if (is_exterior(*node1)) {
      for (const TreeNode* node2 : list2) {
        if (intersect_node_exterior(*node2, check_exterior_intersections, clearance)) {
          return true;
        }
      }
    }
  }
  for (const TreeNode* node2 : list2) {
    if (is_exterior(*node2)) {
      for (const TreeNode* node1 : list1) {
        if (intersect_node_exterior(*node1, check_exterior_intersections, clearance)) {
          return true;
        }
      }
    }
  }
  return false;
}

/// The geometric (non-exterior x non-exterior) half of `any_intersection`:
/// builds the broad-phase grid over @p subtree ∪ @p ancestor_list and
/// queries it. Exterior pairs are @p has_exterior_pair's job, called
/// separately by `any_intersection` before this.
[[nodiscard]] bool broad_phase_geometric_intersection(
    const std::vector<const TreeNode*>& subtree, const std::vector<const TreeNode*>& ancestor_list,
    double clearance) {
  BroadPhaseScratch& scratch = scratch_for_this_thread();
  const double margin = compute_margin(subtree, ancestor_list, clearance);
  collect_geometric_entries(subtree, ancestor_list, margin, scratch.entries);
  if (scratch.entries.size() < 2) {
    return false;
  }
  const double cell_size = compute_cell_size(scratch.entries);
  bucket_entries_into_cells(scratch.entries, cell_size, scratch.cell_hits);
  return query_candidate_pairs(scratch.entries, scratch.cell_hits, clearance);
}

}  // namespace

bool any_intersection(const std::vector<const TreeNode*>& subtree,
                      const std::vector<const TreeNode*>& ancestor_list,
                      bool check_exterior_intersections, double clearance) {
  if (subtree.size() + ancestor_list.size() <= kBruteForceThreshold) {
    return intersect_node_lists(subtree, subtree, check_exterior_intersections, clearance) ||
           intersect_node_lists(subtree, ancestor_list, check_exterior_intersections, clearance);
  }
  if (has_exterior_pair(subtree, subtree, check_exterior_intersections, clearance) ||
      has_exterior_pair(subtree, ancestor_list, check_exterior_intersections, clearance)) {
    return true;
  }
  return broad_phase_geometric_intersection(subtree, ancestor_list, clearance);
}

}  // namespace rna_layout
