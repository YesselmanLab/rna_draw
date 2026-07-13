#pragma once

/**
 * @file overlap_batch.hpp
 * @brief Parallel batch entry point for `check_overlaps` -- the throughput
 *        win for ~1M-structure batch QC (`.claude/plans/
 *        current-plan-checker.md`).
 *
 * `check_overlaps` is a pure, deterministic function of its arguments (no
 * shared mutable state), so it fans out over a `std::thread` pool exactly
 * like `layout_puzzler_batch` (`batch.hpp`) -- see that file's header for
 * the argument. `bindings.cpp`'s `check_overlaps_batch` is the
 * Python-facing entry point, GIL released.
 */

#include <vector>

namespace rna_layout::overlap {

/// One `check_overlaps_batch` element's outcome: the total overlap count
/// (`OverlapResult::num_overlaps()`), or the sentinel `-1` for a
/// per-element failure (e.g. a malformed/asymmetric `pair_map`) -- PER-
/// ELEMENT failure isolation, mirroring `batch.hpp`'s `BatchResult`: one
/// bad structure never aborts the rest of the batch.
inline constexpr int kBatchFailureSentinel = -1;

/**
 * Check every structure in @p xs/@p ys/@p pair_maps for overlaps, fanned
 * out over @p num_threads worker threads. Per-structure `node_r` (matching
 * production, where `node_r` varies per structure); `half_width_factor`
 * reproduces `pipeline_qc.py`'s `backbone_half_width = pair_half_width =
 * 0.75 * node_r` convention uniformly across the batch.
 *
 * @param xs Per-structure nucleotide x-coordinates.
 * @param ys Per-structure nucleotide y-coordinates, same order/length as
 *     @p xs.
 * @param pair_maps Per-structure pair maps, same order/length as @p xs.
 * @param node_rs Per-structure disk radius, same order/length as @p xs.
 * @param half_width_factor Multiplied by each structure's own `node_r` to
 *     get `backbone_half_width`/`pair_half_width`.
 * @param tol Shared tolerance passed to every structure's `OverlapParams`.
 * @param num_threads Worker thread count; `<= 0` uses
 *     `std::thread::hardware_concurrency()` (falling back to `1`),
 *     capped at the number of structures.
 * @return One count per input structure, same order; `kBatchFailureSentinel`
 *     for a structure whose inputs are malformed (never thrown out of this
 *     function).
 */
[[nodiscard]] std::vector<int> check_overlaps_batch(const std::vector<std::vector<double>>& xs,
                                                    const std::vector<std::vector<double>>& ys,
                                                    const std::vector<std::vector<int>>& pair_maps,
                                                    const std::vector<double>& node_rs,
                                                    double half_width_factor, double tol,
                                                    int num_threads);

}  // namespace rna_layout::overlap
