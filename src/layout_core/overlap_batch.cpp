/**
 * @file overlap_batch.cpp
 * @brief Implements `overlap_batch.hpp`'s `check_overlaps_batch`. Mirrors
 *        `batch.cpp`'s atomic-work-index fan-out shape.
 */

#include "rna_layout/overlap_batch.hpp"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <stdexcept>
#include <thread>

#include "rna_layout/overlap_check.hpp"

namespace rna_layout::overlap {

namespace {

/// `num_threads <= 0` -> `hardware_concurrency()` (falling back to `1` if
/// unreported), capped at the number of items. Identical policy to
/// `batch.cpp`'s `resolve_thread_count`.
unsigned resolve_thread_count(int num_threads, std::size_t work_items) {
  unsigned requested;
  if (num_threads > 0) {
    requested = static_cast<unsigned>(num_threads);
  } else {
    requested = std::thread::hardware_concurrency();
    if (requested == 0) {
      requested = 1;
    }
  }
  const auto capped = std::min<std::size_t>(requested, std::max<std::size_t>(work_items, 1));
  return static_cast<unsigned>(capped);
}

/// Checks `xs[i]/ys[i]/pair_maps[i]` into `counts[i]` for whichever
/// indices this call claims from @p next_index -- each index is claimed
/// by exactly one thread, so `counts`' element writes never race.
void run_worker(const std::vector<std::vector<double>>& xs,
                const std::vector<std::vector<double>>& ys,
                const std::vector<std::vector<int>>& pair_maps, const std::vector<double>& node_rs,
                double half_width_factor, double tol, std::atomic<std::size_t>& next_index,
                std::vector<int>& counts) {
  const std::size_t total = xs.size();
  std::size_t index;
  while ((index = next_index.fetch_add(1, std::memory_order_relaxed)) < total) {
    try {
      const double node_r = node_rs[index];
      const OverlapParams params{node_r, half_width_factor * node_r, half_width_factor * node_r,
                                 tol};
      counts[index] = check_overlaps(xs[index], ys[index], pair_maps[index], params).num_overlaps();
    } catch (const std::exception&) {
      // Never let one element's failure escape the batch -- see
      // overlap_batch.hpp's file header.
      counts[index] = kBatchFailureSentinel;
    } catch (...) {
      counts[index] = kBatchFailureSentinel;
    }
  }
}

/// Guards the outer per-structure lists' lengths before any worker thread
/// touches them: `run_worker` indexes `ys[index]`/`pair_maps[index]`/
/// `node_rs[index]` with unchecked `operator[]` (the per-element try/catch
/// only covers exceptions `check_overlaps` itself throws), so a mismatched
/// outer length would otherwise be an out-of-bounds read, not a catchable
/// error. Mirrors `validate_inputs`'s (`overlap_check.cpp`) raise contract:
/// `std::invalid_argument`, which pybind11 maps to Python `ValueError`.
void validate_batch_lengths(const std::vector<std::vector<double>>& xs,
                            const std::vector<std::vector<double>>& ys,
                            const std::vector<std::vector<int>>& pair_maps,
                            const std::vector<double>& node_rs) {
  if (!(xs.size() == ys.size() && ys.size() == pair_maps.size() &&
        pair_maps.size() == node_rs.size())) {
    throw std::invalid_argument("xs, ys, pair_maps, and node_rs must all have equal length");
  }
}

}  // namespace

std::vector<int> check_overlaps_batch(const std::vector<std::vector<double>>& xs,
                                      const std::vector<std::vector<double>>& ys,
                                      const std::vector<std::vector<int>>& pair_maps,
                                      const std::vector<double>& node_rs, double half_width_factor,
                                      double tol, int num_threads) {
  validate_batch_lengths(xs, ys, pair_maps, node_rs);
  std::vector<int> counts(xs.size(), kBatchFailureSentinel);
  if (xs.empty()) {
    return counts;
  }

  const unsigned thread_count = resolve_thread_count(num_threads, xs.size());
  std::atomic<std::size_t> next_index{0};

  if (thread_count <= 1) {
    run_worker(xs, ys, pair_maps, node_rs, half_width_factor, tol, next_index, counts);
    return counts;
  }

  std::vector<std::thread> workers;
  workers.reserve(thread_count);
  for (unsigned t = 0; t < thread_count; ++t) {
    workers.emplace_back([&]() {
      run_worker(xs, ys, pair_maps, node_rs, half_width_factor, tol, next_index, counts);
    });
  }
  for (std::thread& worker : workers) {
    worker.join();
  }
  return counts;
}

}  // namespace rna_layout::overlap
