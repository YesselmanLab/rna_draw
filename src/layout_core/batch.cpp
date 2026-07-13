/**
 * @file batch.cpp
 * @brief Implements `batch.hpp`'s `layout_puzzler_batch`.
 */

#include "rna_layout/batch.hpp"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <thread>
#include <vector>

#include "rna_layout/puzzler.hpp"

namespace rna_layout {

namespace {

/// `num_threads <= 0` -> `hardware_concurrency()` (falling back to `1` if
/// unreported), capped at the number of items -- more threads than work is
/// pure overhead.
[[nodiscard]] unsigned resolve_thread_count(int num_threads, std::size_t work_items) {
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

/// Lays out `structures[i]` into `results[i]` for whichever indices this
/// call claims from @p next_index (an atomic work-queue shared by every
/// worker thread) -- each index is claimed by exactly one thread, so
/// `results`' element writes never race (distinct indices, no shared
/// mutable state elsewhere -- see `batch.hpp`'s file header).
void run_worker(const std::vector<std::string>& structures, const PuzzlerOptions& opts,
                std::atomic<std::size_t>& next_index, std::vector<BatchResult>& results) {
  const std::size_t total = structures.size();
  std::size_t index;
  while ((index = next_index.fetch_add(1, std::memory_order_relaxed)) < total) {
    BatchResult& out = results[index];
    try {
      out.coords = layout_puzzler(structures[index], opts);
      out.ok = true;
    } catch (const std::exception& ex) {
      out.ok = false;
      out.error = ex.what();
    } catch (...) {
      // Never let one element's failure escape the batch -- see this
      // function's doc comment (`batch.hpp`).
      out.ok = false;
      out.error = "unknown error";
    }
  }
}

}  // namespace

std::vector<BatchResult> layout_puzzler_batch(const std::vector<std::string>& structures,
                                              const PuzzlerOptions& opts, int num_threads) {
  std::vector<BatchResult> results(structures.size());
  if (structures.empty()) {
    return results;
  }

  const unsigned thread_count = resolve_thread_count(num_threads, structures.size());
  std::atomic<std::size_t> next_index{0};

  if (thread_count <= 1) {
    run_worker(structures, opts, next_index, results);
    return results;
  }

  std::vector<std::thread> workers;
  workers.reserve(thread_count);
  for (unsigned t = 0; t < thread_count; ++t) {
    workers.emplace_back([&structures, &opts, &next_index, &results]() {
      run_worker(structures, opts, next_index, results);
    });
  }
  for (std::thread& worker : workers) {
    worker.join();
  }
  return results;
}

}  // namespace rna_layout
