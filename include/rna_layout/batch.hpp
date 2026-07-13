#pragma once

/**
 * @file batch.hpp
 * @brief Parallel batch entry point for `layout_puzzler` -- SPEED lever B1
 *        (`.claude/plans/current-plan-speed.md`).
 *
 * `layout_puzzler` is a pure, deterministic function of its arguments: no
 * mutable global/static/thread_local state anywhere in `src/layout_core/`
 * or `include/rna_layout/` is shared ACROSS a call (the one `thread_local`,
 * `broad_phase.cpp`'s per-thread scratch, exists PRECISELY so each thread
 * gets its own independent copy -- see that file's header). `layout_puzzler`
 * is therefore trivially safe to fan out over a `std::thread` pool with no
 * locking. This file is the fan-out; `bindings.cpp`'s
 * `plot_coords_puzzler_batch` is the Python-facing entry point, wrapped in
 * `py::gil_scoped_release` so the threads actually run concurrently.
 */

#include <string>
#include <vector>

#include "rna_layout/turtle.hpp"
#include "rna_layout/types.hpp"

namespace rna_layout {

/// One `layout_puzzler_batch` element's outcome: either @ref coords (on
/// success) or @ref error (the exception message a malformed/degenerate
/// structure raised) -- PER-ELEMENT failure isolation, so one bad structure
/// never aborts the rest of the batch.
struct BatchResult {
  Coords coords;
  bool ok = false;
  std::string error;
};

/**
 * Lay out every one of @p structures with `layout_puzzler`, fanned out over
 * @p num_threads worker threads (`std::thread`, an atomic work-index over
 * @p structures -- no locking needed, see this file's header).
 *
 * @param structures Dot-bracket secondary structures, in any order.
 * @param opts Shared layout options, applied to every structure.
 * @param num_threads Worker thread count; `<= 0` uses
 *     `std::thread::hardware_concurrency()` (falling back to `1` if that
 *     reports `0`), capped at `structures.size()`.
 * @return One `BatchResult` per @p structures element, SAME ORDER, each
 *     independent of every other element and of @p num_threads (a
 *     malformed element's `std::invalid_argument` is caught and recorded in
 *     that element's `error`, never thrown out of this function).
 */
[[nodiscard]] std::vector<BatchResult> layout_puzzler_batch(
    const std::vector<std::string>& structures, const PuzzlerOptions& opts, int num_threads);

}  // namespace rna_layout
