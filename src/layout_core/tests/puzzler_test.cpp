// ctest: `rna_layout::layout_puzzler` (`puzzler.hpp`), Milestone A steps 6-8
// ("finalization with resolver OFF" + the SIBLING and ANCESTOR intersection
// resolvers). Parity against the vendored oracle is a Python-level concern
// (`tests/test_native_parity.py`); this proves the native pipeline runs,
// guards its optimize precondition, and produces sane, deterministic output
// with `check_sibling`/`check_ancestor` on or off (independently or
// together).

#include "rna_layout/puzzler.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

using rna_layout::Coords;
using rna_layout::layout_puzzler;
using rna_layout::PuzzlerOptions;

namespace {

int g_failures = 0;

void expect(bool condition, const char* what) {
  if (!condition) {
    std::cerr << "FAIL " << what << "\n";
    ++g_failures;
  }
}

PuzzlerOptions resolver_off_options() {
  PuzzlerOptions opts;
  opts.check_sibling = false;
  opts.check_ancestor = false;
  opts.optimize = false;
  // `check_exterior` deliberately left at its default (`true`): see
  // `puzzler.cpp`'s `run_config_tree_pipeline` doc comment -- this keeps
  // `update_bounding_boxes` in the pipeline (matching the vendored default
  // options), while `checkAndFixIntersections` is still provably a no-op.
  return opts;
}

void test_default_options_produce_sane_coordinates() {
  // Milestone A step 9: the defaults (sibling+ancestor+optimize all ON,
  // matching the vendored `vrna_plot_options_puzzler()` -- see
  // `puzzler.hpp`'s file header) are now the FULL, real pipeline.
  const std::string structure = "((((...)))((....(((.....))).))(((...))))";
  const Coords coords = layout_puzzler(structure, PuzzlerOptions{});

  expect(coords.x.size() == structure.size(),
         "layout_puzzler (default options): output length matches input length");
  for (std::size_t i = 0; i < coords.x.size(); ++i) {
    expect(std::isfinite(coords.x[i]),
           "layout_puzzler (default options): every x coordinate is finite");
    expect(std::isfinite(coords.y[i]),
           "layout_puzzler (default options): every y coordinate is finite");
  }
}

void test_every_resolver_combination_runs_without_throwing() {
  // Milestone A step 9: every combination of check_sibling/check_ancestor/
  // optimize is now real (none throws) -- unlike pre-step-9, where
  // optimize == true always threw.
  for (bool sibling : {true, false}) {
    for (bool ancestor : {true, false}) {
      for (bool optimize : {true, false}) {
        PuzzlerOptions opts = resolver_off_options();
        opts.check_sibling = sibling;
        opts.check_ancestor = ancestor;
        opts.optimize = optimize;

        const Coords coords = layout_puzzler(std::string("((((....))))"), opts);
        expect(coords.x.size() == 12,
               "layout_puzzler: every check_sibling/check_ancestor/optimize combination "
               "produces output of the expected length without throwing");
      }
    }
  }
}

void test_sibling_and_ancestor_resolver_produces_sane_coordinates() {
  const std::string structure = "((((...)))((....(((.....))).))(((...))))";
  PuzzlerOptions opts = resolver_off_options();
  opts.check_sibling = true;
  opts.check_ancestor = true;
  const Coords coords = layout_puzzler(structure, opts);

  expect(coords.x.size() == structure.size(),
         "layout_puzzler (sibling+ancestor on): output length matches input length");
  for (std::size_t i = 0; i < coords.x.size(); ++i) {
    expect(std::isfinite(coords.x[i]),
           "layout_puzzler (sibling+ancestor on): every x coordinate is finite");
    expect(std::isfinite(coords.y[i]),
           "layout_puzzler (sibling+ancestor on): every y coordinate is finite");
  }
}

void test_sibling_and_ancestor_resolver_is_deterministic() {
  const std::string structure = "((((..((((....))))..))))(((...)))(((...)))(((...)))";
  PuzzlerOptions opts = resolver_off_options();
  opts.check_sibling = true;
  opts.check_ancestor = true;
  const Coords first = layout_puzzler(structure, opts);
  const Coords second = layout_puzzler(structure, opts);
  expect(first.x == second.x,
         "layout_puzzler (sibling+ancestor on): repeated calls produce identical x coordinates");
  expect(first.y == second.y,
         "layout_puzzler (sibling+ancestor on): repeated calls produce identical y coordinates");
}

void test_full_pipeline_produces_sane_coordinates() {
  // All three resolver stages on -- the FULL production config (Milestone A
  // step 9): `resolver_off_options()` + explicit overrides, rather than
  // `PuzzlerOptions{}`, so this stays independent of that struct's own
  // defaults ever changing.
  const std::string structure = "((((...)))((....(((.....))).))(((...))))";
  PuzzlerOptions opts = resolver_off_options();
  opts.check_sibling = true;
  opts.check_ancestor = true;
  opts.optimize = true;
  const Coords coords = layout_puzzler(structure, opts);

  expect(coords.x.size() == structure.size(),
         "layout_puzzler (full pipeline): output length matches input length");
  for (std::size_t i = 0; i < coords.x.size(); ++i) {
    expect(std::isfinite(coords.x[i]),
           "layout_puzzler (full pipeline): every x coordinate is finite");
    expect(std::isfinite(coords.y[i]),
           "layout_puzzler (full pipeline): every y coordinate is finite");
  }
}

void test_full_pipeline_is_deterministic() {
  const std::string structure = "((((..((((....))))..))))(((...)))(((...)))(((...)))";
  PuzzlerOptions opts = resolver_off_options();
  opts.check_sibling = true;
  opts.check_ancestor = true;
  opts.optimize = true;
  const Coords first = layout_puzzler(structure, opts);
  const Coords second = layout_puzzler(structure, opts);
  expect(first.x == second.x,
         "layout_puzzler (full pipeline): repeated calls produce identical x coordinates");
  expect(first.y == second.y,
         "layout_puzzler (full pipeline): repeated calls produce identical y coordinates");
}

void test_ancestor_only_resolver_produces_sane_coordinates() {
  const std::string structure = "((((...)))((....(((.....))).))(((...))))";
  PuzzlerOptions opts = resolver_off_options();
  opts.check_sibling = false;
  opts.check_ancestor = true;
  const Coords coords = layout_puzzler(structure, opts);

  expect(coords.x.size() == structure.size(),
         "layout_puzzler (ancestor on): output length matches input length");
  for (std::size_t i = 0; i < coords.x.size(); ++i) {
    expect(std::isfinite(coords.x[i]),
           "layout_puzzler (ancestor on): every x coordinate is finite");
    expect(std::isfinite(coords.y[i]),
           "layout_puzzler (ancestor on): every y coordinate is finite");
  }
}

void test_sibling_resolver_produces_sane_coordinates() {
  const std::string structure = "((((...)))(((...)))(((...))))";
  PuzzlerOptions opts = resolver_off_options();
  opts.check_sibling = true;
  const Coords coords = layout_puzzler(structure, opts);

  expect(coords.x.size() == structure.size(),
         "layout_puzzler (sibling on): output length matches input length");
  for (std::size_t i = 0; i < coords.x.size(); ++i) {
    expect(std::isfinite(coords.x[i]), "layout_puzzler (sibling on): every x coordinate is finite");
    expect(std::isfinite(coords.y[i]), "layout_puzzler (sibling on): every y coordinate is finite");
  }
}

void test_sibling_resolver_is_deterministic() {
  const std::string structure = "((((..((((....))))..))))(((...)))(((...)))(((...)))";
  PuzzlerOptions opts = resolver_off_options();
  opts.check_sibling = true;
  const Coords first = layout_puzzler(structure, opts);
  const Coords second = layout_puzzler(structure, opts);
  expect(first.x == second.x,
         "layout_puzzler (sibling on): repeated calls produce identical x coordinates");
  expect(first.y == second.y,
         "layout_puzzler (sibling on): repeated calls produce identical y coordinates");
}

void test_resolver_off_produces_sane_coordinates() {
  const std::string structure = "((((...)))(((...)))(((...))))";
  const Coords coords = layout_puzzler(structure, resolver_off_options());

  expect(coords.x.size() == structure.size(), "layout_puzzler: output length matches input length");
  expect(coords.y.size() == structure.size(),
         "layout_puzzler: output length matches input length (y)");

  for (std::size_t i = 0; i < coords.x.size(); ++i) {
    expect(std::isfinite(coords.x[i]), "layout_puzzler: every x coordinate is finite");
    expect(std::isfinite(coords.y[i]), "layout_puzzler: every y coordinate is finite");
  }
}

void test_resolver_off_is_deterministic() {
  const std::string structure = "((((..((((....))))..))))";
  const Coords first = layout_puzzler(structure, resolver_off_options());
  const Coords second = layout_puzzler(structure, resolver_off_options());
  expect(first.x == second.x, "layout_puzzler: repeated calls produce identical x coordinates");
  expect(first.y == second.y, "layout_puzzler: repeated calls produce identical y coordinates");
}

void test_malformed_input_throws() {
  bool threw = false;
  try {
    const Coords unused = layout_puzzler(std::string(""), resolver_off_options());
    (void)unused;
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  expect(threw, "layout_puzzler: empty structure throws invalid_argument");

  threw = false;
  try {
    const Coords unused = layout_puzzler(std::string("()"), resolver_off_options());
    (void)unused;
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  expect(threw, "layout_puzzler: empty-loop structure throws invalid_argument");
}

}  // namespace

int main() {
  test_default_options_produce_sane_coordinates();
  test_every_resolver_combination_runs_without_throwing();
  test_resolver_off_produces_sane_coordinates();
  test_resolver_off_is_deterministic();
  test_sibling_resolver_produces_sane_coordinates();
  test_sibling_resolver_is_deterministic();
  test_ancestor_only_resolver_produces_sane_coordinates();
  test_sibling_and_ancestor_resolver_produces_sane_coordinates();
  test_sibling_and_ancestor_resolver_is_deterministic();
  test_full_pipeline_produces_sane_coordinates();
  test_full_pipeline_is_deterministic();
  test_malformed_input_throws();

  if (g_failures > 0) {
    std::cerr << "puzzler_test: " << g_failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
