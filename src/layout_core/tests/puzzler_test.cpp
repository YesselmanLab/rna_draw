// ctest: `rna_layout::layout_puzzler` (`puzzler.hpp`), Milestone A steps 6-7
// ("finalization with resolver OFF" + the SIBLING intersection resolver).
// Parity against the vendored oracle is a Python-level concern
// (`tests/test_native_parity.py`); this proves the native pipeline runs,
// guards its ancestor/optimize preconditions, and produces sane,
// deterministic output with `check_sibling` on or off.

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

void test_default_options_throw() {
  // Defaults are ancestor-ON and optimize-ON (both not yet ported).
  bool threw = false;
  try {
    const Coords unused = layout_puzzler(std::string("((((....))))"), PuzzlerOptions{});
    (void)unused;
  } catch (const std::logic_error&) {
    threw = true;
  }
  expect(threw, "layout_puzzler: the default PuzzlerOptions (ancestor+optimize ON) throw");
}

void test_ancestor_or_optimize_alone_throws_but_sibling_alone_does_not() {
  for (bool sibling : {true, false}) {
    for (bool ancestor : {true, false}) {
      for (bool optimize : {true, false}) {
        PuzzlerOptions opts = resolver_off_options();
        opts.check_sibling = sibling;
        opts.check_ancestor = ancestor;
        opts.optimize = optimize;

        bool threw = false;
        try {
          const Coords unused = layout_puzzler(std::string("((((....))))"), opts);
          (void)unused;
        } catch (const std::logic_error&) {
          threw = true;
        }
        const bool should_throw = ancestor || optimize;
        expect(threw == should_throw,
               "layout_puzzler: throws iff check_ancestor or optimize is requested "
               "(check_sibling alone is Milestone A step 7, implemented)");
      }
    }
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
  test_default_options_throw();
  test_ancestor_or_optimize_alone_throws_but_sibling_alone_does_not();
  test_resolver_off_produces_sane_coordinates();
  test_resolver_off_is_deterministic();
  test_sibling_resolver_produces_sane_coordinates();
  test_sibling_resolver_is_deterministic();
  test_malformed_input_throws();

  if (g_failures > 0) {
    std::cerr << "puzzler_test: " << g_failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
