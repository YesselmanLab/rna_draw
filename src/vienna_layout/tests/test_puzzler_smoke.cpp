// Minimal `ctest` C++ smoke test: calls the puzzler C API directly (no
// Python/pybind11 in the loop) and checks the returned coordinate count
// matches the structure length. `tests/test_vienna_binding.py` covers the
// full behavior contract (parity, malformed input, memory, ABI); this
// test only proves the CMake build links `libRNA.a` correctly.

#include <cstdlib>
#include <cstring>
#include <iostream>

// clang-format off
extern "C" {
#include <ViennaRNA/vrna_config.h>
#include <ViennaRNA/plotting/layouts.h>
}
// clang-format on

int main() {
  const char* structure = "((((....))))";
  const size_t expected_length = std::strlen(structure);

  float* x = nullptr;
  float* y = nullptr;
  int n = vrna_plot_coords_puzzler(structure, &x, &y, nullptr, nullptr);

  bool ok = (n > 0) && (static_cast<size_t>(n) == expected_length);
  if (!ok) {
    std::cerr << "test_puzzler_smoke: expected length " << expected_length << ", got " << n << "\n";
  }

  std::free(x);
  std::free(y);
  return ok ? 0 : 1;
}
