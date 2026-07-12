cmake_minimum_required(VERSION 3.15...3.27)
# A fixed project name (not `${SKBUILD_PROJECT_NAME}`) so this also
# configures cleanly for a standalone `cmake -B build && cmake --build
# build && ctest --test-dir build` invocation, not only through
# scikit-build-core's `pip install -e .`.
project(rna_draw LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

# pybind11 is a build-system requirement (pyproject.toml [build-system]);
# locate its CMake config via the interpreter that is actually building us,
# rather than relying on scikit-build-core to pre-populate CMAKE_PREFIX_PATH.
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
  OUTPUT_VARIABLE pybind11_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
find_package(pybind11 CONFIG REQUIRED)

# ViennaRNA HEADERS live in the active conda/micromamba environment --
# confirmed present at $CONDA_PREFIX/include (see plan Step 0). Only the
# headers are needed now: the layout core is vendored and compiled
# standalone (see below), so libRNA.a is NO LONGER a link input and there is
# no ViennaRNA runtime dependency.
if(NOT DEFINED ENV{CONDA_PREFIX})
  message(FATAL_ERROR "CONDA_PREFIX is not set -- run the build via 'micromamba run -n py3 ...'")
endif()
set(VIENNA_PREFIX "$ENV{CONDA_PREFIX}")

# Vendored, STANDALONE copy of RNApuzzler/RNAturtle -- see
# src/vienna_layout/vendor/README.md. These two translation units define
# every layout entry point this extension calls
# (vrna_plot_coords_puzzler[_pt], vrna_plot_options_puzzler[_free],
# vrna_plot_coords_turtle[_pt]). They reference EXACTLY TWO ViennaRNA
# library symbols -- vrna_alloc and vrna_ptable -- which are supplied by the
# ~120 LOC compat shim vrna_compat.c (faithful copies of the upstream 2.7.0
# sources; see that file's header). With the shim in the target, the whole
# layout core links WITHOUT libRNA.a: rna_draw has no ViennaRNA runtime dep.
set(VENDOR_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/vienna_layout/vendor")
set(VENDOR_PUZZLER_DIR "${VENDOR_DIR}/RNApuzzler")
set(VENDOR_PUZZLER_SOURCES
  "${VENDOR_PUZZLER_DIR}/RNApuzzler.c"
  "${VENDOR_PUZZLER_DIR}/RNAturtle.c"
  "${VENDOR_DIR}/vrna_compat.c"
)

pybind11_add_module(_vienna_layout
  src/vienna_layout/bindings.cpp
  ${VENDOR_PUZZLER_SOURCES}
)

target_include_directories(_vienna_layout PRIVATE
  "${VIENNA_PREFIX}/include"
  "${VIENNA_PREFIX}/include/ViennaRNA"
  "${VENDOR_PUZZLER_DIR}"
)

# Vendored third-party C: relaxed warnings, and excluded from this
# project's clang-format/clang-tidy passes (see vendor/README.md) -- do
# not reformat it to this repo's C++ style.
set_source_files_properties(${VENDOR_PUZZLER_SOURCES} PROPERTIES
  LANGUAGE C
  COMPILE_OPTIONS "-w"
)

# The archive was built without LTO (see RNAlib2.pc: `-fno-lto -Wl,-fno-lto`);
# matching it avoids an LTO-mismatch link failure.
target_compile_options(_vienna_layout PRIVATE -fno-lto)
target_link_options(_vienna_layout PRIVATE -fno-lto)

# No libRNA.a: the vendored layout core + vrna_compat.c are self-contained,
# needing only libc/libm/libstdc++.
target_link_libraries(_vienna_layout PRIVATE
  stdc++
  m
)

install(TARGETS _vienna_layout DESTINATION rna_draw)

# `ctest` smoke test (no Python/pybind11 in the loop): proves the vendored
# layout core + vrna_compat.c shim link and run STANDALONE, with libRNA.a
# UNLINKED. The full behavior contract (parity, malformed input, memory,
# ABI) lives in `tests/test_vienna_binding.py`.
enable_testing()
add_executable(test_puzzler_smoke
  src/vienna_layout/tests/test_puzzler_smoke.cpp
  ${VENDOR_PUZZLER_SOURCES}
)
target_include_directories(test_puzzler_smoke PRIVATE
  "${VIENNA_PREFIX}/include"
  "${VIENNA_PREFIX}/include/ViennaRNA"
  "${VENDOR_PUZZLER_DIR}"
)
target_compile_options(test_puzzler_smoke PRIVATE -fno-lto)
target_link_options(test_puzzler_smoke PRIVATE -fno-lto)
target_link_libraries(test_puzzler_smoke PRIVATE
  stdc++
  m
)
add_test(NAME puzzler_smoke COMMAND test_puzzler_smoke)
