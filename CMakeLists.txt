cmake_minimum_required(VERSION 3.15...3.27)
# A fixed project name (not `${SKBUILD_PROJECT_NAME}`) so this also
# configures cleanly for a standalone `cmake -B build && cmake --build
# build && ctest --test-dir build` invocation, not only through
# scikit-build-core's `pip install -e .`.
project(rna_draw LANGUAGES CXX)

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

# ViennaRNA (RNAlib) lives in the active conda/micromamba environment --
# confirmed present at $CONDA_PREFIX/{include,lib} (see plan Step 0).
if(NOT DEFINED ENV{CONDA_PREFIX})
  message(FATAL_ERROR "CONDA_PREFIX is not set -- run the build via 'micromamba run -n py3 ...'")
endif()
set(VIENNA_PREFIX "$ENV{CONDA_PREFIX}")

pybind11_add_module(_vienna_layout src/vienna_layout/bindings.cpp)

target_include_directories(_vienna_layout PRIVATE
  "${VIENNA_PREFIX}/include"
  "${VIENNA_PREFIX}/include/ViennaRNA"
)

# The archive was built without LTO (see RNAlib2.pc: `-fno-lto -Wl,-fno-lto`);
# matching it avoids an LTO-mismatch link failure.
target_compile_options(_vienna_layout PRIVATE -fno-lto)
target_link_options(_vienna_layout PRIVATE -fno-lto)

find_library(VIENNA_RNA_LIB RNA PATHS "${VIENNA_PREFIX}/lib" REQUIRED NO_DEFAULT_PATH)
find_package(Threads REQUIRED)

target_link_libraries(_vienna_layout PRIVATE
  "${VIENNA_RNA_LIB}"
  Threads::Threads
  stdc++
  m
)

install(TARGETS _vienna_layout DESTINATION rna_draw)

# `ctest` smoke test (no Python/pybind11 in the loop): proves this CMake
# build links `libRNA.a` correctly. The full behavior contract (parity,
# malformed input, memory, ABI) lives in `tests/test_vienna_binding.py`.
enable_testing()
add_executable(test_puzzler_smoke src/vienna_layout/tests/test_puzzler_smoke.cpp)
target_include_directories(test_puzzler_smoke PRIVATE
  "${VIENNA_PREFIX}/include"
  "${VIENNA_PREFIX}/include/ViennaRNA"
)
target_compile_options(test_puzzler_smoke PRIVATE -fno-lto)
target_link_options(test_puzzler_smoke PRIVATE -fno-lto)
target_link_libraries(test_puzzler_smoke PRIVATE
  "${VIENNA_RNA_LIB}"
  Threads::Threads
  stdc++
  m
)
add_test(NAME puzzler_smoke COMMAND test_puzzler_smoke)
