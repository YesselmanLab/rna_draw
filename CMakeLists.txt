cmake_minimum_required(VERSION 3.15...3.27)
# A fixed project name (not `${SKBUILD_PROJECT_NAME}`) so this also
# configures cleanly for a standalone `cmake -B build && cmake --build
# build && ctest --test-dir build` invocation, not only through
# scikit-build-core's `pip install -e .`.
project(rna_draw LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

# Milestone A (cpp-layout-core) step 11 -- RETIRED from the default build:
# the vendored ViennaRNA C (`_vienna_layout` + its parity-instrumentation
# TU) now builds ONLY when this option is explicitly turned ON. The owned
# `rna_layout` native port (`include/rna_layout`, `src/layout_core/`) is
# the shipped production engine (`rna_draw.layout.production`); the
# default/shipped build needs no ViennaRNA header or runtime dependency at
# all. Turn this ON for parity/regression testing against the vendored
# oracle (`tests/test_native_parity.py`, `tests/test_vienna_binding.py`) --
# see `.claude/plans/current-plan.md`, step 11.
option(RNA_DRAW_BUILD_ORACLE
  "Compile the vendored-C parity oracle (_vienna_layout + vendor_instrument.c) \
for differential testing against _layout_core. OFF by default -- the shipped \
build does not need ViennaRNA." OFF)

# FP REPRODUCIBILITY (cpp-reviewer refinement #1): pin matched
# floating-point flags across BOTH the native library/module and the
# vendored oracle build, so the parity gates (turtle/box coordinate
# agreement <=1e-6; criterion-3 detection-set equality, a thresholded
# compare at the epsilonRecognize boundary, definitions.inc:59) cannot
# FALSE-FAIL from FMA contraction or expression-reassociation differences
# between the two "sides" of a same-compiler, same-flags build. Applied to
# every native/oracle target below via `RNA_DRAW_FP_FLAGS`.
# Residual note: `asin`/`acos`/`sqrt`/`cos`/`sin` calls are still routed
# through the platform libm on both sides (same libm here, since both
# targets compile in the same toolchain invocation) -- a genuine
# cross-libm build would need a documented small transcendental tolerance
# (a few ULP, i.e. ~1e-15 relative) on top of `-ffp-contract=off`; the
# parity harness's 1e-6 turtle gate already has ample headroom over that.
set(RNA_DRAW_FP_FLAGS -ffp-contract=off)

# pybind11 is a build-system requirement (pyproject.toml [build-system]);
# locate its CMake config via the interpreter that is actually building us,
# rather than relying on scikit-build-core to pre-populate CMAKE_PREFIX_PATH.
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
  OUTPUT_VARIABLE pybind11_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
find_package(pybind11 CONFIG REQUIRED)

enable_testing()

# ---------------------------------------------------------------------------
# _vienna_layout: the vendored ViennaRNA parity ORACLE. RNA_DRAW_BUILD_ORACLE
# -only (Milestone A step 11) -- not part of the default/shipped build. Every
# target, include path, and the CONDA_PREFIX/ViennaRNA-header requirement
# below is scoped to this option so a default `cmake -B build` needs neither
# ViennaRNA nor an active conda/micromamba environment.
# ---------------------------------------------------------------------------
if(RNA_DRAW_BUILD_ORACLE)
  # ViennaRNA HEADERS live in the active conda/micromamba environment --
  # confirmed present at $CONDA_PREFIX/include (see plan Step 0). Only the
  # headers are needed: the layout core is vendored and compiled standalone
  # (see below), so libRNA.a is NOT a link input and there is no ViennaRNA
  # runtime dependency even for this oracle-only target.
  if(NOT DEFINED ENV{CONDA_PREFIX})
    message(FATAL_ERROR
      "CONDA_PREFIX is not set -- RNA_DRAW_BUILD_ORACLE=ON needs ViennaRNA "
      "headers; run the build via 'micromamba run -n py3 ...', or turn this "
      "option OFF for the default (no-ViennaRNA) build.")
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
  # layout core links WITHOUT libRNA.a: this oracle target has no ViennaRNA
  # runtime dependency either, only a build-time header dependency.
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

  # The parity-instrumentation TU (see that file's header for the
  # macro-interposition mechanism it names). Never edits vendored .inc
  # logic; only compiled into this same module so `_vienna_layout` can
  # expose oracle-only dumps.
  target_sources(_vienna_layout PRIVATE src/vienna_layout/vendor_instrument.c)
  target_compile_definitions(_vienna_layout PRIVATE RNA_DRAW_BUILD_ORACLE)
  set_source_files_properties(src/vienna_layout/vendor_instrument.c PROPERTIES
    LANGUAGE C
    COMPILE_OPTIONS "-w"
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
  # matching it avoids an LTO-mismatch link failure. `RNA_DRAW_FP_FLAGS`: see
  # this file's top-of-file comment (FP reproducibility for the parity gates).
  target_compile_options(_vienna_layout PRIVATE -fno-lto ${RNA_DRAW_FP_FLAGS})
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
  add_executable(test_puzzler_smoke
    src/vienna_layout/tests/test_puzzler_smoke.cpp
    ${VENDOR_PUZZLER_SOURCES}
  )
  target_include_directories(test_puzzler_smoke PRIVATE
    "${VIENNA_PREFIX}/include"
    "${VIENNA_PREFIX}/include/ViennaRNA"
    "${VENDOR_PUZZLER_DIR}"
  )
  target_compile_options(test_puzzler_smoke PRIVATE -fno-lto ${RNA_DRAW_FP_FLAGS})
  target_link_options(test_puzzler_smoke PRIVATE -fno-lto)
  target_link_libraries(test_puzzler_smoke PRIVATE
    stdc++
    m
  )
  add_test(NAME puzzler_smoke COMMAND test_puzzler_smoke)
endif()

# ---------------------------------------------------------------------------
# rna_layout: the owned modern-C++ port (Milestone A, cpp-layout-core).
# A parallel, isolated tree (include/rna_layout, src/layout_core) -- no
# ViennaRNA header/symbol dependency, validated against `_vienna_layout`
# above as a differential parity oracle. See
# `.claude/plans/current-plan.md` for the full port plan; this slice adds
# the config-tree/bounding-box/wedge construction (Milestone A step 4) on
# top of the foundation value types + geometry/config helpers + the
# turtle-base engine (`rna_layout::layout_turtle`).
# ---------------------------------------------------------------------------
set(RNA_LAYOUT_CORE_SOURCES
  src/layout_core/geometry.cpp
  src/layout_core/config.cpp
  src/layout_core/pair_table.cpp
  src/layout_core/turtle.cpp
  src/layout_core/bounding_boxes.cpp
  src/layout_core/config_tree.cpp
  src/layout_core/bounding_wedge.cpp
  src/layout_core/intersections.cpp
  src/layout_core/intersect_boxes.cpp
  src/layout_core/intersect_tree.cpp
  src/layout_core/nucleotide_coords.cpp
  src/layout_core/resolve_exterior.cpp
  src/layout_core/deltas.cpp
  src/layout_core/config_changes.cpp
  src/layout_core/resolve_siblings.cpp
  src/layout_core/rotation_angle.cpp
  src/layout_core/rotation_angle2.cpp
  src/layout_core/exterior_boxes.cpp
  src/layout_core/resolve_ancestors.cpp
  src/layout_core/optimize.cpp
  src/layout_core/optimize2.cpp
  src/layout_core/resolve.cpp
  src/layout_core/puzzler.cpp
  src/layout_core/debug_dump.cpp
)

add_library(rna_layout_core STATIC ${RNA_LAYOUT_CORE_SOURCES})
target_include_directories(rna_layout_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_compile_options(rna_layout_core PRIVATE ${RNA_DRAW_FP_FLAGS})

pybind11_add_module(_layout_core src/layout_core/bindings.cpp)
target_link_libraries(_layout_core PRIVATE rna_layout_core)
target_compile_options(_layout_core PRIVATE ${RNA_DRAW_FP_FLAGS})

install(TARGETS _layout_core DESTINATION rna_draw)

# ctest unit tests: pure C++, no Python/pybind11 in the loop (parity-vs-
# oracle comparisons run at the Python level, `tests/test_native_parity.py`,
# since that is where both compiled extensions are importable together).
foreach(rna_layout_test_name geometry_test config_test turtle_test tree_test
                             bounding_boxes_test config_tree_test bounding_wedge_test
                             intersections_test intersect_tree_test puzzler_test
                             deltas_test resolve_siblings_test resolve_ancestors_test
                             optimize_test)
  add_executable(${rna_layout_test_name} src/layout_core/tests/${rna_layout_test_name}.cpp)
  target_link_libraries(${rna_layout_test_name} PRIVATE rna_layout_core)
  target_compile_options(${rna_layout_test_name} PRIVATE ${RNA_DRAW_FP_FLAGS})
  add_test(NAME ${rna_layout_test_name} COMMAND ${rna_layout_test_name})
endforeach()
