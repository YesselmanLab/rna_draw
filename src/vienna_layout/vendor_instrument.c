/*
 * vendor_instrument.c -- oracle-side parity instrumentation for the
 * vendored RNApuzzler/RNAturtle layout core (vendor/RNApuzzler/), compiled
 * ONLY when RNA_DRAW_BUILD_ORACLE is on (see CMakeLists.txt).
 *
 * MECHANISM (named per the plan review's refinement #2 -- see
 * .claude/plans/current-plan.md, "cpp-reviewer APPROVED" section): the
 * vendored layout files are `#include`d `.inc` amalgams, each compiled
 * into ONE translation unit per engine (RNApuzzler.c / RNAturtle.c already
 * do this). Their tree/box/config-change internals (configtree.inc,
 * handleConfigChanges.inc, ...) are PRIVATE (`static`) functions, invisible
 * from any other TU -- a plain wrapper TU (like this one) CANNOT call them
 * directly. To dump them WITHOUT editing vendored logic, later steps
 * compile a SECOND copy of the relevant .inc(s) into THIS TU using
 * MACRO INTERPOSITION: `#define <privateFunctionName> <shimName>` before
 * `#include`-ing the .inc, so every call the .inc's own code makes to that
 * function resolves to the shim instead, which can record its arguments
 * and then invoke the original body (still compiled, just under the shim's
 * name). This is the same technique the rna_draw fork already uses for
 * `rnadraw_clearance_value` (definitions.inc:39) and the
 * `max_config_changes` patch in RNApuzzler.c -- no vendored .inc file is
 * edited; only an oracle-only TU that never ships in a production build
 * includes it differently.
 *
 * THIS SLICE (Milestone A, turtle-base, steps 1-3): turtle has no private
 * tree/box state to intercept -- `vrna_plot_coords_turtle`/`_pt` are
 * already PUBLIC entry points, dumped directly by bindings.cpp's
 * `dump_turtle` (a plain alias of `plot_coords_turtle`, wired there, not
 * here). This file is a placeholder that (a) proves the
 * RNA_DRAW_BUILD_ORACLE CMake wiring compiles and links, and (b) names the
 * mechanism `dump_tree`/`dump_detections`/`dump_change_trace` will use once
 * the native side has a tree/detections/trace to compare against
 * (Milestone A steps 4, 5, and 7 respectively).
 */

const char *
rnadraw_oracle_instrumentation_version(void)
{
  return "vendor_instrument v0: turtle dump only (see this file's header "
         "for the macro-interposition mechanism planned tree/box/detection/"
         "change-trace dumps will use)";
}
