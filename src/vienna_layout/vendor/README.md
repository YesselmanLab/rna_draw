# Vendored: RNApuzzler / RNAturtle (standalone, no libRNA.a)

Source: [ViennaRNA 2.7.0](https://github.com/ViennaRNA/ViennaRNA/releases/tag/v2.7.0),
`src/ViennaRNA/plotting/RNApuzzler/` (`RNApuzzler.c`, `RNAturtle.c` and
their `includes/*.inc` / `headers/*.h`).

License waived for this project by the project owner (upstream ViennaRNA is
GPLv2+; this in-tree copy is used and modified here with the owner's
permission rather than under a redistributable license grant).

This copy is intentionally editable and will diverge from upstream over
time as `rna_draw` modifies the layout algorithm (see
`src/vienna_layout/bindings.cpp` for the two resolver levers exposed on top
of it: `allow_flipping` and `max_config_changes`).

## Standalone build (libRNA.a UNLINKED)

`RNApuzzler.c` + `RNAturtle.c` define every layout entry point the
extension calls (`vrna_plot_coords_puzzler[_pt]`,
`vrna_plot_options_puzzler[_free]`, `vrna_plot_coords_turtle[_pt]`). They
reference EXACTLY TWO ViennaRNA library symbols -- `vrna_alloc` and
`vrna_ptable` -- which are supplied by `vrna_compat.c` in this directory
(a ~120 LOC compat shim: faithful copies of upstream 2.7.0's `vrna_alloc`
and `vrna_ptable`/`extract_pairs`, plus a trivial `vrna_log` stderr stub
for the OOM/malformed-input reporter that the binding's pre-validation
never actually reaches). With that shim compiled in, `_vienna_layout`
links WITHOUT `libRNA.a` (see `CMakeLists.txt`) -- rna_draw has **no
ViennaRNA runtime dependency** for layout. Output is byte-identical to the
former libRNA.a-linked build (verified on `((((....))))` for puzzler and
turtle). Only the ViennaRNA *headers* (`$CONDA_PREFIX/include`) are still
needed, at compile time.

`naview` was removed: it lived only in `libRNA.a`'s non-reentrant
`naview.o` and per-structure benchmarking showed it never rescues a
puzzler-dirty structure -- a documented dead-end.

This is third-party C, kept as close to upstream as practical: it is
excluded from `clang-format`/`clang-tidy` and compiled with relaxed
warnings (see `CMakeLists.txt`). Do not reformat it wholesale; keep diffs
against upstream minimal and documented.
