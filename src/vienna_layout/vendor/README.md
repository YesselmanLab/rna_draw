# Vendored: RNApuzzler / RNAturtle

Source: [ViennaRNA 2.7.0](https://github.com/ViennaRNA/ViennaRNA/releases/tag/v2.7.0),
`src/ViennaRNA/plotting/RNApuzzler/`.

License waived for this project by the project owner (upstream ViennaRNA is
GPLv2+; this in-tree copy is used and modified here with the owner's
permission rather than under a redistributable license grant).

This copy is intentionally editable and will diverge from upstream over
time as `rna_draw` modifies the layout algorithm (see
`src/vienna_layout/bindings.cpp` for the two resolver levers exposed on top
of it: `allow_flipping` and `max_config_changes`). It is compiled directly
into the `_vienna_layout` extension (see `CMakeLists.txt`) instead of being
linked from the prebuilt `libRNA.a`, so that its behavior can be changed.
Because it defines the same public symbols as the prebuilt archive
(`vrna_plot_coords_puzzler`, `vrna_plot_coords_puzzler_pt`,
`vrna_plot_options_puzzler`, `vrna_plot_options_puzzler_free`,
`vrna_plot_coords_turtle`, `vrna_plot_coords_turtle_pt`), the static linker
never pulls the stock `RNApuzzler.o`/`RNAturtle.o` objects out of
`libRNA.a` for those symbols -- this vendored copy always wins.

This is third-party C, kept as close to upstream as practical: it is
excluded from `clang-format`/`clang-tidy` and compiled with relaxed
warnings (see `CMakeLists.txt`). Do not reformat it wholesale; keep diffs
against upstream minimal and documented.
