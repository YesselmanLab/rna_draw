#pragma once

// Shared ctest helper: build a fully-updated (T1) config tree from a
// dot-bracket structure, for `config_tree_test.cpp`/`bounding_boxes_test.cpp`/
// `bounding_wedge_test.cpp`. Not a production entry point -- the real one is
// `bindings.cpp`'s `dump_config_tree_binding` (`_layout_core.dump_tree`).

#include <memory>
#include <string>

#include "rna_layout/config_tree.hpp"
#include "rna_layout/pair_table.hpp"
#include "rna_layout/tree.hpp"
#include "rna_layout/turtle.hpp"

namespace rna_layout::test {

inline std::unique_ptr<TreeNode> build_updated_tree(const std::string& dot_bracket, double paired,
                                                    double unpaired) {
  const std::vector<int> pair_table = make_pair_table(dot_bracket);
  const TurtleLayout turtle = run_turtle_layout(pair_table, paired, unpaired);
  const double bulge_dist = stem_bulge_distance(unpaired);

  std::unique_ptr<TreeNode> tree =
      build_config_tree(pair_table, turtle.base_info, turtle.configs, turtle.coords, bulge_dist);
  update_bounding_boxes(*tree, paired, unpaired);
  return tree;
}

}  // namespace rna_layout::test
