/**
 * @file debug_dump.cpp
 * @brief Implements `debug_dump.hpp`.
 */

#include "rna_layout/debug_dump.hpp"

namespace rna_layout {

namespace {

void dump_tree_recursive(const TreeNode& node, int parent_id, std::vector<DumpTreeNode>& out) {
  const int id = static_cast<int>(out.size());

  DumpTreeNode entry;
  entry.id = id;
  entry.parent_id = parent_id;
  entry.loop_start = node.loop_start;
  entry.stem_start = node.stem_start;

  if (node.cfg.has_value()) {
    DumpConfig cfg;
    cfg.radius = node.cfg->radius;
    cfg.min_radius = node.cfg->min_radius;
    cfg.default_radius = node.cfg->default_radius;
    cfg.arcs.reserve(node.cfg->arcs.size());
    for (const ConfigArc& arc : node.cfg->arcs) {
      cfg.arcs.push_back(DumpConfigArc{arc.segments, arc.angle});
    }
    entry.cfg = std::move(cfg);
  }

  if (node.lbox.has_value()) {
    entry.lbox = DumpLoopBox{node.lbox->center.x, node.lbox->center.y, node.lbox->radius};
  }

  if (node.sbox.has_value()) {
    const StemBox& s = *node.sbox;
    entry.sbox = DumpStemBox{s.a.x,       s.a.y, s.b.x,
                             s.b.y,       s.c.x, s.c.y,
                             s.e.x,       s.e.y, static_cast<int>(s.bulges.size()),
                             s.bulge_dist};
  }

  out.push_back(std::move(entry));

  for (const auto& child : node.children) {
    dump_tree_recursive(*child, id, out);
  }
}

}  // namespace

std::vector<DumpTreeNode> dump_tree(const TreeNode& root) {
  std::vector<DumpTreeNode> out;
  dump_tree_recursive(root, /*parent_id=*/-1, out);
  return out;
}

}  // namespace rna_layout
