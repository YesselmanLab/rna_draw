/**
 * @file config_changes.cpp
 * @brief `check_and_apply_config_changes`, ported from
 *        `checkAndApplyConfigChanges` (`handleConfigChanges.inc:55`) -- the
 *        resolver's sole gate between "a config delta was computed" and
 *        "the tree reflects it".
 */

#include <cmath>

#include "resolve_internal.hpp"
#include "rna_layout/config.hpp"
#include "rna_layout/geometry.hpp"

namespace rna_layout {

namespace {

/// Fix too-small changes in place: if every `delta_cfg` entry has magnitude
/// below `geom::kEpsilon3`, double every entry, up to 100 times. Ported from
/// `checkAndApplyConfigChanges`'s inline "fix deltas if changes are too
/// small" block (`handleConfigChanges.inc:76-98`) -- micro-changes grow a
/// loop's radius without gaining enough separation to actually resolve the
/// intersection, so this nudges them up to a resolvable magnitude
/// (0.1 degrees, min) before the validity/apply step ever sees them.
void fix_too_small_changes(std::vector<double>& delta_cfg) {
  constexpr int kMaxIterations = 100;
  for (int attempt = 0; attempt < kMaxIterations; ++attempt) {
    bool valid = false;
    for (double delta : delta_cfg) {
      if (std::fabs(delta) >= geom::kEpsilon3) {
        valid = true;
        break;
      }
    }
    if (valid) {
      break;
    }
    for (double& delta : delta_cfg) {
      delta *= 2.0;
    }
  }
}

}  // namespace

bool check_and_apply_config_changes(TreeNode& tree, std::vector<double>& delta_cfg,
                                    IntersectionType it, double unpaired, double paired,
                                    ResolverState& state) {
  fix_too_small_changes(delta_cfg);

  const bool accepted =
      cfg_is_valid(*tree.cfg, delta_cfg);  // NOLINT(bugprone-unchecked-optional-access)
  ++state.changes_applied;
  if (accepted) {
    // `radius_new = -1.0`: "unknown -- calculate optimal radius" (grow-only
    // sentinel), matching `checkAndApplyConfigChanges`'s own call
    // (`handleConfigChanges.inc:106-107`).
    apply_changes_to_config_and_bounding_boxes(tree, delta_cfg, /*radius_new=*/-1.0, paired,
                                               unpaired);
  }

  state.trace.push_back(ChangeTraceEntry{tree.id, it, delta_cfg, accepted});
  return accepted;
}

}  // namespace rna_layout
