#ifndef SPR_PHASE_H
#define SPR_PHASE_H

#include "spr_utils.h"
#include "phylonode.h"
#include <vector>

class PhyloTree;

/** Max consecutive radii without improvement before bailing on the escalator. */
static const int RADIUS_STALL_LIMIT = 2;

/** Alias kept for historical readability in ratchet_phase / tbr_phase. */
using NeighborSave = SPRNeighborSave;

/**
 * Orients the tree by setting neighbors[0] = parent for every non-root node,
 * enabling O(1) parent lookup.
 */
void orientTreeToRoot(PhyloTree* tree, int max_id);

/** BFS-collects every reachable node from the tree's root. */
std::vector<PhyloNode*> collectAllNodes(PhyloTree* tree, int max_id);

/**
 * Run one SPR radius pass (rounds of find-select-apply, strict improvement).
 * Caller must keep the tree oriented-to-root and depths precomputed.
 *
 * @param tree         Tree to optimize.
 * @param radius       Max SPR radius (0 = unbounded).
 * @param known_score  Current score (skip recompute if > 0).
 * @param wall_seconds Soft wall-clock cap (0 = no cap).
 * @return Score after the pass.
 */
int sprOptimizeAtRadius(PhyloTree* tree, int radius,
                        int known_score = 0, double wall_seconds = 0.0);

#endif // SPR_PHASE_H
