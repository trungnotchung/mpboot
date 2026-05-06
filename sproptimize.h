#ifndef SPROPTIMIZE_H
#define SPROPTIMIZE_H

#include <vector>
#include <map>

class Node;
class PhyloTree;
class PhyloNode;
class PhyloNeighbor;

struct SPROptimizeOptions {
    int max_passes    = 10;
    int max_radius    = 32;   // 0 = unbounded
    int ratchet_iters = 0;
    int ratchet_seed  = 42;
    int ratchet_runs  = 1;
    double wall_seconds = 0.0;  // <=0 = no cap (default unbounded)
};

class SPROptimizer {
public:
    /**
     * @param tree Initialized tree
     */
    SPROptimizer(PhyloTree* tree);

    ~SPROptimizer();

    /**
     * Runs the main radius-escalating SPR pass loop, then an optional
     * parsimony ratchet (Nixon 1999) phase.
     * @param opts Configuration (see SPROptimizeOptions).
     * @return Final parsimony score.
     */
    int optimizeTree(const SPROptimizeOptions& opts = SPROptimizeOptions());

    /**
     * Runs one radius pass (rounds of find-select-apply, strict improvement).
     * @param radius        Max SPR radius (0 = unbounded).
     * @param known_score   Current score (skip recompute if >0).
     * @param wall_seconds  Soft wall-clock cap (0 = no cap).
     * @return Score after the pass.
     */
    int optimizeAtRadius(int radius, int known_score = 0, double wall_seconds = 0.0);

private:
    PhyloTree* tree;
    int current_parsimony_score;
    int best_score_seen;
};

#endif
