#ifndef SPROPTIMIZE_H
#define SPROPTIMIZE_H

#include <vector>

class PhyloTree;

struct SPROptimizeOptions {
    int max_passes    = 10;
    int max_radius    = 32;   // 0 = unbounded
    int ratchet_iters = 0;
    int ratchet_seed  = 42;
    int ratchet_runs  = 1;
    double wall_seconds = 0.0;  // <=0 = no cap (default unbounded)
    // TBR (Tree Bisection and Reconnection) post-SPR pass.
    // Disabled when tbr_iters == 0. tbr_max_radius is the BFS depth from the
    // bisection scar (used for both subtrees); >=2 enables actual TBR moves.
    int tbr_iters      = 0;
    int tbr_max_radius = 5;
};

class SPROptimizer {
public:
    /** @param tree Initialized tree (must have root and aln set). */
    SPROptimizer(PhyloTree* tree);
    ~SPROptimizer();

    /**
     * Runs the main radius-escalating SPR pass loop, then optional parsimony
     * ratchet (Nixon 1999) and TBR phases.
     * @param opts Configuration (see SPROptimizeOptions).
     * @return Final parsimony score.
     */
    int optimizeTree(const SPROptimizeOptions& opts = SPROptimizeOptions());

    /**
     * Runs one radius pass (rounds of find-select-apply, strict improvement).
     * Public for test_spr_unit.cpp; callers must keep the tree oriented-to-root
     * and depths precomputed between calls.
     * @param radius        Max SPR radius (0 = unbounded).
     * @param known_score   Current score (skip recompute if >0).
     * @param wall_seconds  Soft wall-clock cap (0 = no cap).
     * @return Score after the pass.
     */
    int optimizeAtRadius(int radius, int known_score = 0, double wall_seconds = 0.0);

private:
    PhyloTree* tree;
};

#endif
