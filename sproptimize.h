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
    int drift_iters   = 0;
    int drift_radius  = 0;    // 0 = inherit max_radius
    int ratchet_iters = 0;
    int ratchet_seed  = 42;
    int ratchet_runs  = 1;
    int cycles        = 1;
    int sector_size   = 0;    // 0 = disabled
    int sector_count  = 8;
    int sector_seed   = 42;
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
     * Runs main pass loop, then optional drift / ratchet / sectorial cycles.
     * @param opts Configuration (see SPROptimizeOptions).
     * @return Final parsimony score.
     */
    int optimizeTree(const SPROptimizeOptions& opts = SPROptimizeOptions());

    /**
     * Runs sector-restricted SPR (Goloboff sectorial search).
     * @param K_leaves     Target sector size in leaves.
     * @param n_sectors    Number of sector evaluations.
     * @param seed         RNG seed for sector centers.
     * @param max_radius   Max SPR radius (clamped to K_leaves).
     * @param wall_seconds Soft phase wall-clock budget (0 = no cap).
     * @param known_score  Current tree score (skip recompute if >0).
     * @return Score after all sectors processed.
     */
    int optimizeSectorial(int K_leaves, int n_sectors, int seed,
                           int max_radius, double wall_seconds = 0.0,
                           int known_score = 0);

    /**
     * Runs one radius pass (rounds of find-select-apply).
     * @param radius        Max SPR radius (0 = unbounded).
     * @param allow_drift   Accept zero-delta moves if true.
     * @param known_score   Current score (skip recompute if >0).
     * @param wall_seconds  Soft wall-clock cap (0 = no cap).
     * @param sector_member If non-null, restricts moves to sector members.
     * @return Score after the pass.
     */
    int optimizeAtRadius(int radius, bool allow_drift = false, int known_score = 0,
                          double wall_seconds = 0.0,
                          const std::vector<bool>* sector_member = nullptr);

private:
    PhyloTree* tree;
    int current_parsimony_score;
    int best_score_seen;
};

#endif
