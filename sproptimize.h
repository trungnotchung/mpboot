#ifndef SPROPTIMIZE_H
#define SPROPTIMIZE_H

#include <vector>
#include <map>

class Node;
class PhyloTree;
class PhyloNode;
class PhyloNeighbor;

class SPROptimizer {
public:
    SPROptimizer(PhyloTree* tree);
    ~SPROptimizer();

    int optimizeTree(int max_passes = 10, int max_radius = 32, int drift_iters = 0,
                     int drift_radius = 0, int ratchet_iters = 0, int ratchet_seed = 42,
                     int cycles = 1, int sector_size = 0, int sector_count = 8,
                     int sector_seed = 42, int ratchet_runs = 1);

    int optimizeSectorial(int K_leaves, int n_sectors, int seed,
                           int max_radius, double wall_seconds = 0.0,
                           int known_score = 0);

    int optimizeAtRadius(int radius, bool allow_drift = false, int known_score = 0,
                          double wall_seconds = 0.0,
                          const std::vector<bool>* sector_member = nullptr);

private:
    PhyloTree* tree;
    int current_parsimony_score;
    int best_score_seen;
};

#endif
