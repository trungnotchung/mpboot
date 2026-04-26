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

    /** Run full optimization. Returns final parsimony score. */
    int optimizeTree(int max_passes = 10, int max_radius = 32, int drift_iters = 0);

    /** Run SPR rounds at a given radius. Returns score after optimization.
     *  If known_score > 0, skips the initial recompute (caller guarantees fitch state is valid). */
    int optimizeAtRadius(int radius, bool allow_drift = false, int known_score = 0);

private:
    PhyloTree* tree;
    int current_parsimony_score;
    int best_score_seen;
};

#endif
