#ifndef SPROPTIMIZE_H
#define SPROPTIMIZE_H

#include <vector>
#include <map>

class Node;
class PhyloTree;
class PhyloNode;
class PhyloNeighbor;
class Fitch;

class SPROptimizer {
public:
    SPROptimizer(PhyloTree* tree);
    ~SPROptimizer();

    /** Run full optimization. Returns final parsimony score. */
    int optimizeTree(int max_passes = 1, int max_radius = 32);

    /** Run SPR rounds at a given radius. Returns score after optimization. */
    int optimizeAtRadius(int radius, Fitch& fitch);

private:
    PhyloTree* tree;
    int current_parsimony_score;
    int best_score_seen;
};

#endif
