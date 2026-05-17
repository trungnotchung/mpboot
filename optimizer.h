#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include <vector>

class PhyloTree;
struct BenchmarkStats;

struct PlacementOptimizeOptions {
    int    max_passes     = 10;    // Max SPR passes (0 = unlimited).
    int    max_radius     = 32;    // Max SPR radius (0 = unbounded).
    double wall_seconds   = 0.0;   // Wall-clock cap (<= 0 = no cap).
    int    ratchet_iters  = 0;     // Parsimony ratchet iterations (0 = disabled).
    int    ratchet_seed   = 42;    // RNG seed for ratchet reweighting.
    int    ratchet_runs   = 1;     // Best-of-K independent ratchet runs.
    int    tbr_iters      = 0;     // TBR rounds after SPR converges (0 = disabled).
    int    tbr_max_radius = 5;     // BFS depth from the bisection scar (used both sides).
};

class PlacementOptimizer {
public:
    /** @param tree Initialized tree. */
    PlacementOptimizer(PhyloTree* tree);
    ~PlacementOptimizer();

    /**
     * Runs the main radius-escalating SPR pass loop, then optional parsimony
     * ratchet (Nixon 1999) and TBR phases.
     * @param opts Configuration.
     * @param bench_stats Optional benchmark stats to track phase timing.
     * @return Final parsimony score.
     */
    int optimizeTree(const PlacementOptimizeOptions& opts = PlacementOptimizeOptions(),
                     BenchmarkStats* bench_stats = nullptr);

    /**
     * Run one radius pass.
     * @return Score after the pass.
     */
    int optimizeAtRadius(int radius, int known_score = 0, double wall_seconds = 0.0);

private:
    PhyloTree* tree;
};

#endif
