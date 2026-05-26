#ifndef RATCHET_PHASE_H
#define RATCHET_PHASE_H

#include <chrono>

class PhyloTree;
class PlacementOptimizer;
struct PlacementOptimizeOptions;

/**
 * Parsimony ratchet (Nixon 1999): perturb the alignment weights, run SPR
 * to convergence on the weighted dataset, then unweight and SPR-converge
 * again. Best-of-K runs supported via opts.ratchet_runs.
 * @param optimizer    PlacementOptimizer instance for calling optimizeAtRadius.
 * @param tree         Tree to optimize.
 * @param opts         Configuration (uses ratchet_iters, ratchet_seed,
 *                     ratchet_runs, max_radius, wall_seconds).
 * @param entry_score  Tree's score on entry; ratchet cannot return worse.
 * @param max_id       Max node id (passed through to scratch helpers).
 * @param start_time   Wall-budget anchor: opts.wall_seconds is relative to this.
 * @return Score after ratchet (<= entry_score).
 */
int runRatchetPhase(PlacementOptimizer* optimizer,
                    PhyloTree* tree,
                    const PlacementOptimizeOptions& opts,
                    int entry_score,
                    int max_id,
                    std::chrono::high_resolution_clock::time_point start_time);

#endif // RATCHET_PHASE_H
