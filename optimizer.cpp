#include "phylotree.h"
#include "phylonode.h"
#include "optimizer.h"
#include "spr_phase.h"
#include "ratchet_phase.h"
#include "tbr_phase.h"
#include "spr_context.h"
#include "spr_delta_exact.h"

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace std;
using namespace std::chrono;

PlacementOptimizer::PlacementOptimizer(PhyloTree* tree) : tree(tree) {
    if (!tree)        throw std::invalid_argument("PlacementOptimizer: tree cannot be NULL");
    if (!tree->root)  throw std::invalid_argument("PlacementOptimizer: tree->root must be set");
    if (!tree->aln)   throw std::invalid_argument("PlacementOptimizer: tree->aln must be set");
}

PlacementOptimizer::~PlacementOptimizer() {}

int PlacementOptimizer::optimizeAtRadius(int radius, int known_score, double wall_seconds) {
    return sprOptimizeAtRadius(tree, radius, known_score, wall_seconds);
}

int PlacementOptimizer::optimizeTree(const PlacementOptimizeOptions& opts) {
    int max_passes      = opts.max_passes;
    int max_radius      = opts.max_radius;
    int ratchet_iters   = opts.ratchet_iters;
    int ratchet_runs    = opts.ratchet_runs < 1 ? 1 : opts.ratchet_runs;
    double wall_seconds = opts.wall_seconds;

    setActiveSPRTree(tree);

    auto t0 = high_resolution_clock::now();
    int initial_score = tree->fitchRunForSPR();
    auto t1 = high_resolution_clock::now();

    cout << "Fitch: score=" << initial_score
         << " (max_passes=" << max_passes
         << ", max_radius=" << max_radius
         << ", init=" << duration_cast<milliseconds>(t1 - t0).count() << "ms)" << endl;

    int tracked_score = initial_score;
    int max_id = tree->fitchMaxNodeId();
    orientTreeToRoot(tree, max_id);
    SPRDeltaExact::precomputeDepths(tree);

    auto optimizer_total_start = high_resolution_clock::now();
    auto elapsed = [&]() {
        return duration_cast<milliseconds>(
            high_resolution_clock::now() - optimizer_total_start).count() / 1000.0;
    };
    auto wall_exceeded = [&]() { return wall_seconds > 0.0 && elapsed() > wall_seconds; };
    auto remaining = [&]() -> double {
        if (wall_seconds <= 0.0) return 0.0;
        double remaining_seconds = wall_seconds - elapsed();
        return (remaining_seconds < 0.5) ? 0.5 : remaining_seconds;
    };

    // SPR phase: radius-escalating passes, strict improvement.
    for (int pass = 0; pass < max_passes; pass++) {
        if (wall_exceeded()) break;
        int start = tracked_score;
        int consecutive_empty = 0;
        int cur = start;

        for (int r = 1; ; r *= 2) {
            r = min(r, max_radius);
            int after = optimizeAtRadius(r, cur, remaining());
            if (after >= cur) {
                consecutive_empty++;
                if (consecutive_empty >= RADIUS_STALL_LIMIT && r < max_radius) {
                    cout << "  " << consecutive_empty
                         << " consecutive radii with no improvement, skipping higher" << endl;
                    break;
                }
            } else {
                consecutive_empty = 0;
                cur = after;
            }
            if (r == max_radius) break;
            if (wall_exceeded()) break;
        }

        int end = cur;
        tracked_score = end;
        double improvement = (start > 0) ? (double)(start - end) / start : 0.0;
        cout << "Pass " << pass + 1 << ": " << start << " -> " << end
             << " (delta=" << (start - end)
             << ", improvement=" << fixed << setprecision(4)
             << (improvement * 100) << "%)" << endl;
        if (end >= start) break;
    }

    // Ratchet phase (Nixon 1999) — disabled when ratchet_iters == 0.
    if (ratchet_iters > 0 && !wall_exceeded()) {
        tracked_score = runRatchetPhase(this, tree, opts, tracked_score, max_id,
                                         optimizer_total_start);
    }

    // TBR phase — disabled when tbr_iters == 0.
    if (opts.tbr_iters > 0 && !wall_exceeded()) {
        tracked_score = runTBRPhase(this, tree, opts, tracked_score, max_id,
                                     optimizer_total_start);
    }

    cout << "\nPlacement optimize complete: " << initial_score << " -> " << tracked_score
         << " (total delta=" << (initial_score - tracked_score) << ")" << endl;

    setActiveSPRTree(nullptr);
    return tracked_score;
}
