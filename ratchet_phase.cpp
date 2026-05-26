#include "phylotree.h"
#include "phylonode.h"
#include "ratchet_phase.h"
#include "optimizer.h"
#include "spr_phase.h"
#include "spr_utils.h"
#include "spr_delta_exact.h"
#include "nucleotide_utils.h"

#include <chrono>
#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;
using namespace std::chrono;

// Ratchet acceptance + pattern-reweight tunables.
static const double   RATCHET_INV_TEMPERATURE      = 2.0;
static const int      RATCHET_ZERO_ACCEPT_PCT      = 50;
// Each pattern has a 1-in-(MASK+1) chance of being scaled by REWEIGHT_FACTOR.
static const unsigned RATCHET_REWEIGHT_DENOM_MASK  = 3u;
static const int      RATCHET_REWEIGHT_FACTOR      = 2;

int runRatchetPhase(PlacementOptimizer* optimizer,
                    PhyloTree* tree,
                    const PlacementOptimizeOptions& opts,
                    int entry_score,
                    int max_id,
                    high_resolution_clock::time_point start_time) {
    int ratchet_iters = opts.ratchet_iters;
    int ratchet_seed  = opts.ratchet_seed;
    int ratchet_runs  = opts.ratchet_runs;
    if (ratchet_runs < 1) ratchet_runs = 1;
    int max_radius    = opts.max_radius;
    double wall_seconds = opts.wall_seconds;

    auto optimizer_total_elapsed = [&]() {
        return duration_cast<milliseconds>(
            high_resolution_clock::now() - start_time).count() / 1000.0;
    };
    auto wall_exceeded = [&]() {
        return wall_seconds > 0.0 && optimizer_total_elapsed() > wall_seconds;
    };
    auto remaining_budget = [&]() -> double {
        if (wall_seconds <= 0.0) return 0.0;
        double remaining_seconds = wall_seconds - optimizer_total_elapsed();
        return (remaining_seconds < 0.5) ? 0.5 : remaining_seconds;
    };

    cout << "\n=== Ratchet phase (" << ratchet_iters << " iters"
         << (ratchet_runs > 1 ? " x " + std::to_string(ratchet_runs) + " runs" : "")
         << ", seed=" << ratchet_seed << ") ===" << endl;
    int pre_ratchet = entry_score;
    int tracked_score = entry_score;

    std::vector<int> orig_freq;
    tree->fitchSnapshotPatternFreq(orig_freq);
    int nptn = (int)orig_freq.size();

    vector<PhyloNode*> all_nodes_for_save = collectAllNodes(tree, max_id);
    auto saveAll = [&](vector<NeighborSave>& out) {
        out.clear();
        out.reserve(all_nodes_for_save.size() * 3);
        for (PhyloNode* n : all_nodes_for_save) sprSaveNodeTopology(n, out);
    };

    vector<NeighborSave> pre_ratchet_topology;
    saveAll(pre_ratchet_topology);
    std::vector<nuc_one_hot> pre_ratchet_major = tree->fitchSaveNodeMajor();
    PhyloTree::FitchAuxSnapshot pre_ratchet_aux = tree->fitchSaveAux();

    vector<NeighborSave> overall_best_topology;
    std::vector<nuc_one_hot> overall_best_major;
    PhyloTree::FitchAuxSnapshot overall_best_aux;
    int overall_best_score = tracked_score;
    if (ratchet_runs > 1) {
        overall_best_topology = pre_ratchet_topology;
        overall_best_major = pre_ratchet_major;
        overall_best_aux = pre_ratchet_aux;
    }

    for (int run = 0; run < ratchet_runs; run++) {
        if (ratchet_runs > 1 && run > 0) {
            cout << "-- Run " << (run + 1) << "/" << ratchet_runs << " --" << endl;
            sprUndoTopology(pre_ratchet_topology);
            tree->fitchRestoreNodeMajor(pre_ratchet_major);
            tree->fitchRestoreAux(pre_ratchet_aux);
            orientTreeToRoot(tree, max_id);
            SPRDeltaExact::precomputeDepths(tree);
            tracked_score = pre_ratchet;
        }

        vector<NeighborSave> best_saves;
        saveAll(best_saves);
        int best_score = tracked_score;

        std::mt19937 rng((unsigned)(ratchet_seed + run * 101));

        bool tree_at_best = true;

        for (int it = 0; it < ratchet_iters; it++) {
            if (wall_exceeded()) {
                cout << "Ratchet stopping at iter " << it
                     << ": optimizer wall cap (" << wall_seconds << "s) reached" << endl;
                break;
            }
            int reweighted = 0;
            for (int ptn = 0; ptn < nptn; ptn++) {
                if ((rng() & RATCHET_REWEIGHT_DENOM_MASK) == 0u) {
                    tree->fitchScalePatternFreq(ptn, RATCHET_REWEIGHT_FACTOR);
                    reweighted++;
                }
            }

            int weighted_score = tree->fitchRecomputeWithDiffs();
            orientTreeToRoot(tree, max_id);
            SPRDeltaExact::precomputeDepths(tree);
            bool capped = false;
            for (int r = 1; ; r *= 2) {
                r = min(r, max_radius);
                int after = optimizer->optimizeAtRadius(r, weighted_score, remaining_budget());
                if (after < weighted_score) weighted_score = after;
                if (r == max_radius) break;
                if (wall_exceeded()) { capped = true; break; }
            }

            tree->fitchRestorePatternFreq(orig_freq);
            int unweighted_score = tree->fitchRecomputeWithDiffs();
            orientTreeToRoot(tree, max_id);
            SPRDeltaExact::precomputeDepths(tree);

            for (int r = 1; ; r *= 2) {
                r = min(r, max_radius);
                int after = optimizer->optimizeAtRadius(r, unweighted_score, remaining_budget());
                if (after < unweighted_score) unweighted_score = after;
                if (r == max_radius) break;
                if (wall_exceeded()) { capped = true; break; }
            }

            cout << "Ratchet " << (it + 1) << ": reweighted=" << reweighted
                 << " unweighted=" << unweighted_score
                 << " best=" << best_score
                 << (capped ? " [CAPPED]" : "") << endl;

            int delta = unweighted_score - best_score;
            bool accept;
            if (delta < 0) {
                // improvement: always accept (greedy half).
                accept = true;
            } else if (delta == 0) {
                // plateau: 50/50 coin.
                accept = ((int)(rng() % 100u) < RATCHET_ZERO_ACCEPT_PCT);
            } else {
                // worse: accept with exp(-delta * beta). delta=1 ~13.5%, delta=2 ~1.8%,
                // delta>=3 effectively never. Lets the chain forget tiny mistakes.
                double prob = std::exp(-(double)delta * RATCHET_INV_TEMPERATURE);
                double draw = (double)rng() / (double)std::mt19937::max();
                accept = (draw < prob);
            }

            if (delta < 0) {
                best_score = unweighted_score;
                saveAll(best_saves);
                tracked_score = unweighted_score;
                tree_at_best = true;
            } else if (accept) {
                tracked_score = unweighted_score;
                tree_at_best = false;
            } else {
                sprUndoTopology(best_saves);
                tree->fitchRecomputeWithDiffs();
                orientTreeToRoot(tree, max_id);
                SPRDeltaExact::precomputeDepths(tree);
                tracked_score = best_score;
                tree_at_best = true;
            }
        }
        if (!tree_at_best) {
            sprUndoTopology(best_saves);
            tree->fitchRecomputeWithDiffs();
            orientTreeToRoot(tree, max_id);
            SPRDeltaExact::precomputeDepths(tree);
        }
        if (ratchet_runs > 1) {
            cout << "Run " << (run + 1) << " result: " << pre_ratchet
                 << " -> " << best_score << endl;
        } else {
            cout << "Ratchet phase: " << pre_ratchet << " -> " << best_score
                 << " (total delta=" << (pre_ratchet - best_score) << ")" << endl;
        }
        tracked_score = best_score;

        if (best_score < overall_best_score) {
            overall_best_score = best_score;
            saveAll(overall_best_topology);
            overall_best_major = tree->fitchSaveNodeMajor();
            overall_best_aux = tree->fitchSaveAux();
        }
    }  // end runs loop

    if (ratchet_runs > 1) {
        sprUndoTopology(overall_best_topology);
        tree->fitchRestoreNodeMajor(overall_best_major);
        tree->fitchRestoreAux(overall_best_aux);
        orientTreeToRoot(tree, max_id);
        SPRDeltaExact::precomputeDepths(tree);
        tracked_score = overall_best_score;
        cout << "Ratchet best across " << ratchet_runs << " runs: "
             << pre_ratchet << " -> " << overall_best_score
             << " (delta=" << (pre_ratchet - overall_best_score) << ")" << endl;
    }

    return tracked_score;
}
