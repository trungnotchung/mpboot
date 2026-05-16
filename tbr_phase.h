#ifndef TBR_PHASE_H
#define TBR_PHASE_H

#include "phylonode.h"
#include <chrono>
#include <vector>

class PhyloTree;
class PlacementOptimizer;
struct PlacementOptimizeOptions;

/**
 * State produced by bisectInternalEdge() and consumed by reconnectBisected()
 * / undoReconnect(). Records the four siblings short-circuited on each side
 * of the bisected edge so the original topology can be rebuilt.
 *
 * Layout after bisectInternalEdge:
 *
 *      bisect_a_sib1                 bisect_b_sib1
 *           \                      /
 *           bisect_a [removed] bisect_b
 *           /                      \
 *      bisect_a_sib2                 bisect_b_sib2
 *
 * After bisection: bisect_a_sib1<->bisect_a_sib2 directly; bisect_b_sib1<->bisect_b_sib2
 * directly. bisect_a / bisect_b keep their internal Neighbor objects (so we
 * can later re-insert them) but are no longer reachable from the rest of the tree.
 */
struct BisectionState {
    PhyloNode* bisect_a;
    PhyloNode* bisect_b;
    PhyloNode* bisect_a_sib1;
    PhyloNode* bisect_a_sib2;
    PhyloNode* bisect_b_sib1;
    PhyloNode* bisect_b_sib2;
    bool valid;               ///< False if the requested edge is not bisectable.

    BisectionState() : bisect_a(nullptr), bisect_b(nullptr),
                       bisect_a_sib1(nullptr), bisect_a_sib2(nullptr),
                       bisect_b_sib1(nullptr), bisect_b_sib2(nullptr), valid(false) {}
};

/**
 * Bisect on the internal edge between bisect_a and bisect_b. Both endpoints
 * must be non-root, non-leaf, binary internal nodes that are currently
 * connected. Orientation-agnostic.
 * @return BisectionState with valid=true on success.
 */
BisectionState bisectInternalEdge(PhyloTree* tree, PhyloNode* bisect_a,
                                  PhyloNode* bisect_b);

/**
 * Re-insert bisect_a on (reattach_a_u, reattach_a_v) in subtree A and bisect_b on
 * (reattach_b_u, reattach_b_v) in subtree B, restoring the bisect_a<->bisect_b link.
 * Pass (s.bisect_a_sib1, s.bisect_a_sib2, s.bisect_b_sib1, s.bisect_b_sib2) for a
 * no-op round-trip.
 */
void reconnectBisected(PhyloTree* tree, const BisectionState& state,
                       PhyloNode* reattach_a_u, PhyloNode* reattach_a_v,
                       PhyloNode* reattach_b_u, PhyloNode* reattach_b_v);

/**
 * Inverse of reconnectBisected. The (reattach_a_u, reattach_a_v, reattach_b_u, reattach_b_v)
 * arguments MUST match the most recent reconnectBisected call.
 */
void undoReconnect(PhyloTree* tree, const BisectionState& state,
                   PhyloNode* reattach_a_u, PhyloNode* reattach_a_v,
                   PhyloNode* reattach_b_u, PhyloNode* reattach_b_v);

/**
 * Score a candidate TBR move via apply→Fitch→undo (correctness-first;
 * no closed-form delta yet — see plan/PLAN_TBR_OPTIMIZER.md).
 * @param pre_score  Score before the move (caller-supplied; not recomputed).
 * @return post_score - pre_score, or INT_MAX if the bisection is invalid.
 */
int evaluateTBRMove(PhyloTree* tree, int pre_score,
                    PhyloNode* bisect_a, PhyloNode* bisect_b,
                    PhyloNode* reattach_a_u, PhyloNode* reattach_a_v,
                    PhyloNode* reattach_b_u, PhyloNode* reattach_b_v);

/**
 * "Best single improving move per round" TBR hill climber. Each round
 * enumerates internal binary edges, bisects, BFS-enumerates alt edges within
 * tbr_max_radius hops in each subtree, scores via evaluateTBRMove, applies
 * the single best improving move, and repeats.
 * @param tbr_max_radius  Max BFS depth (used for both subtrees A and B).
 *                        2+ enables real TBR moves.
 */
int optimizeTBRAtRadius(PhyloTree* tree, int tbr_max_radius,
                        int known_score = 0, double wall_seconds = 0.0,
                        int max_rounds = 1000000);

/**
 * TBR phase driver: alternates one best TBR move at radius opts.tbr_max_radius
 * with full SPR re-convergence (radius-doubling 1..max_radius). Repeats
 * opts.tbr_iters times or until neither phase improves.
 * @param optimizer    PlacementOptimizer instance for calling optimizeAtRadius.
 * @param tree         Tree to optimize.
 * @param opts         Configuration (uses tbr_iters, tbr_max_radius,
 *                     max_radius, wall_seconds).
 * @param entry_score  Score on entry.
 * @param max_id       Max node id (passed through to scratch helpers).
 * @param start_time   Wall-budget anchor.
 * @return Score after TBR phase (<= entry_score).
 */
int runTBRPhase(PlacementOptimizer* optimizer,
                PhyloTree* tree,
                const PlacementOptimizeOptions& opts,
                int entry_score,
                int max_id,
                std::chrono::high_resolution_clock::time_point start_time);

#endif // TBR_PHASE_H
