#ifndef TBR_PRIMITIVES_H
#define TBR_PRIMITIVES_H

#include "phylonode.h"
#include <vector>

class PhyloTree;

/**
 * State produced by bisectInternalEdge() and consumed by reconnectBisected()
 * / undoReconnect(). Records the four siblings short-circuited on each side
 * of the bisected edge so the original topology can be rebuilt.
 *
 * Layout after bisectInternalEdge:
 *
 *      edge_p_sib1                 edge_q_sib1
 *           \                      /
 *           edge_p [removed] edge_p_back
 *           /                      \
 *      edge_p_sib2                 edge_q_sib2
 *
 * After bisection: edge_p_sib1<->edge_p_sib2 directly; edge_q_sib1<->edge_q_sib2
 * directly. edge_p / edge_p_back keep their internal Neighbor objects (so we
 * can later re-insert them) but are no longer reachable from the rest of the tree.
 */
struct BisectionState {
    PhyloNode* edge_p;        ///< One endpoint of the bisected edge.
    PhyloNode* edge_p_back;   ///< Other endpoint.
    PhyloNode* edge_p_sib1;
    PhyloNode* edge_p_sib2;
    PhyloNode* edge_q_sib1;
    PhyloNode* edge_q_sib2;
    bool valid;               ///< False if the requested edge is not bisectable.

    BisectionState() : edge_p(nullptr), edge_p_back(nullptr),
                       edge_p_sib1(nullptr), edge_p_sib2(nullptr),
                       edge_q_sib1(nullptr), edge_q_sib2(nullptr), valid(false) {}
};

/**
 * Bisect on the internal edge between edge_p and edge_p_back. Both endpoints
 * must be non-root, non-leaf, binary internal nodes that are currently
 * connected. Orientation-agnostic.
 * @return BisectionState with valid=true on success.
 */
BisectionState bisectInternalEdge(PhyloTree* tree, PhyloNode* edge_p,
                                  PhyloNode* edge_p_back);

/**
 * Re-insert edge_p on (edgeA_u, edgeA_v) in subtree A and edge_p_back on
 * (edgeB_u, edgeB_v) in subtree B, restoring the edge_p<->edge_p_back link.
 * Pass (s.edge_p_sib1, s.edge_p_sib2, s.edge_q_sib1, s.edge_q_sib2) for a
 * no-op round-trip.
 */
void reconnectBisected(PhyloTree* tree, const BisectionState& state,
                       PhyloNode* edgeA_u, PhyloNode* edgeA_v,
                       PhyloNode* edgeB_u, PhyloNode* edgeB_v);

/**
 * Inverse of reconnectBisected. The (edgeA_u, edgeA_v, edgeB_u, edgeB_v)
 * arguments MUST match the most recent reconnectBisected call.
 */
void undoReconnect(PhyloTree* tree, const BisectionState& state,
                   PhyloNode* edgeA_u, PhyloNode* edgeA_v,
                   PhyloNode* edgeB_u, PhyloNode* edgeB_v);

/**
 * Score a candidate TBR move via apply→Fitch→undo (correctness-first;
 * no closed-form delta yet — see plan/PLAN_TBR_OPTIMIZER.md).
 * @param pre_score  Score before the move (caller-supplied; not recomputed).
 * @return post_score - pre_score, or INT_MAX if the bisection is invalid.
 */
int evaluateTBRMove(PhyloTree* tree, int pre_score,
                    PhyloNode* edge_p, PhyloNode* edge_p_back,
                    PhyloNode* edgeA_u, PhyloNode* edgeA_v,
                    PhyloNode* edgeB_u, PhyloNode* edgeB_v);

/** Self-contained TBR primitive unit tests. @return # failed checks. */
int runTBRUnitTests(PhyloTree* tree);

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

#endif // TBR_PRIMITIVES_H
