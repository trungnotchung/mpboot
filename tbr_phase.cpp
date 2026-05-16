#include "optimizer.h"
#include "phylotree.h"
#include "phylonode.h"
#include "spr_utils.h"
#include "spr_phase.h"
#include "spr_delta_exact.h"
#include "tbr_phase.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <climits>
#include <queue>
#include <chrono>
#include <utility>

using namespace std;
using namespace std::chrono;

BisectionState bisectInternalEdge(PhyloTree* tree, PhyloNode* bisect_a,
                                  PhyloNode* bisect_b) {
    BisectionState state;  // valid = false by default

    // Both endpoints must be non-root, non-leaf, binary internal nodes.
    if (!tree || !bisect_a || bisect_a == tree->root) return state;
    if (bisect_a->isLeaf()) return state;
    if (bisect_a->degree() != BINARY_NODE_DEGREE) return state;
    if (!bisect_b || bisect_b == tree->root || bisect_b == bisect_a) return state;
    if (bisect_b->isLeaf()) return state;
    if (bisect_b->degree() != BINARY_NODE_DEGREE) return state;

    // Endpoints must currently share an edge.
    bool connected = false;
    FOR_NEIGHBOR_IT(bisect_a, NULL, it) {
        if ((PhyloNode*)(*it)->node == bisect_b) { connected = true; break; }
    }
    if (!connected) return state;

    // Capture the two siblings on each side (neighbors other than the partner).
    PhyloNode* bisect_a_sib1 = NULL; PhyloNode* bisect_a_sib2 = NULL;
    FOR_NEIGHBOR_IT(bisect_a, bisect_b, it) {
        PhyloNode* neighbor = (PhyloNode*)(*it)->node;
        if (!bisect_a_sib1) bisect_a_sib1 = neighbor;
        else if (!bisect_a_sib2) bisect_a_sib2 = neighbor;
    }
    PhyloNode* bisect_b_sib1 = NULL; PhyloNode* bisect_b_sib2 = NULL;
    FOR_NEIGHBOR_IT(bisect_b, bisect_a, it) {
        PhyloNode* neighbor = (PhyloNode*)(*it)->node;
        if (!bisect_b_sib1) bisect_b_sib1 = neighbor;
        else if (!bisect_b_sib2) bisect_b_sib2 = neighbor;
    }
    if (!bisect_a_sib1 || !bisect_a_sib2 || !bisect_b_sib1 || !bisect_b_sib2) return state;

    // Short-circuit each side. We only rewire the sibling-side pointers; the
    // Neighbor objects on bisect_a / bisect_b themselves stay untouched so
    // we can reverse the bisection later.
    bisect_a_sib1->updateNeighbor(bisect_a, bisect_a_sib2);
    bisect_a_sib2->updateNeighbor(bisect_a, bisect_a_sib1);
    bisect_b_sib1->updateNeighbor(bisect_b, bisect_b_sib2);
    bisect_b_sib2->updateNeighbor(bisect_b, bisect_b_sib1);

    state.bisect_a = bisect_a;
    state.bisect_b = bisect_b;
    state.bisect_a_sib1 = bisect_a_sib1; state.bisect_a_sib2 = bisect_a_sib2;
    state.bisect_b_sib1 = bisect_b_sib1; state.bisect_b_sib2 = bisect_b_sib2;
    state.valid = true;
    return state;
}

/**
 * Subdivide edge (u, v) by inserting free_node between them. free_node must
 * still hold (old_sib1, old_sib2) in its non-back slots from the prior bisect.
 */
static void spliceFreeNodeOnEdge(PhyloNode* u, PhyloNode* v, PhyloNode* free_node,
                                  PhyloNode* old_sib1, PhyloNode* old_sib2) {
    u->updateNeighbor(v, free_node);
    v->updateNeighbor(u, free_node);
    free_node->updateNeighbor(old_sib1, u);
    free_node->updateNeighbor(old_sib2, v);
}

void reconnectBisected(PhyloTree* /*tree*/, const BisectionState& s,
                       PhyloNode* reattach_a_u, PhyloNode* reattach_a_v,
                       PhyloNode* reattach_b_u, PhyloNode* reattach_b_v) {
    if (!s.valid) return;
    spliceFreeNodeOnEdge(reattach_a_u, reattach_a_v, s.bisect_a,      s.bisect_a_sib1, s.bisect_a_sib2);
    spliceFreeNodeOnEdge(reattach_b_u, reattach_b_v, s.bisect_b, s.bisect_b_sib1, s.bisect_b_sib2);
}

/** Inverse of spliceFreeNodeOnEdge. */
static void unspliceFreeNode(PhyloNode* u, PhyloNode* v, PhyloNode* free_node,
                              PhyloNode* old_sib1, PhyloNode* old_sib2) {
    u->updateNeighbor(free_node, v);
    v->updateNeighbor(free_node, u);
    free_node->updateNeighbor(u, old_sib1);
    free_node->updateNeighbor(v, old_sib2);
}

void undoReconnect(PhyloTree* /*tree*/, const BisectionState& s,
                   PhyloNode* reattach_a_u, PhyloNode* reattach_a_v,
                   PhyloNode* reattach_b_u, PhyloNode* reattach_b_v) {
    if (!s.valid) return;
    unspliceFreeNode(reattach_a_u, reattach_a_v, s.bisect_a,      s.bisect_a_sib1, s.bisect_a_sib2);
    unspliceFreeNode(reattach_b_u, reattach_b_v, s.bisect_b, s.bisect_b_sib1, s.bisect_b_sib2);
}

int evaluateTBRMove(PhyloTree* tree, int pre_score,
                    PhyloNode* bisect_a, PhyloNode* bisect_b,
                    PhyloNode* reattach_a_u, PhyloNode* reattach_a_v,
                    PhyloNode* reattach_b_u, PhyloNode* reattach_b_v) {
    BisectionState s = bisectInternalEdge(tree, bisect_a, bisect_b);
    if (!s.valid) return INT_MAX;

    reconnectBisected(tree, s, reattach_a_u, reattach_a_v, reattach_b_u, reattach_b_v);
    int post_score = tree->fitchRecomputeScore();

    // Restore: undo the alt reconnect, then reconnect at the original siblings.
    undoReconnect(tree, s, reattach_a_u, reattach_a_v, reattach_b_u, reattach_b_v);
    reconnectBisected(tree, s, s.bisect_a_sib1, s.bisect_a_sib2, s.bisect_b_sib1, s.bisect_b_sib2);
    int restored = tree->fitchRecomputeScore();

    if (restored != pre_score) {
        cerr << "[TBR] WARN: evaluateTBRMove restore mismatch: pre="
             << pre_score << " restored=" << restored << endl;
    }

    return post_score - pre_score;
}


/**
 * Enumerate every internal binary edge (non-root, non-leaf endpoints on both
 * sides) in the tree, orientation-agnostic; each edge once. De-duplicated via
 * a sparse 2D bitset keyed on (min_id, max_id) of the endpoints.
 */
static vector<pair<PhyloNode*, PhyloNode*>>
enumerateAllInternalBinaryEdges(PhyloTree* tree) {
    vector<pair<PhyloNode*, PhyloNode*>> out;
    if (!tree || !tree->root) return out;
    int max_id = tree->fitchMaxNodeId();
    vector<bool> visited(max_id + 1, false);
    vector<vector<bool>> edge_seen(max_id + 1);  // sparse, lazy-allocated.

    vector<PhyloNode*> stack;
    PhyloNode* start = (PhyloNode*)tree->root->neighbors[0]->node;
    stack.push_back(start);
    while (!stack.empty()) {
        PhyloNode* node = stack.back(); stack.pop_back();
        if (node->id < 0 || node->id > max_id) continue;
        if (visited[node->id]) continue;
        visited[node->id] = true;
        if (!node->isLeaf() && node->degree() == BINARY_NODE_DEGREE && node != tree->root) {
            FOR_NEIGHBOR_IT(node, NULL, it) {
                PhyloNode* neighbor = (PhyloNode*)(*it)->node;
                if (!neighbor || neighbor == tree->root || neighbor->isLeaf()
                    || neighbor->degree() != BINARY_NODE_DEGREE) continue;
                if (neighbor->id < 0 || neighbor->id > max_id) continue;
                int lo = min(node->id, neighbor->id), hi = max(node->id, neighbor->id);
                if ((int)edge_seen[lo].size() <= hi) edge_seen[lo].resize(hi + 1, false);
                if (edge_seen[lo][hi]) continue;
                edge_seen[lo][hi] = true;
                out.push_back({node, neighbor});
            }
        }
        FOR_NEIGHBOR_IT(node, NULL, it) {
            PhyloNode* neighbor = (PhyloNode*)(*it)->node;
            if (neighbor && neighbor->id >= 0 && neighbor->id <= max_id && !visited[neighbor->id]) {
                stack.push_back(neighbor);
            }
        }
    }
    return out;
}

/**
 * BFS-enumerate edges within `max_depth` hops of the bisection scar in one
 * subtree. The (seed1, seed2) trivial edge — reconnect-at-original — is
 * excluded by blocking each seed from crossing into the other.
 */
static void enumerateEdgesInSubtree(PhyloNode* seed1, PhyloNode* seed2,
                                     int max_depth,
                                     vector<pair<PhyloNode*, PhyloNode*>>& out) {
    if (max_depth <= 0) return;
    queue<tuple<PhyloNode*, PhyloNode*, int>> q;
    q.push(make_tuple(seed1, seed2, 0));
    q.push(make_tuple(seed2, seed1, 0));
    while (!q.empty()) {
        PhyloNode* node = get<0>(q.front());
        PhyloNode* parent = get<1>(q.front());
        int depth = get<2>(q.front());
        q.pop();
        if (depth >= max_depth) continue;
        FOR_NEIGHBOR_IT(node, parent, it) {
            PhyloNode* child = (PhyloNode*)(*it)->node;
            if (!child) continue;
            out.push_back({node, child});
            q.push(make_tuple(child, node, depth + 1));
        }
    }
}

int optimizeTBRAtRadius(PhyloTree* tree, int tbr_max_radius,
                        int known_score, double wall_seconds, int max_rounds) {
    int cur_score = (known_score > 0) ? known_score : tree->fitchRecomputeScore();
    int initial_score = cur_score;

    auto wall_start = high_resolution_clock::now();
    auto wall_exceeded = [&]() {
        if (wall_seconds <= 0.0) return false;
        double el = duration_cast<milliseconds>(
            high_resolution_clock::now() - wall_start).count() / 1000.0;
        return el > wall_seconds;
    };

    cout << "=== TBR optimize (tbr_max_radius=" << tbr_max_radius
         << ", start=" << cur_score << ") ===" << endl;

    int round = 0;
    bool improved = true;
    while (improved && round < max_rounds) {
        if (wall_exceeded()) {
            cout << "  TBR wall budget reached after " << round << " rounds" << endl;
            break;
        }
        improved = false;
        round++;

        vector<pair<PhyloNode*, PhyloNode*>> internal_edges =
            enumerateAllInternalBinaryEdges(tree);

        int best_delta = 0;
        PhyloNode *best_bisect_a = NULL, *best_bisect_b = NULL;
        PhyloNode *best_reattach_a_u = NULL, *best_reattach_a_v = NULL;
        PhyloNode *best_reattach_b_u = NULL, *best_reattach_b_v = NULL;
        long candidates_evaluated = 0;

        for (auto& edge : internal_edges) {
            if (wall_exceeded()) break;
            PhyloNode* bisect_a = edge.first;
            PhyloNode* bisect_b = edge.second;

            BisectionState s = bisectInternalEdge(tree, bisect_a, bisect_b);
            if (!s.valid) continue;

            vector<pair<PhyloNode*, PhyloNode*>> alt_a, alt_b;
            enumerateEdgesInSubtree(s.bisect_a_sib1, s.bisect_a_sib2, tbr_max_radius, alt_a);
            enumerateEdgesInSubtree(s.bisect_b_sib1, s.bisect_b_sib2, tbr_max_radius, alt_b);

            // Reconnect at originals so the tree is whole again before each
            // candidate's own bisect/reconnect/undo cycle.
            reconnectBisected(tree, s, s.bisect_a_sib1, s.bisect_a_sib2, s.bisect_b_sib1, s.bisect_b_sib2);

            for (auto& ea : alt_a) {
                if (wall_exceeded()) break;
                for (auto& eb : alt_b) {
                    if (wall_exceeded()) break;
                    int delta = evaluateTBRMove(tree, cur_score,
                                                 bisect_a, bisect_b,
                                                 ea.first, ea.second,
                                                 eb.first, eb.second);
                    candidates_evaluated++;
                    if (delta < best_delta) {
                        best_delta = delta;
                        best_bisect_a = bisect_a; best_bisect_b = bisect_b;
                        best_reattach_a_u = ea.first; best_reattach_a_v = ea.second;
                        best_reattach_b_u = eb.first; best_reattach_b_v = eb.second;
                    }
                }
            }
        }

        if (best_delta < 0 && best_bisect_a) {
            // Apply the best improving move permanently.
            BisectionState s = bisectInternalEdge(tree, best_bisect_a, best_bisect_b);
            if (s.valid) {
                reconnectBisected(tree, s, best_reattach_a_u, best_reattach_a_v,
                                  best_reattach_b_u, best_reattach_b_v);
                int new_score = tree->fitchRecomputeScore();
                cout << "  TBR round " << round << ": " << cur_score << " -> "
                     << new_score << " (best_delta=" << best_delta
                     << ", " << candidates_evaluated << " candidates eval'd)"
                     << endl;
                cur_score = new_score;
                improved = true;
            }
        } else {
            cout << "  TBR round " << round << ": no improving move ("
                 << candidates_evaluated << " candidates eval'd)" << endl;
        }
    }

    cout << "TBR complete: " << initial_score << " -> " << cur_score
         << " (total delta=" << (initial_score - cur_score) << ", "
         << round << " rounds)" << endl;
    return cur_score;
}

static const int TBR_PHASE_RADIUS_STALL_LIMIT    = 2;
static const int TBR_PHASE_INTERLEAVE_SPR_PASSES = 10;

int runTBRPhase(PlacementOptimizer* optimizer,
                PhyloTree* tree,
                const PlacementOptimizeOptions& opts,
                int entry_score,
                int max_id,
                high_resolution_clock::time_point start_time) {
    int tracked_score = entry_score;
    int max_radius = opts.max_radius;
    double wall_seconds = opts.wall_seconds;

    auto optimizer_total_elapsed = [&]() {
        return duration_cast<milliseconds>(
            high_resolution_clock::now() - start_time).count() / 1000.0;
    };
    auto wall_exceeded = [&]() {
        return wall_seconds > 0.0 && optimizer_total_elapsed() > wall_seconds;
    };

    cout << "\n=== TBR phase (iters=" << opts.tbr_iters
         << ", tbr_max_radius=" << opts.tbr_max_radius
         << ", interleaved with SPR) ===" << endl;

    for (int it = 0; it < opts.tbr_iters; it++) {
        if (wall_exceeded()) break;
        int pre_iter = tracked_score;

        double tbr_remaining = (wall_seconds > 0.0)
            ? std::max(0.5, wall_seconds - optimizer_total_elapsed()) : 0.0;
        int after_tbr = optimizeTBRAtRadius(tree, opts.tbr_max_radius,
                                             tracked_score, tbr_remaining,
                                             /*max_rounds=*/1);
        if (after_tbr >= pre_iter) {
            cout << "  TBR phase converged at " << pre_iter
                 << " (no improving TBR move past SPR-local optimum)" << endl;
            break;
        }
        tracked_score = after_tbr;

        orientTreeToRoot(tree, max_id);
        SPRDeltaExact::precomputeDepths(tree);

        auto spr_remaining = [&]() -> double {
            if (wall_seconds <= 0.0) return 0.0;
            double remaining_seconds = wall_seconds - optimizer_total_elapsed();
            return (remaining_seconds < 0.5) ? 0.5 : remaining_seconds;
        };
        for (int pass = 0; pass < TBR_PHASE_INTERLEAVE_SPR_PASSES; pass++) {
            if (wall_exceeded()) break;
            int spr_pre = tracked_score;
            int consecutive_empty = 0;
            int cur = spr_pre;
            for (int r = 1; ; r *= 2) {
                r = std::min(r, max_radius);
                int after = optimizer->optimizeAtRadius(r, cur, spr_remaining());
                if (after >= cur) {
                    consecutive_empty++;
                    if (consecutive_empty >= TBR_PHASE_RADIUS_STALL_LIMIT && r < max_radius)
                        break;
                } else {
                    consecutive_empty = 0;
                    cur = after;
                }
                if (r == max_radius) break;
                if (wall_exceeded()) break;
            }
            tracked_score = cur;
            if (cur >= spr_pre) break;
        }
        cout << "  TBR+SPR iter " << (it + 1) << ": " << pre_iter
             << " -> " << tracked_score << " (delta=" << (pre_iter - tracked_score)
             << ")" << endl;
    }

    return tracked_score;
}
