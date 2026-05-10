#include "sproptimize.h"
#include "phylotree.h"
#include "phylonode.h"
#include "spr_utils.h"
#include "tbr_primitives.h"
#include <iostream>
#include <vector>
#include <climits>
#include <queue>
#include <chrono>
#include <utility>

using namespace std;
using namespace std::chrono;

/**
 * Find any internal edge (u, v) where both endpoints are non-root, non-leaf,
 * binary internal nodes. With require_internal_siblings, both subtrees must
 * also have ≥1 non-leaf sibling (so reconnect-elsewhere has alternatives).
 */
static bool findAnyInternalBinaryEdgePair(PhyloTree* tree,
                                          PhyloNode** out_u, PhyloNode** out_v,
                                          bool require_internal_siblings = false) {
    if (out_u) *out_u = NULL;
    if (out_v) *out_v = NULL;
    if (!tree || !tree->root) return false;
    int max_id = tree->fitchMaxNodeId();
    vector<PhyloNode*> stack;
    vector<bool> visited(max_id + 1, false);
    PhyloNode* start = (PhyloNode*)tree->root->neighbors[0]->node;
    stack.push_back(start);
    while (!stack.empty()) {
        PhyloNode* node = stack.back(); stack.pop_back();
        if (node->id < 0 || node->id > max_id) continue;
        if (visited[node->id]) continue;
        visited[node->id] = true;
        if (!node->isLeaf() && node->degree() == BINARY_NODE_DEGREE) {
            PhyloNode* parent = (PhyloNode*)node->neighbors[0]->node;
            if (parent && parent != tree->root && !parent->isLeaf()
                && parent->degree() == BINARY_NODE_DEGREE) {
                bool ok = !require_internal_siblings;
                if (require_internal_siblings) {
                    bool p_side_has_internal = false;
                    FOR_NEIGHBOR_IT(node, parent, it) {
                        PhyloNode* sib = (PhyloNode*)(*it)->node;
                        if (sib && !sib->isLeaf()) { p_side_has_internal = true; break; }
                    }
                    bool q_side_has_internal = false;
                    FOR_NEIGHBOR_IT(parent, node, it) {
                        PhyloNode* sib = (PhyloNode*)(*it)->node;
                        if (sib && !sib->isLeaf()) { q_side_has_internal = true; break; }
                    }
                    ok = p_side_has_internal && q_side_has_internal;
                }
                if (ok) {
                    if (out_u) *out_u = node;
                    if (out_v) *out_v = parent;
                    return true;
                }
            }
        }
        FOR_NEIGHBOR_IT(node, NULL, it) {
            PhyloNode* nb = (PhyloNode*)(*it)->node;
            if (nb && nb != tree->root && nb->id >= 0 && nb->id <= max_id
                && !visited[nb->id]) {
                stack.push_back(nb);
            }
        }
    }
    return false;
}

/**
 * Test-only wrapper returning just one endpoint. Not safe across topology
 * mutations because the partner is recovered via neighbors[0] later — use the
 * *Pair helper to capture both endpoints atomically.
 */
static PhyloNode* findAnyInternalBinaryEdge(PhyloTree* tree,
                                             bool require_internal_siblings = false) {
    PhyloNode* u = NULL; PhyloNode* v = NULL;
    findAnyInternalBinaryEdgePair(tree, &u, &v, require_internal_siblings);
    return u;
}

/** Recover edge_p_back from a known edge_p without relying on neighbors[0]. */
static PhyloNode* findInternalPartner(PhyloTree* tree, PhyloNode* u) {
    if (!u || u == tree->root || u->isLeaf()) return NULL;
    if (u->degree() != BINARY_NODE_DEGREE) return NULL;
    FOR_NEIGHBOR_IT(u, NULL, it) {
        PhyloNode* nb = (PhyloNode*)(*it)->node;
        if (nb && nb != tree->root && !nb->isLeaf()
            && nb->degree() == BINARY_NODE_DEGREE) {
            return nb;
        }
    }
    return NULL;
}

BisectionState bisectInternalEdge(PhyloTree* tree, PhyloNode* edge_p,
                                  PhyloNode* edge_p_back) {
    BisectionState s;  // s.valid = false by default

    // Both endpoints must be non-root, non-leaf, binary internal nodes.
    if (!tree || !edge_p || edge_p == tree->root) return s;
    if (edge_p->isLeaf()) return s;
    if (edge_p->degree() != BINARY_NODE_DEGREE) return s;
    if (!edge_p_back || edge_p_back == tree->root || edge_p_back == edge_p) return s;
    if (edge_p_back->isLeaf()) return s;
    if (edge_p_back->degree() != BINARY_NODE_DEGREE) return s;

    // Endpoints must currently share an edge.
    bool connected = false;
    FOR_NEIGHBOR_IT(edge_p, NULL, it) {
        if ((PhyloNode*)(*it)->node == edge_p_back) { connected = true; break; }
    }
    if (!connected) return s;

    // Capture the two siblings on each side (neighbors other than the partner).
    PhyloNode* edge_p_sib1 = NULL; PhyloNode* edge_p_sib2 = NULL;
    FOR_NEIGHBOR_IT(edge_p, edge_p_back, it) {
        PhyloNode* nb = (PhyloNode*)(*it)->node;
        if (!edge_p_sib1) edge_p_sib1 = nb;
        else if (!edge_p_sib2) edge_p_sib2 = nb;
    }
    PhyloNode* edge_q_sib1 = NULL; PhyloNode* edge_q_sib2 = NULL;
    FOR_NEIGHBOR_IT(edge_p_back, edge_p, it) {
        PhyloNode* nb = (PhyloNode*)(*it)->node;
        if (!edge_q_sib1) edge_q_sib1 = nb;
        else if (!edge_q_sib2) edge_q_sib2 = nb;
    }
    if (!edge_p_sib1 || !edge_p_sib2 || !edge_q_sib1 || !edge_q_sib2) return s;

    // Short-circuit each side. We only rewire the sibling-side pointers; the
    // Neighbor objects on edge_p / edge_p_back themselves stay untouched so
    // we can reverse the bisection later.
    edge_p_sib1->updateNeighbor(edge_p, edge_p_sib2);
    edge_p_sib2->updateNeighbor(edge_p, edge_p_sib1);
    edge_q_sib1->updateNeighbor(edge_p_back, edge_q_sib2);
    edge_q_sib2->updateNeighbor(edge_p_back, edge_q_sib1);

    s.edge_p = edge_p;
    s.edge_p_back = edge_p_back;
    s.edge_p_sib1 = edge_p_sib1; s.edge_p_sib2 = edge_p_sib2;
    s.edge_q_sib1 = edge_q_sib1; s.edge_q_sib2 = edge_q_sib2;
    s.valid = true;
    return s;
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
                       PhyloNode* edgeA_u, PhyloNode* edgeA_v,
                       PhyloNode* edgeB_u, PhyloNode* edgeB_v) {
    if (!s.valid) return;
    spliceFreeNodeOnEdge(edgeA_u, edgeA_v, s.edge_p,      s.edge_p_sib1, s.edge_p_sib2);
    spliceFreeNodeOnEdge(edgeB_u, edgeB_v, s.edge_p_back, s.edge_q_sib1, s.edge_q_sib2);
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
                   PhyloNode* edgeA_u, PhyloNode* edgeA_v,
                   PhyloNode* edgeB_u, PhyloNode* edgeB_v) {
    if (!s.valid) return;
    unspliceFreeNode(edgeA_u, edgeA_v, s.edge_p,      s.edge_p_sib1, s.edge_p_sib2);
    unspliceFreeNode(edgeB_u, edgeB_v, s.edge_p_back, s.edge_q_sib1, s.edge_q_sib2);
}

int evaluateTBRMove(PhyloTree* tree, int pre_score,
                    PhyloNode* edge_p, PhyloNode* edge_p_back,
                    PhyloNode* edgeA_u, PhyloNode* edgeA_v,
                    PhyloNode* edgeB_u, PhyloNode* edgeB_v) {
    BisectionState s = bisectInternalEdge(tree, edge_p, edge_p_back);
    if (!s.valid) return INT_MAX;

    reconnectBisected(tree, s, edgeA_u, edgeA_v, edgeB_u, edgeB_v);
    int post_score = tree->fitchRecomputeScore();

    // Restore: undo the alt reconnect, then reconnect at the original siblings.
    undoReconnect(tree, s, edgeA_u, edgeA_v, edgeB_u, edgeB_v);
    reconnectBisected(tree, s, s.edge_p_sib1, s.edge_p_sib2, s.edge_q_sib1, s.edge_q_sib2);
    int restored = tree->fitchRecomputeScore();

    if (restored != pre_score) {
        cerr << "[TBR] WARN: evaluateTBRMove restore mismatch: pre="
             << pre_score << " restored=" << restored << endl;
    }

    return post_score - pre_score;
}

// ----- test harness ----------------------------------------------------------

static int tbr_test_pass = 0;
static int tbr_test_fail = 0;

#define TBR_CHECK(cond, msg) do { \
    if (cond) { tbr_test_pass++; cout << "  PASS: " << msg << endl; } \
    else { tbr_test_fail++; cerr << "  FAIL: " << msg << endl; } \
} while (0)

#define TBR_TEST(name) cout << "\n--- TBR test: " << name << " ---" << endl

static void testBisectionShortCircuit(PhyloTree* tree) {
    TBR_TEST("bisectInternalEdge collapses siblings and unhooks edge_p<->edge_p_back");
    PhyloNode* edge_p = findAnyInternalBinaryEdge(tree);
    if (!edge_p) {
        cout << "  SKIP: no internal binary edge available" << endl;
        return;
    }

    BisectionState s = bisectInternalEdge(tree, edge_p, findInternalPartner(tree, edge_p));
    TBR_CHECK(s.valid, "bisection returned valid=true");
    if (!s.valid) return;

    TBR_CHECK(s.edge_p == edge_p, "edge_p preserved in state");
    TBR_CHECK(s.edge_p_back != nullptr && s.edge_p_back != edge_p,
              "edge_p_back is the other endpoint");
    TBR_CHECK(s.edge_p_sib1 && s.edge_p_sib2 && s.edge_q_sib1 && s.edge_q_sib2,
              "all four siblings recorded");

    bool edge_q_sib1_still_points_to_back = false;
    FOR_NEIGHBOR_IT(s.edge_q_sib1, NULL, it) {
        if ((PhyloNode*)(*it)->node == s.edge_p_back) {
            edge_q_sib1_still_points_to_back = true; break;
        }
    }
    TBR_CHECK(!edge_q_sib1_still_points_to_back,
              "edge_q_sib1 no longer points to edge_p_back");

    bool edge_q_sib1_now_points_to_edge_q_sib2 = false;
    FOR_NEIGHBOR_IT(s.edge_q_sib1, NULL, it) {
        if ((PhyloNode*)(*it)->node == s.edge_q_sib2) {
            edge_q_sib1_now_points_to_edge_q_sib2 = true; break;
        }
    }
    TBR_CHECK(edge_q_sib1_now_points_to_edge_q_sib2,
              "edge_q_sib1 now connects directly to edge_q_sib2");

    // Restore tree to its original state for subsequent tests by reconnecting
    // at the same siblings.
    reconnectBisected(tree, s, s.edge_p_sib1, s.edge_p_sib2, s.edge_q_sib1, s.edge_q_sib2);
}

// Bisect then reconnect at originals; Fitch score must be unchanged.
static void testBisectReconnectRoundTrip(PhyloTree* tree) {
    TBR_TEST("bisect -> reconnect at original siblings preserves Fitch score");
    int before = tree->fitchRunForSPR();

    PhyloNode* edge_p = findAnyInternalBinaryEdge(tree);
    if (!edge_p) { cout << "  SKIP: no internal binary edge" << endl; return; }

    BisectionState s = bisectInternalEdge(tree, edge_p, findInternalPartner(tree, edge_p));
    if (!s.valid) { cout << "  SKIP: bisection rejected" << endl; return; }

    reconnectBisected(tree, s, s.edge_p_sib1, s.edge_p_sib2, s.edge_q_sib1, s.edge_q_sib2);

    int after = tree->fitchRecomputeScore();
    TBR_CHECK(after == before,
              std::string("Fitch score unchanged after round-trip (before=")
              + std::to_string(before) + " after=" + std::to_string(after) + ")");
}

// bisect → reconnect → undoReconnect → reconnect; Fitch score must be unchanged.
static void testUndoReconnectInverse(PhyloTree* tree) {
    TBR_TEST("undoReconnect is the inverse of reconnectBisected (round-trip)");
    int before = tree->fitchRunForSPR();

    PhyloNode* edge_p = findAnyInternalBinaryEdge(tree);
    if (!edge_p) { cout << "  SKIP: no internal binary edge" << endl; return; }

    BisectionState s = bisectInternalEdge(tree, edge_p, findInternalPartner(tree, edge_p));
    if (!s.valid) { cout << "  SKIP" << endl; return; }

    reconnectBisected(tree, s, s.edge_p_sib1, s.edge_p_sib2, s.edge_q_sib1, s.edge_q_sib2);
    undoReconnect(tree, s, s.edge_p_sib1, s.edge_p_sib2, s.edge_q_sib1, s.edge_q_sib2);
    reconnectBisected(tree, s, s.edge_p_sib1, s.edge_p_sib2, s.edge_q_sib1, s.edge_q_sib2);

    int after = tree->fitchRecomputeScore();
    TBR_CHECK(after == before,
              std::string("Fitch score unchanged after bisect/reconn/undo/reconn (before=")
              + std::to_string(before) + " after=" + std::to_string(after) + ")");
}

/**
 * Find an alt edge in the subtree on one side of the bisection scar (anything
 * other than the trivial (inside, block) edge). Returns false on cherry
 * subtrees with no alternatives.
 */
static bool findAltEdgeInSubtree(PhyloNode* inside, PhyloNode* block,
                                  PhyloNode** out_u, PhyloNode** out_v) {
    if (out_u) *out_u = NULL;
    if (out_v) *out_v = NULL;
    FOR_NEIGHBOR_IT(inside, NULL, it) {
        PhyloNode* nb = (PhyloNode*)(*it)->node;
        if (nb && nb != block) {
            if (out_u) *out_u = inside;
            if (out_v) *out_v = nb;
            return true;
        }
    }
    if (block) {
        FOR_NEIGHBOR_IT(block, NULL, it) {
            PhyloNode* nb = (PhyloNode*)(*it)->node;
            if (nb && nb != inside) {
                if (out_u) *out_u = block;
                if (out_v) *out_v = nb;
                return true;
            }
        }
    }
    return false;
}

// Stronger inverse property: reconnect at a NON-original location, undo, then
// reconnect at original. Catches a class of bugs where undo only "appears" to
// restore but leaves stale neighbor pointers behind.
static void testBisectElsewhereThenUndo(PhyloTree* tree) {
    TBR_TEST("bisect -> reconnect elsewhere -> undo -> reconnect original preserves Fitch");
    int before = tree->fitchRunForSPR();

    PhyloNode* edge_p = findAnyInternalBinaryEdge(tree, /*require_internal_siblings=*/true);
    if (!edge_p) { cout << "  SKIP: no internal binary edge with internal siblings" << endl; return; }

    BisectionState s = bisectInternalEdge(tree, edge_p, findInternalPartner(tree, edge_p));
    if (!s.valid) { cout << "  SKIP" << endl; return; }

    // Find alternative edges in subtree A and subtree B.
    PhyloNode *aA_u = NULL, *aA_v = NULL, *aB_u = NULL, *aB_v = NULL;
    bool ok_a = findAltEdgeInSubtree(s.edge_p_sib1, s.edge_p_sib2, &aA_u, &aA_v);
    bool ok_b = findAltEdgeInSubtree(s.edge_q_sib1, s.edge_q_sib2, &aB_u, &aB_v);
    if (!ok_a || !ok_b) {
        reconnectBisected(tree, s, s.edge_p_sib1, s.edge_p_sib2, s.edge_q_sib1, s.edge_q_sib2);
        cout << "  SKIP: no alternative edges reachable" << endl;
        return;
    }

    reconnectBisected(tree, s, aA_u, aA_v, aB_u, aB_v);
    undoReconnect(tree, s, aA_u, aA_v, aB_u, aB_v);
    reconnectBisected(tree, s, s.edge_p_sib1, s.edge_p_sib2, s.edge_q_sib1, s.edge_q_sib2);

    int after = tree->fitchRecomputeScore();
    TBR_CHECK(after == before,
              std::string("Fitch score restored after bisect/reconn-elsewhere/undo/reconn-original (before=")
              + std::to_string(before) + " after=" + std::to_string(after) + ")");
}

// Up to N internal binary edges whose endpoints have non-leaf siblings on
// both sides — i.e. usable as TBR test candidates with at least one alt edge
// available in each subtree.
static vector<PhyloNode*> collectInternalBinaryEdges(PhyloTree* tree, int max_count) {
    vector<PhyloNode*> result;
    if (!tree || !tree->root) return result;
    int max_id = tree->fitchMaxNodeId();
    vector<bool> visited(max_id + 1, false);
    vector<PhyloNode*> stack;
    PhyloNode* start = (PhyloNode*)tree->root->neighbors[0]->node;
    stack.push_back(start);
    while (!stack.empty() && (int)result.size() < max_count) {
        PhyloNode* node = stack.back(); stack.pop_back();
        if (node->id < 0 || node->id > max_id) continue;
        if (visited[node->id]) continue;
        visited[node->id] = true;
        if (!node->isLeaf() && node->degree() == BINARY_NODE_DEGREE) {
            PhyloNode* parent = (PhyloNode*)node->neighbors[0]->node;
            if (parent && parent != tree->root && !parent->isLeaf()
                && parent->degree() == BINARY_NODE_DEGREE) {
                bool p_internal = false, q_internal = false;
                FOR_NEIGHBOR_IT(node, parent, it) {
                    if (!((PhyloNode*)(*it)->node)->isLeaf()) { p_internal = true; break; }
                }
                FOR_NEIGHBOR_IT(parent, node, it) {
                    if (!((PhyloNode*)(*it)->node)->isLeaf()) { q_internal = true; break; }
                }
                if (p_internal && q_internal) {
                    result.push_back(node);
                }
            }
        }
        FOR_NEIGHBOR_IT(node, NULL, it) {
            PhyloNode* nb = (PhyloNode*)(*it)->node;
            if (nb && nb != tree->root && nb->id >= 0 && nb->id <= max_id
                && !visited[nb->id]) {
                stack.push_back(nb);
            }
        }
    }
    return result;
}

// Per-candidate round trip: isolates primitive bugs from evaluateTBRMove bugs.
static void testRoundTripPerCandidate(PhyloTree* tree) {
    TBR_TEST("bisect+reconnect-original round trip preserves Fitch (per candidate)");
    int pre = tree->fitchRecomputeScore();
    vector<PhyloNode*> candidates = collectInternalBinaryEdges(tree, 10);
    int failed = 0;
    for (size_t i = 0; i < candidates.size(); i++) {
        PhyloNode* ep = candidates[i];
        PhyloNode* partner = findInternalPartner(tree, ep);
        if (!partner) continue;
        BisectionState s = bisectInternalEdge(tree, ep, partner);
        if (!s.valid) continue;
        reconnectBisected(tree, s, s.edge_p_sib1, s.edge_p_sib2, s.edge_q_sib1, s.edge_q_sib2);
        int got = tree->fitchRecomputeScore();
        if (got != pre) {
            cerr << "  candidate[" << i << "] round-trip mismatch: got " << got
                 << " expected " << pre << " (ep=" << ep->id
                 << " partner=" << partner->id << ")" << endl;
            failed++;
        }
    }
    TBR_CHECK(failed == 0, std::string("Bare round-trip preserves Fitch on ")
              + std::to_string(candidates.size()) + " candidates ("
              + std::to_string(failed) + " failed)");
}

// Central correctness gate: evaluateTBRMove's delta == direct apply-and-recompute.
static void testTBRDeltaCorrectness(PhyloTree* tree) {
    TBR_TEST("evaluateTBRMove delta == direct apply-and-recompute (up to 10 cases)");
    int pre = tree->fitchRecomputeScore();

    vector<PhyloNode*> candidates = collectInternalBinaryEdges(tree, 10);
    if (candidates.empty()) {
        cout << "  SKIP: no internal binary edges with internal siblings" << endl;
        return;
    }

    int checked = 0, mismatches = 0;
    for (PhyloNode* edge_p : candidates) {
        PhyloNode* partner = findInternalPartner(tree, edge_p);
        if (!partner) continue;

        // Probe-bisect to discover the sibling layout, then restore.
        BisectionState probe = bisectInternalEdge(tree, edge_p, partner);
        if (!probe.valid) continue;
        PhyloNode *aA_u = NULL, *aA_v = NULL, *aB_u = NULL, *aB_v = NULL;
        bool ok_a = findAltEdgeInSubtree(probe.edge_p_sib1, probe.edge_p_sib2, &aA_u, &aA_v);
        bool ok_b = findAltEdgeInSubtree(probe.edge_q_sib1, probe.edge_q_sib2, &aB_u, &aB_v);
        if (!ok_a) { aA_u = probe.edge_p_sib1; aA_v = probe.edge_p_sib2; }
        if (!ok_b) { aB_u = probe.edge_q_sib1; aB_v = probe.edge_q_sib2; }
        reconnectBisected(tree, probe, probe.edge_p_sib1, probe.edge_p_sib2,
                          probe.edge_q_sib1, probe.edge_q_sib2);
        int restored_pre = tree->fitchRecomputeScore();
        if (restored_pre != pre) {
            cerr << "  setup-mismatch on candidate " << checked
                 << " (got " << restored_pre << ", expected " << pre << ")" << endl;
            mismatches++; checked++; continue;
        }

        // Ground truth: bisect, reconnect, recompute, restore.
        BisectionState s = bisectInternalEdge(tree, edge_p, partner);
        reconnectBisected(tree, s, aA_u, aA_v, aB_u, aB_v);
        int gt_post = tree->fitchRecomputeScore();
        int gt_delta = gt_post - pre;
        undoReconnect(tree, s, aA_u, aA_v, aB_u, aB_v);
        reconnectBisected(tree, s, s.edge_p_sib1, s.edge_p_sib2, s.edge_q_sib1, s.edge_q_sib2);
        int gt_restored = tree->fitchRecomputeScore();
        if (gt_restored != pre) {
            cerr << "  GT-restore mismatch on candidate " << checked << endl;
            mismatches++; checked++; continue;
        }

        int reported_delta = evaluateTBRMove(tree, pre, edge_p, partner,
                                              aA_u, aA_v, aB_u, aB_v);

        if (reported_delta != gt_delta) {
            cerr << "  delta mismatch on candidate " << checked
                 << ": reported=" << reported_delta << " gt=" << gt_delta << endl;
            mismatches++;
        }
        int final_score = tree->fitchRecomputeScore();
        if (final_score != pre) {
            cerr << "  evaluateTBRMove did not restore tree on candidate "
                 << checked << ": " << final_score << " vs " << pre << endl;
            mismatches++;
        }
        checked++;
    }

    TBR_CHECK(mismatches == 0,
              std::string("All ") + std::to_string(checked)
              + " TBR delta computations match ground truth ("
              + std::to_string(mismatches) + " failed)");
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
        PhyloNode* n = stack.back(); stack.pop_back();
        if (n->id < 0 || n->id > max_id) continue;
        if (visited[n->id]) continue;
        visited[n->id] = true;
        if (!n->isLeaf() && n->degree() == BINARY_NODE_DEGREE && n != tree->root) {
            FOR_NEIGHBOR_IT(n, NULL, it) {
                PhyloNode* m = (PhyloNode*)(*it)->node;
                if (!m || m == tree->root || m->isLeaf()
                    || m->degree() != BINARY_NODE_DEGREE) continue;
                if (m->id < 0 || m->id > max_id) continue;
                int lo = min(n->id, m->id), hi = max(n->id, m->id);
                if ((int)edge_seen[lo].size() <= hi) edge_seen[lo].resize(hi + 1, false);
                if (edge_seen[lo][hi]) continue;
                edge_seen[lo][hi] = true;
                out.push_back({n, m});
            }
        }
        FOR_NEIGHBOR_IT(n, NULL, it) {
            PhyloNode* nb = (PhyloNode*)(*it)->node;
            if (nb && nb->id >= 0 && nb->id <= max_id && !visited[nb->id]) {
                stack.push_back(nb);
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
    int initial = cur_score;

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
        PhyloNode *best_ep = NULL, *best_ebp = NULL;
        PhyloNode *best_aA_u = NULL, *best_aA_v = NULL;
        PhyloNode *best_aB_u = NULL, *best_aB_v = NULL;
        long candidates_evaluated = 0;

        for (auto& edge : internal_edges) {
            if (wall_exceeded()) break;
            PhyloNode* edge_p = edge.first;
            PhyloNode* edge_p_back = edge.second;

            BisectionState s = bisectInternalEdge(tree, edge_p, edge_p_back);
            if (!s.valid) continue;

            vector<pair<PhyloNode*, PhyloNode*>> alt_a, alt_b;
            enumerateEdgesInSubtree(s.edge_p_sib1, s.edge_p_sib2, tbr_max_radius, alt_a);
            enumerateEdgesInSubtree(s.edge_q_sib1, s.edge_q_sib2, tbr_max_radius, alt_b);

            // Reconnect at originals so the tree is whole again before each
            // candidate's own bisect/reconnect/undo cycle.
            reconnectBisected(tree, s, s.edge_p_sib1, s.edge_p_sib2, s.edge_q_sib1, s.edge_q_sib2);

            for (auto& ea : alt_a) {
                if (wall_exceeded()) break;
                for (auto& eb : alt_b) {
                    if (wall_exceeded()) break;
                    int delta = evaluateTBRMove(tree, cur_score,
                                                 edge_p, edge_p_back,
                                                 ea.first, ea.second,
                                                 eb.first, eb.second);
                    candidates_evaluated++;
                    if (delta < best_delta) {
                        best_delta = delta;
                        best_ep = edge_p; best_ebp = edge_p_back;
                        best_aA_u = ea.first; best_aA_v = ea.second;
                        best_aB_u = eb.first; best_aB_v = eb.second;
                    }
                }
            }
        }

        if (best_delta < 0 && best_ep) {
            // Apply the best improving move permanently.
            BisectionState s = bisectInternalEdge(tree, best_ep, best_ebp);
            if (s.valid) {
                reconnectBisected(tree, s, best_aA_u, best_aA_v,
                                  best_aB_u, best_aB_v);
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

    cout << "TBR complete: " << initial << " -> " << cur_score
         << " (total delta=" << (initial - cur_score) << ", "
         << round << " rounds)" << endl;
    return cur_score;
}

int runTBRUnitTests(PhyloTree* tree) {
    cout << "\n========================================" << endl;
    cout << "  TBR Primitives Unit Tests" << endl;
    cout << "========================================" << endl;
    tbr_test_pass = 0; tbr_test_fail = 0;

    // Each test self-restores the tree to its original topology before
    // returning, so the order is independent.
    testBisectionShortCircuit(tree);
    testBisectReconnectRoundTrip(tree);
    testUndoReconnectInverse(tree);
    testBisectElsewhereThenUndo(tree);
    testRoundTripPerCandidate(tree);
    testTBRDeltaCorrectness(tree);

    cout << "\nTBR Tests: " << tbr_test_pass << " passed, "
         << tbr_test_fail << " failed" << endl;
    return tbr_test_fail;
}
