#include "sproptimize.h"
#include "fitch.h"
#include "spr_delta_exact.h"
#include "spr_mutation_ops.h"
#include "phylotree.h"
#include "phylonode.h"
#include <iostream>
#include <queue>
#include <set>
#include <chrono>
#include <algorithm>
#include <cassert>
#include <functional>

using namespace std;
using namespace std::chrono;

// Radii to search in order; 0 = unbounded (explore entire tree).
static const int SPR_RADII[] = {1, 2, 4, 8, 16, 0};
static const int NUM_SPR_RADII = sizeof(SPR_RADII) / sizeof(SPR_RADII[0]);
static const int DEFAULT_MAX_PASSES = 1;
static const int MAX_ROUNDS_PER_RADIUS = 100;
static const int BINARY_NODE_DEGREE = 3;

static inline nuc_one_hot fitchMerge(nuc_one_hot a, nuc_one_hot b) {
    nuc_one_hot inter = a & b;
    return inter ? inter : (a | b);
}

static inline int penaltyDelta(nuc_one_hot old_child, nuc_one_hot new_child,
                                nuc_one_hot other, nuc_one_hot& new_node_out) {
    int old_pen = (old_child & other) ? 0 : 1;
    nuc_one_hot inter = new_child & other;
    new_node_out = inter ? inter : (new_child | other);
    int new_pen = inter ? 0 : 1;
    return new_pen - old_pen;
}

static inline int rootEdgeDelta(nuc_one_hot root_f, nuc_one_hot old_rn, nuc_one_hot new_rn) {
    int old_re = (root_f & old_rn) ? 0 : 1;
    int new_re = (root_f & new_rn) ? 0 : 1;
    return new_re - old_re;
}

static PhyloNode* findOtherChild(PhyloNode* node, PhyloNode* parent,
                                  PhyloNode* known_child) {
    if (!node) return nullptr;
    FOR_NEIGHBOR_IT(node, parent, nit) {
        PhyloNode* c = (PhyloNode*)(*nit)->node;
        if (c != known_child) return c;
    }
    return nullptr;
}

static inline void collectDiffs(Fitch* cf, PhyloNode* node, vector<int>& out) {
    const vector<int>* d = cf->getFitchDiffs(node);
    if (d && !d->empty()) out.insert(out.end(), d->begin(), d->end());
}

class SPRSourceState {
public:
    SPRSourceState(Fitch* cf, PhyloNode* src, PhyloNode* src_parent);

    int evaluate(PhyloNode* dst, PhyloNode* dst_parent) const;
    bool isValid() const { return valid; }

private:
    Fitch* cf;
    PhyloNode* src;
    PhyloNode* sp;  // = src_parent
    PhyloNode* gp;
    PhyloNode* sib;
    bool valid;

    int nptn;

    const nuc_one_hot* src_arr;
    const nuc_one_hot* sib_arr;
    const nuc_one_hot* sp_arr;
    const nuc_one_hot* gp_arr;

    PhyloNode* gp_other;
    const nuc_one_hot* gp_other_arr;

    struct PathNode {
        PhyloNode* node;
        const nuc_one_hot* other_arr;
        const nuc_one_hot* node_arr;
    };
    vector<PathNode> src_to_root;
    vector<PathNode> above_gp;
    vector<int> base_diffs;

    int evalCaseA(PhyloNode* dst, PhyloNode* dst_parent) const;
    int evalCaseB1(PhyloNode* dst) const;
    int evalCaseB2(PhyloNode* dst, PhyloNode* dst_parent, PhyloNode* lca) const;

    static void propagate(const vector<PathNode>& path, int p, int freq,
                          nuc_one_hot& old_f, nuc_one_hot& new_f, int& score);
    void propagateRange(size_t start, size_t end, int p, int freq,
                        nuc_one_hot& old_f, nuc_one_hot& new_f, int& score) const;
};

void SPRSourceState::propagate(const vector<PathNode>& path, int p, int freq,
                                nuc_one_hot& old_f, nuc_one_hot& new_f, int& score) {
    for (size_t i = 0; i < path.size(); i++) {
        if (old_f == new_f) break;
        const auto& pn = path[i];
        nuc_one_hot new_node_f;
        score += penaltyDelta(old_f, new_f, pn.other_arr[p], new_node_f) * freq;
        old_f = pn.node_arr[p];
        new_f = new_node_f;
    }
}

void SPRSourceState::propagateRange(size_t start, size_t end, int p, int freq,
                                     nuc_one_hot& old_f, nuc_one_hot& new_f, int& score) const {
    for (size_t i = start; i < end && i < src_to_root.size(); i++) {
        if (old_f == new_f) break;
        const auto& pn = src_to_root[i];
        nuc_one_hot new_node_f;
        score += penaltyDelta(old_f, new_f, pn.other_arr[p], new_node_f) * freq;
        old_f = pn.node_arr[p];
        new_f = new_node_f;
    }
}

SPRSourceState::SPRSourceState(Fitch* cf, PhyloNode* src, PhyloNode* src_parent)
    : cf(cf), src(src), sp(src_parent), valid(false) {
    nptn = cf->getNumPatterns();

    gp = SPRMutationOps::getParent(sp);
    sib = findOtherChild(sp, gp, src);
    if (!sib || !gp) return;
    valid = true;

    src_arr = cf->getMajorArrayForNode(src);
    sib_arr = cf->getMajorArrayForNode(sib);
    sp_arr  = cf->getMajorArrayForNode(sp);
    gp_arr  = cf->getMajorArrayForNode(gp);

    gp_other = findOtherChild(gp, SPRMutationOps::getParent(gp), sp);
    gp_other_arr = gp_other ? cf->getMajorArrayForNode(gp_other) : nullptr;

    // Build above-gp path (for Case A)
    {
        PhyloNode* prev = gp;
        PhyloNode* cur = SPRMutationOps::getParent(gp);
        while (cur) {
            PhyloNode* par = SPRMutationOps::getParent(cur);
            PhyloNode* other = findOtherChild(cur, par, prev);
            const nuc_one_hot* other_a = other
                ? cf->getMajorArrayForNode(other)
                : cf->getMajorArrayForNode(cur);
            above_gp.push_back({cur, other_a, cf->getMajorArrayForNode(cur)});
            prev = cur;
            cur = par;
        }
    }

    // Build src_to_root path (for Case B): gp → parent(gp) → ... → root
    {
        PhyloNode* prev = sp;
        PhyloNode* cur = gp;
        while (cur) {
            PhyloNode* par = SPRMutationOps::getParent(cur);
            PhyloNode* other = findOtherChild(cur, par, prev);
            const nuc_one_hot* other_a = other
                ? cf->getMajorArrayForNode(other)
                : cf->getMajorArrayForNode(cur);
            src_to_root.push_back({cur, other_a, cf->getMajorArrayForNode(cur)});
            prev = cur;
            cur = par;
        }
    }

    collectDiffs(cf, src, base_diffs);
    collectDiffs(cf, sib, base_diffs);
    collectDiffs(cf, sp, base_diffs);
    sort(base_diffs.begin(), base_diffs.end());
    base_diffs.erase(unique(base_diffs.begin(), base_diffs.end()), base_diffs.end());
}

int SPRSourceState::evaluate(PhyloNode* dst, PhyloNode* dst_parent) const {
    if (!valid) return 0;

    PhyloNode* lca = SPRDeltaExact::findLCA(sp, dst, nullptr);
    if (!lca) return 0;

    if (lca == sp) {
        return evalCaseA(dst, dst_parent);
    } else if (dst == lca) {
        return evalCaseB1(dst);
    } else {
        return evalCaseB2(dst, dst_parent, lca);
    }
}

// ========== CASE A: dst in sib's subtree ==========
int SPRSourceState::evalCaseA(PhyloNode* dst, PhyloNode* dst_parent) const {
    // Build dp_to_sib path (dst_parent → ... → sp, exclusive)
    vector<PathNode> dp_to_sib;
    {
        PhyloNode* prev = dst;
        PhyloNode* cur = dst_parent;
        while (cur && cur != sp) {
            PhyloNode* par = SPRMutationOps::getParent(cur);
            PhyloNode* other = findOtherChild(cur, par, prev);
            const nuc_one_hot* other_a = other
                ? cf->getMajorArrayForNode(other)
                : cf->getMajorArrayForNode(cur);
            dp_to_sib.push_back({cur, other_a, cf->getMajorArrayForNode(cur)});
            prev = cur;
            cur = par;
        }
    }

    // Collect affected patterns: base_diffs + dst + dp_to_sib edges
    vector<int> affected(base_diffs);
    collectDiffs(cf, dst, affected);
    for (const auto& pn : dp_to_sib) collectDiffs(cf, pn.node, affected);
    sort(affected.begin(), affected.end());
    affected.erase(unique(affected.begin(), affected.end()), affected.end());

    const nuc_one_hot* dst_arr = cf->getMajorArrayForNode(dst);
    int score = 0;

    for (int p : affected) {
        int freq = cf->getPatternFreq(p);
        nuc_one_hot src_f = src_arr[p], sib_f = sib_arr[p];
        nuc_one_hot dst_f = dst_arr[p], sp_f = sp_arr[p];

        // 1. sp penalty change
        score += (((src_f & dst_f) ? 0 : 1) - ((src_f & sib_f) ? 0 : 1)) * freq;

        // 2. Propagate dp→sib
        nuc_one_hot old_f = dst_f;
        nuc_one_hot new_f = fitchMerge(src_f, dst_f);
        propagate(dp_to_sib, p, freq, old_f, new_f, score);

        nuc_one_hot new_sib_f = (!dp_to_sib.empty() && old_f != new_f) ? new_f : sib_f;

        // 3. At gp
        if (gp->isLeaf()) {
            score += rootEdgeDelta(gp_arr[p], sp_f, new_sib_f) * freq;
        } else {
            nuc_one_hot gp_other_f = gp_other_arr ? gp_other_arr[p] : 0xF;
            nuc_one_hot new_gp_f;
            score += penaltyDelta(sp_f, new_sib_f, gp_other_f, new_gp_f) * freq;

            if (new_gp_f != gp_arr[p]) {
                nuc_one_hot ao = gp_arr[p], an = new_gp_f;
                propagate(above_gp, p, freq, ao, an, score);
            }
        }
    }

    return score;
}

// ========== CASE B1: dst == lca ==========
int SPRSourceState::evalCaseB1(PhyloNode* dst) const {
    PhyloNode* lca = dst;

    // Find src_path_len: how many nodes in src_to_root before lca
    size_t src_path_len = 0;
    for (size_t i = 0; i < src_to_root.size(); i++) {
        if (src_to_root[i].node == lca) { src_path_len = i; break; }
        src_path_len = i + 1;
    }

    PhyloNode* src_branch = (src_path_len == 0)
        ? sp : src_to_root[src_path_len - 1].node;
    const nuc_one_hot* src_branch_arr = cf->getMajorArrayForNode(src_branch);

    PhyloNode* lca_parent = SPRMutationOps::getParent(lca);
    PhyloNode* lca_other = findOtherChild(lca, lca_parent, src_branch);
    const nuc_one_hot* lca_other_arr = lca_other ? cf->getMajorArrayForNode(lca_other) : nullptr;
    const nuc_one_hot* lca_arr = cf->getMajorArrayForNode(lca);

    // Collect affected patterns: base_diffs + src_path nodes + lca
    vector<int> affected(base_diffs);
    for (size_t i = 0; i < src_path_len; i++) collectDiffs(cf, src_to_root[i].node, affected);
    collectDiffs(cf, lca, affected);
    sort(affected.begin(), affected.end());
    affected.erase(unique(affected.begin(), affected.end()), affected.end());

    int score = 0;

    for (int p : affected) {
        int freq = cf->getPatternFreq(p);
        nuc_one_hot src_f = src_arr[p], sib_f = sib_arr[p], sp_f = sp_arr[p];
        int old_sp_pen = (src_f & sib_f) ? 0 : 1;

        // 1. src-side propagation: sp_f → sib_f
        nuc_one_hot src_old = sp_f, src_new = sib_f;
        for (size_t i = 0; i < src_path_len; i++) {
            if (src_old == src_new) break;
            const auto& pn = src_to_root[i];
            nuc_one_hot nf;
            score += penaltyDelta(src_old, src_new, pn.other_arr[p], nf) * freq;
            src_old = pn.node_arr[p];
            src_new = nf;
        }

        nuc_one_hot src_child_old = src_branch_arr[p];
        nuc_one_hot src_child_new = (src_old != src_new) ? src_new : src_child_old;

        // 2. sp_new penalty
        score += (((src_f & src_child_new) ? 0 : 1) - old_sp_pen) * freq;

        // 3. sp_new Fitch
        nuc_one_hot sp_new_f = fitchMerge(src_f, src_child_new);

        // 4. At lca: child changed from src_child_old to sp_new_f
        nuc_one_hot lca_other_f = lca_other_arr ? lca_other_arr[p] : 0xF;
        nuc_one_hot old_lca_f = lca_arr[p];
        nuc_one_hot new_lca_f;
        score += penaltyDelta(src_child_old, sp_new_f, lca_other_f, new_lca_f) * freq;

        // 5. Propagate above lca
        if (!lca_parent && lca->isLeaf()) {
            score += rootEdgeDelta(lca_arr[p], src_child_old, sp_new_f) * freq;
        } else if (new_lca_f != old_lca_f && lca_parent) {
            nuc_one_hot ao = old_lca_f, an = new_lca_f;
            if (lca_parent->isLeaf()) {
                const nuc_one_hot* plca_arr = cf->getMajorArrayForNode(lca_parent);
                score += rootEdgeDelta(plca_arr[p], old_lca_f, new_lca_f) * freq;
            } else {
                propagateRange(src_path_len + 1, src_to_root.size(), p, freq, ao, an, score);
            }
        }
    }

    return score;
}

// ========== CASE B2: general ==========
int SPRSourceState::evalCaseB2(PhyloNode* dst, PhyloNode* dst_parent, PhyloNode* lca) const {
    // Find src_path_len: gp → ... → node before lca
    size_t src_path_len = 0;
    for (size_t i = 0; i < src_to_root.size(); i++) {
        if (src_to_root[i].node == lca) { src_path_len = i; break; }
        src_path_len = i + 1;
    }

    PhyloNode* src_branch = (src_path_len == 0)
        ? sp : src_to_root[src_path_len - 1].node;
    const nuc_one_hot* src_branch_arr = cf->getMajorArrayForNode(src_branch);

    // Build dst path (dst_parent → ... → lca, exclusive)
    vector<PathNode> dst_path;
    {
        PhyloNode* prev = dst;
        PhyloNode* cur = dst_parent;
        while (cur && cur != lca) {
            PhyloNode* par = SPRMutationOps::getParent(cur);
            PhyloNode* other = findOtherChild(cur, par, prev);
            const nuc_one_hot* other_a = other
                ? cf->getMajorArrayForNode(other)
                : cf->getMajorArrayForNode(cur);
            dst_path.push_back({cur, other_a, cf->getMajorArrayForNode(cur)});
            prev = cur;
            cur = par;
        }
    }

    PhyloNode* dst_branch = dst_path.empty() ? dst : dst_path.back().node;
    const nuc_one_hot* dst_branch_arr = cf->getMajorArrayForNode(dst_branch);
    const nuc_one_hot* lca_arr = cf->getMajorArrayForNode(lca);
    const nuc_one_hot* dst_arr = cf->getMajorArrayForNode(dst);

    vector<int> affected(base_diffs);
    collectDiffs(cf, dst, affected);
    for (size_t i = 0; i < src_path_len; i++) collectDiffs(cf, src_to_root[i].node, affected);
    for (const auto& pn : dst_path) collectDiffs(cf, pn.node, affected);
    sort(affected.begin(), affected.end());
    affected.erase(unique(affected.begin(), affected.end()), affected.end());

    int score = 0;

    for (int p : affected) {
        int freq = cf->getPatternFreq(p);
        nuc_one_hot src_f = src_arr[p], sib_f = sib_arr[p];
        nuc_one_hot dst_f = dst_arr[p], sp_f = sp_arr[p];

        // 1. sp penalty change
        score += (((src_f & dst_f) ? 0 : 1) - ((src_f & sib_f) ? 0 : 1)) * freq;

        // 2. src-side: sp_f → sib_f
        nuc_one_hot src_old = sp_f, src_new = sib_f;
        for (size_t i = 0; i < src_path_len; i++) {
            if (src_old == src_new) break;
            const auto& pn = src_to_root[i];
            nuc_one_hot nf;
            score += penaltyDelta(src_old, src_new, pn.other_arr[p], nf) * freq;
            src_old = pn.node_arr[p];
            src_new = nf;
        }

        // 3. dst-side: dst_f → sp_new_f
        nuc_one_hot sp_new_f = fitchMerge(src_f, dst_f);
        nuc_one_hot dst_old = dst_f, dst_new = sp_new_f;
        propagate(dst_path, p, freq, dst_old, dst_new, score);

        // 4. At LCA: combine changes
        bool src_ch = (src_old != src_new), dst_ch = (dst_old != dst_new);
        if (!src_ch && !dst_ch) continue;

        nuc_one_hot old_lca_f = lca_arr[p];
        nuc_one_hot sc_old = src_branch_arr[p];
        nuc_one_hot dc_old = dst_branch_arr[p];
        nuc_one_hot sc_new = src_ch ? src_new : sc_old;
        nuc_one_hot dc_new = dst_ch ? dst_new : dc_old;

        int old_pen = (sc_old & dc_old) ? 0 : 1;
        nuc_one_hot new_lca_f = fitchMerge(sc_new, dc_new);
        int new_pen = (sc_new & dc_new) ? 0 : 1;
        score += (new_pen - old_pen) * freq;

        // 5. Propagate above LCA
        if (new_lca_f != old_lca_f) {
            nuc_one_hot ao = old_lca_f, an = new_lca_f;
            propagateRange(src_path_len + 1, src_to_root.size(), p, freq, ao, an, score);
        }
    }

    return score;
}

struct SPRCandidate {
    PhyloNode* src;
    PhyloNode* src_parent;
    PhyloNode* sib1;
    PhyloNode* sib2;
    PhyloNode* dst;
    PhyloNode* dst_parent;
    int delta;
};

struct NeiSave { Neighbor* nei; Node* orig; };

static vector<PhyloNode*> collectAllNodes(PhyloTree* tree) {
    vector<PhyloNode*> nodes;
    set<PhyloNode*> visited;
    queue<PhyloNode*> q;
    q.push((PhyloNode*)tree->root);
    while (!q.empty()) {
        PhyloNode* n = q.front(); q.pop();
        if (visited.count(n)) continue;
        visited.insert(n);
        nodes.push_back(n);
        FOR_NEIGHBOR_IT(n, nullptr, it)
            if (!visited.count((PhyloNode*)(*it)->node))
                q.push((PhyloNode*)(*it)->node);
    }
    return nodes;
}

static void orientTreeToRoot(PhyloTree* tree) {
    PhyloNode* root = (PhyloNode*)tree->root;

    set<PhyloNode*> visited;
    queue<pair<PhyloNode*, PhyloNode*>> q;
    q.push({root, nullptr});

    while (!q.empty()) {
        auto [node, parent] = q.front(); q.pop();
        if (visited.count(node)) continue;
        visited.insert(node);

        if (parent) {
            for (size_t i = 0; i < node->neighbors.size(); i++) {
                if (node->neighbors[i]->node == parent) {
                    if (i != 0) swap(node->neighbors[0], node->neighbors[i]);
                    break;
                }
            }
        }

        FOR_NEIGHBOR_IT(node, nullptr, it) {
            PhyloNode* child = (PhyloNode*)(*it)->node;
            if (!visited.count(child)) q.push({child, node});
        }
    }

    SPRMutationOps::setRoot(root);
}

static void applySPRMove(PhyloNode* sp, PhyloNode* sib1, PhyloNode* sib2,
                         PhyloNode* dst, PhyloNode* dp) {
    if (dp == sib1) {
        sib1->updateNeighbor(sp, sib2);  sib2->updateNeighbor(sp, sib1);
        sib1->updateNeighbor(dst, sp);   sp->updateNeighbor(sib2, dst);
        dst->updateNeighbor(sib1, sp);
    } else if (dp == sib2) {
        sib2->updateNeighbor(sp, sib1);  sib1->updateNeighbor(sp, sib2);
        sib2->updateNeighbor(dst, sp);   sp->updateNeighbor(sib1, dst);
        dst->updateNeighbor(sib2, sp);
    } else if (dst == sib1) {
        dst->updateNeighbor(sp, sib2);   sib2->updateNeighbor(sp, dst);
        dst->updateNeighbor(dp, sp);     sp->updateNeighbor(sib2, dp);
        dp->updateNeighbor(dst, sp);
    } else if (dst == sib2) {
        dst->updateNeighbor(sp, sib1);   sib1->updateNeighbor(sp, dst);
        dst->updateNeighbor(dp, sp);     sp->updateNeighbor(sib1, dp);
        dp->updateNeighbor(dst, sp);
    } else {
        sib1->updateNeighbor(sp, sib2);  sib2->updateNeighbor(sp, sib1);
        sp->updateNeighbor(sib1, dst);   sp->updateNeighbor(sib2, dp);
        dst->updateNeighbor(dp, sp);     dp->updateNeighbor(dst, sp);
    }
}

static void saveTopology(Node* n, vector<NeiSave>& saves) {
    for (auto it = n->neighbors.begin(); it != n->neighbors.end(); it++)
        saves.push_back({*it, (*it)->node});
}

static void undoTopology(vector<NeiSave>& saves) {
    for (auto& s : saves) s.nei->node = s.orig;
}

static void markDirty(const vector<SPRCandidate>& moves, set<PhyloNode*>& dirty,
                      int max_nodes) {
    dirty.clear();
    for (const auto& m : moves) {
        dirty.insert(m.src);
        dirty.insert(m.src_parent);
        dirty.insert(m.sib1);
        dirty.insert(m.sib2);
        dirty.insert(m.dst);
        dirty.insert(m.dst_parent);
        int walk = 0;
        for (PhyloNode* n = SPRMutationOps::getParent(m.src_parent);
             n && walk < max_nodes; n = SPRMutationOps::getParent(n), walk++)
            dirty.insert(n);
        walk = 0;
        for (PhyloNode* n = SPRMutationOps::getParent(m.dst_parent);
             n && walk < max_nodes; n = SPRMutationOps::getParent(n), walk++)
            dirty.insert(n);
    }
}

static bool isNeighborhoodDirty(PhyloNode* node, int radius, const set<PhyloNode*>& dirty) {
    if (dirty.empty()) return true;
    if (dirty.count(node)) return true;

    queue<pair<PhyloNode*, int>> q;
    set<PhyloNode*> visited;
    q.push({node, 0});
    while (!q.empty()) {
        auto [n, d] = q.front(); q.pop();
        if (visited.count(n)) continue;
        visited.insert(n);
        if (dirty.count(n)) return true;
        if (d < radius) {
            FOR_NEIGHBOR_IT(n, nullptr, it) {
                PhyloNode* nb = (PhyloNode*)(*it)->node;
                if (!visited.count(nb)) q.push({nb, d + 1});
            }
        }
    }
    return false;
}

SPROptimizer::SPROptimizer(PhyloTree* tree) : tree(tree), current_parsimony_score(0),
    best_score_seen(0) {
    if (!tree) throw std::invalid_argument("SPROptimizer: tree cannot be NULL");
}

SPROptimizer::~SPROptimizer() {}

// Select non-conflicting moves greedily (best delta first).
static vector<SPRCandidate> selectMoves(vector<SPRCandidate>& candidates) {
    sort(candidates.begin(), candidates.end(),
         [](const SPRCandidate& a, const SPRCandidate& b) { return a.delta < b.delta; });

    set<Node*> used;
    vector<SPRCandidate> selected;
    for (const auto& m : candidates) {
        if (used.count(m.src) || used.count(m.src_parent) ||
            used.count(m.sib1) || used.count(m.sib2) ||
            used.count(m.dst) || used.count(m.dst_parent))
            continue;
        selected.push_back(m);
        used.insert(m.src);    used.insert(m.src_parent);
        used.insert(m.sib1);   used.insert(m.sib2);
        used.insert(m.dst);    used.insert(m.dst_parent);
    }
    return selected;
}

int SPROptimizer::batchSPROptimize(int radius, Fitch& cf) {
    int cur_score = cf.recompute();
    int initial_score = cur_score;

    orientTreeToRoot(tree);
    SPRDeltaExact::precomputeDepths(tree);

    string radius_str = (radius == 0) ? "unbounded" : to_string(radius);
    cout << "=== Batch SPR (radius " << radius_str << ") ===" << endl;
    cout << "  Starting score: " << cur_score << endl;

    set<PhyloNode*> dirty_nodes;
    int dirty_check_radius = (radius == 0) ? 32 : radius;

    for (int round = 0; round < MAX_ROUNDS_PER_RADIUS; round++) {
        auto t0 = high_resolution_clock::now();
        vector<PhyloNode*> all_nodes = collectAllNodes(tree);
        vector<SPRCandidate> candidates;
        int moves_evaluated = 0, src_skipped = 0;

        for (PhyloNode* src : all_nodes) {
            if (src == tree->root || src->degree() < 2) continue;

            PhyloNode* sp = SPRMutationOps::getParent(src);
            if (!sp || sp->degree() != BINARY_NODE_DEGREE) continue;

            if (!dirty_nodes.empty() &&
                !isNeighborhoodDirty(src, dirty_check_radius, dirty_nodes)) {
                src_skipped++;
                continue;
            }

            PhyloNode* sib1 = NULL, *sib2 = NULL;
            FOR_NEIGHBOR_IT(sp, src, it) {
                PhyloNode* nb = (PhyloNode*)(*it)->node;
                if (!sib1) sib1 = nb; else if (!sib2) sib2 = nb;
            }
            if (!sib1 || !sib2) continue;

            SPRSourceState state(&cf, src, sp);
            if (!state.isValid()) continue;

            set<PhyloNode*> visited;
            visited.insert(src);
            visited.insert(sp);

            function<void(PhyloNode*, PhyloNode*, int)> dfs =
                [&](PhyloNode* node, PhyloNode* from, int dist) {
                if (visited.count(node)) return;
                visited.insert(node);

                if (from != sp && from != src && node != src && node != sp) {
                    int delta = state.evaluate(node, from);
                    moves_evaluated++;
                    if (delta < 0)
                        candidates.push_back({src, sp, sib1, sib2, node, from, delta});
                }

                if (radius > 0 && dist >= radius) return;
                FOR_NEIGHBOR_IT(node, nullptr, nit) {
                    PhyloNode* next = (PhyloNode*)(*nit)->node;
                    if (!visited.count(next)) dfs(next, node, dist + 1);
                }
            };

            dfs(sib1, sp, 0);
            dfs(sib2, sp, 0);
        }

        long long elapsed = duration_cast<milliseconds>(high_resolution_clock::now() - t0).count();

        if (candidates.empty()) {
            cout << "  Round " << round + 1 << ": no improving moves ("
                 << moves_evaluated << " eval";
            if (src_skipped > 0) cout << ", " << src_skipped << " skip";
            cout << ", " << elapsed << "ms)" << endl << flush;
            break;
        }

        vector<SPRCandidate> selected = selectMoves(candidates);

        vector<NeiSave> saves;
        for (const auto& m : selected) {
            saveTopology(m.src_parent, saves);
            saveTopology(m.sib1, saves);
            saveTopology(m.sib2, saves);
            saveTopology(m.dst, saves);
            saveTopology(m.dst_parent, saves);
        }

        int expected = 0;
        for (const auto& m : selected)
            { applySPRMove(m.src_parent, m.sib1, m.sib2, m.dst, m.dst_parent); expected += m.delta; }

        int new_score = cf.recompute();
        orientTreeToRoot(tree);
        SPRDeltaExact::precomputeDepths(tree);

        cout << "  Round " << round + 1 << ": " << candidates.size() << " found, "
             << selected.size() << " applied (exp=" << expected
             << " act=" << (new_score - cur_score) << ") -> " << new_score
             << " (" << moves_evaluated << " eval";
        if (src_skipped > 0) cout << ", " << src_skipped << " skip";
        cout << ", " << elapsed << "ms)" << endl << flush;

        if (new_score >= cur_score) {
            undoTopology(saves);
            cf.recompute();
            orientTreeToRoot(tree);
            SPRDeltaExact::precomputeDepths(tree);
            cout << "  No improvement - reverted" << endl;
            break;
        }

        markDirty(selected, dirty_nodes, (int)all_nodes.size());
        cur_score = new_score;
        SPRDeltaExact::setCustomFitch(&cf, cur_score);
    }

    cout << "  Result: " << initial_score << " -> " << cur_score
         << " (delta=" << (initial_score - cur_score) << ")" << endl;
    return cur_score;
}

int SPROptimizer::optimizeTree(int max_passes) {
    Fitch cf(tree);
    int initial_score = cf.run();
    SPRMutationOps::setCustomFitch(&cf);
    SPRDeltaExact::setCustomFitch(&cf, initial_score);

    cout << "Fitch: score=" << initial_score
         << " (max_passes=" << max_passes << ")" << endl;

    for (int pass = 0; pass < max_passes; pass++) {
        int start = cf.recompute();
        SPRDeltaExact::setCustomFitch(&cf, start);

        for (int i = 0; i < NUM_SPR_RADII; i++)
            batchSPROptimize(SPR_RADII[i], cf);

        int end = cf.recompute();
        cout << "Pass " << pass + 1 << ": " << start << " -> " << end
             << " (delta=" << (start - end) << ")" << endl;
        if (end >= start) break;
    }

    current_parsimony_score = cf.recompute();
    best_score_seen = current_parsimony_score;

    cout << "\nSPR complete: " << initial_score << " -> " << current_parsimony_score
         << " (total delta=" << (initial_score - current_parsimony_score) << ")" << endl;
    return current_parsimony_score;
}
