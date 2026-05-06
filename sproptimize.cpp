#include "sproptimize.h"
#include "spr_context.h"
#include "spr_delta_exact.h"
#include "spr_mutation_ops.h"
#include "spr_utils.h"
#include "phylotree.h"
#include "phylonode.h"
#include <iostream>
#include <queue>
#include <set>
#include <chrono>
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <random>
#include <cmath>

using namespace std;
using namespace std::chrono;

static const int MAX_ROUNDS_PER_RADIUS = 100;
static const int RADIUS_STALL_LIMIT = 2;
static const int BINARY_NODE_DEGREE = 3;
static const double CONVERGENCE_THRESHOLD = 0.001;
static const double RATCHET_INV_TEMPERATURE = 2.0;
static const int    RATCHET_ZERO_ACCEPT_PCT = 50;

namespace {
struct OptimizerScratch {
    PhyloNode*           root = nullptr;
    vector<bool>         affected_bitset;
    vector<int>          affected_list;
    int                  bitset_nptn = 0;
    vector<PathStep>     scratch_path;
    vector<bool>         orient_visited;
    vector<int>          orient_visited_ids;
    vector<bool>         bfs_visited;
    vector<int>          bfs_visited_ids;
    vector<bool>         select_used;
    vector<int>          select_used_ids;
};
} // namespace

static OptimizerScratch g_scratch;

/**
 * O(1) parent lookup. Requires orientTreeToRoot() since last topology change.
 * @param node Node to query.
 * @return Parent node, or nullptr for the root or a null input.
 */
static inline PhyloNode* getParent(PhyloNode* node) {
    if (!node || node == g_scratch.root) return nullptr;
    return (PhyloNode*)node->neighbors[0]->node;
}

/**
 * Resizes affected bitset if pattern count changed.
 * @param nptn Required pattern count.
 */
static inline void bitsetInit(int nptn) {
    if (g_scratch.bitset_nptn != nptn) {
        g_scratch.affected_bitset.assign(nptn, false);
        g_scratch.affected_list.reserve(nptn);
        g_scratch.bitset_nptn = nptn;
    }
}

/** Clears the affected bitset in O(|set|). */
static inline void bitsetClear() {
    for (int p : g_scratch.affected_list) g_scratch.affected_bitset[p] = false;
    g_scratch.affected_list.clear();
}

/**
 * Adds a node's diff list to the affected bitset.
 * @param diffs Sorted-unique diff list, or nullptr.
 */
static inline void bitsetAddDiffs(const vector<int>* diffs) {
    if (!diffs) return;
    for (int p : *diffs) {
        if (!g_scratch.affected_bitset[p]) {
            g_scratch.affected_bitset[p] = true;
            g_scratch.affected_list.push_back(p);
        }
    }
}

/**
 * Per-source precomputed state for O(M) delta evaluation.
 */
class SPRSourceState {
public:
    SPRSourceState(PhyloTree* tree, PhyloNode* src, PhyloNode* src_parent);

    int evaluate(PhyloNode* dst, PhyloNode* dst_parent) const;
    bool isValid() const { return valid; }

    bool hasImprovementPotential() const { return has_nonzero_diffs; }

    void setupBaseBitset() const;
    void restoreBaseBitset() const;

private:
    PhyloTree* tree;
    PhyloNode* src;
    PhyloNode* src_parent;
    PhyloNode* grandparent;
    PhyloNode* sibling;
    bool valid;
    bool has_nonzero_diffs;

    const nuc_one_hot* src_states;
    const nuc_one_hot* sibling_states;
    const nuc_one_hot* src_parent_states;
    const nuc_one_hot* grandparent_states;

    PhyloNode* gp_sibling;
    const nuc_one_hot* gp_sibling_states;
    vector<PathStep> src_to_root;

    mutable int base_bitset_size;  // number of base_diffs entries in g_scratch.affected_list

    int evalCaseA(PhyloNode* dst, PhyloNode* dst_parent) const;
    int evalCaseB1(PhyloNode* dst) const;
    int evalCaseB2(PhyloNode* dst, PhyloNode* dst_parent, PhyloNode* lca) const;
};

SPRSourceState::SPRSourceState(PhyloTree* tree, PhyloNode* src, PhyloNode* src_parent)
    : tree(tree), src(src), src_parent(src_parent), valid(false),
      has_nonzero_diffs(false) {
    grandparent = getParent(src_parent);
    sibling = findOtherChild(src_parent, grandparent, src);
    if (!sibling || !grandparent) return;
    valid = true;

    src_states = tree->fitchMajorArrayFor(src);
    sibling_states = tree->fitchMajorArrayFor(sibling);
    src_parent_states = tree->fitchMajorArrayFor(src_parent);
    grandparent_states = tree->fitchMajorArrayFor(grandparent);

    gp_sibling = findOtherChild(grandparent, getParent(grandparent), src_parent);
    gp_sibling_states = gp_sibling ? tree->fitchMajorArrayFor(gp_sibling) : nullptr;

    {
        PhyloNode* prev = src_parent;
        PhyloNode* cur = grandparent;
        while (cur) {
            PhyloNode* cur_parent = getParent(cur);
            PhyloNode* other = findOtherChild(cur, cur_parent, prev);
            const nuc_one_hot* other_states = other
                ? tree->fitchMajorArrayFor(other)
                : tree->fitchMajorArrayFor(cur);
            src_to_root.push_back({cur, other_states, tree->fitchMajorArrayFor(cur)});
            prev = cur;
            cur = cur_parent;
        }
    }

    {
        const vector<int>* d;
        d = tree->fitchDiffsFor(src);       if (d && !d->empty()) has_nonzero_diffs = true;
        d = tree->fitchDiffsFor(sibling);   if (d && !d->empty()) has_nonzero_diffs = true;
        d = tree->fitchDiffsFor(src_parent); if (d && !d->empty()) has_nonzero_diffs = true;
    }
    base_bitset_size = 0;
}

void SPRSourceState::setupBaseBitset() const {
    bitsetClear();
    bitsetAddDiffs(tree->fitchDiffsFor(src));
    bitsetAddDiffs(tree->fitchDiffsFor(sibling));
    bitsetAddDiffs(tree->fitchDiffsFor(src_parent));
    base_bitset_size = (int)g_scratch.affected_list.size();
}

void SPRSourceState::restoreBaseBitset() const {
    // Remove dst-specific entries (those added after base setup)
    for (int i = base_bitset_size; i < (int)g_scratch.affected_list.size(); i++) {
        g_scratch.affected_bitset[g_scratch.affected_list[i]] = false;
    }
    g_scratch.affected_list.resize(base_bitset_size);
}

int SPRSourceState::evaluate(PhyloNode* dst, PhyloNode* dst_parent) const {
    if (!valid) return 0;

    PhyloNode* lca = SPRDeltaExact::findLCA(src_parent, dst, nullptr);
    if (!lca) return 0;

    if (lca == src_parent) {
        return evalCaseA(dst, dst_parent);
    } else if (dst == lca) {
        return evalCaseB1(dst);
    } else {
        return evalCaseB2(dst, dst_parent, lca);
    }
}

// ========== CASE A: dst in sibling's subtree ==========
int SPRSourceState::evalCaseA(PhyloNode* dst, PhyloNode* dst_parent) const {
    g_scratch.scratch_path.clear();
    auto& scratch_path = g_scratch.scratch_path;
    {
        PhyloNode* prev = dst;
        PhyloNode* cur = dst_parent;
        while (cur && cur != src_parent) {
            PhyloNode* cur_parent = getParent(cur);
            PhyloNode* other = findOtherChild(cur, cur_parent, prev);
            const nuc_one_hot* other_states = other
                ? tree->fitchMajorArrayFor(other)
                : tree->fitchMajorArrayFor(cur);
            scratch_path.push_back({cur, other_states, tree->fitchMajorArrayFor(cur)});
            prev = cur;
            cur = cur_parent;
        }
    }

    // Add dst-specific diffs on top of persistent base bitset
    bitsetAddDiffs(tree->fitchDiffsFor(dst));
    for (const auto& step : scratch_path) bitsetAddDiffs(tree->fitchDiffsFor(step.node));

    const nuc_one_hot* dst_states = tree->fitchMajorArrayFor(dst);
    int score = 0;

    for (int ptn : g_scratch.affected_list) {
        int freq = tree->fitchPatternFreq(ptn);
        nuc_one_hot src_fitch = src_states[ptn], sibling_fitch = sibling_states[ptn];
        nuc_one_hot dst_fitch = dst_states[ptn], sp_fitch = src_parent_states[ptn];

        score += (((src_fitch & dst_fitch) ? 0 : 1) - ((src_fitch & sibling_fitch) ? 0 : 1)) * freq;

        nuc_one_hot old_fitch = dst_fitch;
        nuc_one_hot new_fitch = fitchMerge(src_fitch, dst_fitch);
        propagatePath(scratch_path, ptn, freq, old_fitch, new_fitch, score);

        nuc_one_hot new_sibling_fitch = (!scratch_path.empty() && old_fitch != new_fitch)
            ? new_fitch : sibling_fitch;

        if (grandparent->isLeaf()) {
            score += rootEdgeDelta(grandparent_states[ptn], sp_fitch, new_sibling_fitch) * freq;
        } else {
            nuc_one_hot gp_sibling_fitch = gp_sibling_states ? gp_sibling_states[ptn] : 0xF;
            nuc_one_hot new_gp_fitch;
            score += penaltyDelta(sp_fitch, new_sibling_fitch, gp_sibling_fitch, new_gp_fitch) * freq;

            if (new_gp_fitch != grandparent_states[ptn]) {
                nuc_one_hot old_prop = grandparent_states[ptn], new_prop = new_gp_fitch;
                propagatePath(src_to_root, 1, src_to_root.size(), ptn, freq, old_prop, new_prop, score);
            }
        }
    }

    restoreBaseBitset();
    return score;
}

// ========== CASE B1: dst == lca ==========
int SPRSourceState::evalCaseB1(PhyloNode* dst) const {
    PhyloNode* lca = dst;

    size_t src_path_len = 0;
    for (size_t i = 0; i < src_to_root.size(); i++) {
        if (src_to_root[i].node == lca) { src_path_len = i; break; }
        src_path_len = i + 1;
    }

    PhyloNode* src_branch = (src_path_len == 0)
        ? src_parent : src_to_root[src_path_len - 1].node;
    const nuc_one_hot* src_branch_states = tree->fitchMajorArrayFor(src_branch);

    PhyloNode* lca_parent = getParent(lca);
    PhyloNode* lca_sibling = findOtherChild(lca, lca_parent, src_branch);
    const nuc_one_hot* lca_sibling_states = lca_sibling ? tree->fitchMajorArrayFor(lca_sibling) : nullptr;
    const nuc_one_hot* lca_states = tree->fitchMajorArrayFor(lca);

    // Add dst-specific diffs on top of persistent base bitset
    for (size_t i = 0; i < src_path_len; i++) bitsetAddDiffs(tree->fitchDiffsFor(src_to_root[i].node));
    bitsetAddDiffs(tree->fitchDiffsFor(lca));

    int score = 0;

    for (int ptn : g_scratch.affected_list) {
        int freq = tree->fitchPatternFreq(ptn);
        nuc_one_hot src_fitch = src_states[ptn], sibling_fitch = sibling_states[ptn];
        nuc_one_hot sp_fitch = src_parent_states[ptn];
        int old_sp_penalty = (src_fitch & sibling_fitch) ? 0 : 1;

        nuc_one_hot src_old = sp_fitch, src_new = sibling_fitch;
        propagatePath(src_to_root, 0, src_path_len, ptn, freq, src_old, src_new, score);

        nuc_one_hot src_child_old = src_branch_states[ptn];
        nuc_one_hot src_child_new = (src_old != src_new) ? src_new : src_child_old;

        score += (((src_fitch & src_child_new) ? 0 : 1) - old_sp_penalty) * freq;

        nuc_one_hot new_sp_fitch = fitchMerge(src_fitch, src_child_new);

        nuc_one_hot lca_sibling_fitch = lca_sibling_states ? lca_sibling_states[ptn] : 0xF;
        nuc_one_hot old_lca_fitch = lca_states[ptn];
        nuc_one_hot new_lca_fitch;
        score += penaltyDelta(src_child_old, new_sp_fitch, lca_sibling_fitch, new_lca_fitch) * freq;

        if (!lca_parent && lca->isLeaf()) {
            score += rootEdgeDelta(lca_states[ptn], src_child_old, new_sp_fitch) * freq;
        } else if (new_lca_fitch != old_lca_fitch && lca_parent) {
            nuc_one_hot old_prop = old_lca_fitch, new_prop = new_lca_fitch;
            if (lca_parent->isLeaf()) {
                const nuc_one_hot* lca_parent_states = tree->fitchMajorArrayFor(lca_parent);
                score += rootEdgeDelta(lca_parent_states[ptn], old_lca_fitch, new_lca_fitch) * freq;
            } else {
                propagatePath(src_to_root, src_path_len + 1, src_to_root.size(), ptn, freq, old_prop, new_prop, score);
            }
        }
    }

    restoreBaseBitset();
    return score;
}

// ========== CASE B2: general ==========
int SPRSourceState::evalCaseB2(PhyloNode* dst, PhyloNode* dst_parent, PhyloNode* lca) const {
    size_t src_path_len = 0;
    for (size_t i = 0; i < src_to_root.size(); i++) {
        if (src_to_root[i].node == lca) { src_path_len = i; break; }
        src_path_len = i + 1;
    }

    PhyloNode* src_branch = (src_path_len == 0)
        ? src_parent : src_to_root[src_path_len - 1].node;
    const nuc_one_hot* src_branch_states = tree->fitchMajorArrayFor(src_branch);

    g_scratch.scratch_path.clear();
    auto& scratch_path = g_scratch.scratch_path;
    {
        PhyloNode* prev = dst;
        PhyloNode* cur = dst_parent;
        while (cur && cur != lca) {
            PhyloNode* cur_parent = getParent(cur);
            PhyloNode* other = findOtherChild(cur, cur_parent, prev);
            const nuc_one_hot* other_states = other
                ? tree->fitchMajorArrayFor(other)
                : tree->fitchMajorArrayFor(cur);
            scratch_path.push_back({cur, other_states, tree->fitchMajorArrayFor(cur)});
            prev = cur;
            cur = cur_parent;
        }
    }

    PhyloNode* dst_branch = scratch_path.empty() ? dst : scratch_path.back().node;
    const nuc_one_hot* dst_branch_states = tree->fitchMajorArrayFor(dst_branch);
    const nuc_one_hot* lca_states = tree->fitchMajorArrayFor(lca);
    const nuc_one_hot* dst_states = tree->fitchMajorArrayFor(dst);

    // Add dst-specific diffs on top of persistent base bitset
    bitsetAddDiffs(tree->fitchDiffsFor(dst));
    for (size_t i = 0; i < src_path_len; i++) bitsetAddDiffs(tree->fitchDiffsFor(src_to_root[i].node));
    for (const auto& step : scratch_path) bitsetAddDiffs(tree->fitchDiffsFor(step.node));

    int score = 0;

    for (int ptn : g_scratch.affected_list) {
        int freq = tree->fitchPatternFreq(ptn);
        nuc_one_hot src_fitch = src_states[ptn], sibling_fitch = sibling_states[ptn];
        nuc_one_hot dst_fitch = dst_states[ptn], sp_fitch = src_parent_states[ptn];

        score += (((src_fitch & dst_fitch) ? 0 : 1) - ((src_fitch & sibling_fitch) ? 0 : 1)) * freq;

        nuc_one_hot src_old = sp_fitch, src_new = sibling_fitch;
        propagatePath(src_to_root, 0, src_path_len, ptn, freq, src_old, src_new, score);

        nuc_one_hot new_sp_fitch = fitchMerge(src_fitch, dst_fitch);
        nuc_one_hot dst_old = dst_fitch, dst_new = new_sp_fitch;
        propagatePath(scratch_path, ptn, freq, dst_old, dst_new, score);

        bool src_changed = (src_old != src_new), dst_changed = (dst_old != dst_new);
        if (!src_changed && !dst_changed) continue;

        nuc_one_hot old_lca_fitch = lca_states[ptn];
        nuc_one_hot src_child_old_fitch = src_branch_states[ptn];
        nuc_one_hot dst_child_old_fitch = dst_branch_states[ptn];
        nuc_one_hot src_child_new_fitch = src_changed ? src_new : src_child_old_fitch;
        nuc_one_hot dst_child_new_fitch = dst_changed ? dst_new : dst_child_old_fitch;

        int old_penalty = (src_child_old_fitch & dst_child_old_fitch) ? 0 : 1;
        nuc_one_hot new_lca_fitch = fitchMerge(src_child_new_fitch, dst_child_new_fitch);
        int new_penalty = (src_child_new_fitch & dst_child_new_fitch) ? 0 : 1;
        score += (new_penalty - old_penalty) * freq;

        if (new_lca_fitch != old_lca_fitch) {
            nuc_one_hot old_prop = old_lca_fitch, new_prop = new_lca_fitch;
            propagatePath(src_to_root, src_path_len + 1, src_to_root.size(), ptn, freq, old_prop, new_prop, score);
        }
    }

    restoreBaseBitset();
    return score;
}

/**
 * One candidate SPR move.
 */
struct SPRCandidate {
    PhyloNode* src;
    PhyloNode* src_parent;
    PhyloNode* sibling1;
    PhyloNode* sibling2;
    PhyloNode* dst;
    PhyloNode* dst_parent;
    int delta;          ///< <0 = improvement, ==0 = plateau move (drift).
};

using NeighborSave = SPRNeighborSave;

/**
 * BFS from root, returning every reachable node once.
 * @param tree   Tree to traverse.
 * @param max_id Maximum node id.
 * @return Vector of all nodes reachable from tree->root.
 */
static vector<PhyloNode*> collectAllNodes(PhyloTree* tree, int max_id) {
    vector<PhyloNode*> nodes;
    vector<bool> visited(max_id + 1, false);
    queue<PhyloNode*> q;
    q.push((PhyloNode*)tree->root);
    while (!q.empty()) {
        PhyloNode* node = q.front(); q.pop();
        if (visited[node->id]) continue;
        visited[node->id] = true;
        nodes.push_back(node);
        FOR_NEIGHBOR_IT(node, nullptr, nit) {
            PhyloNode* child = (PhyloNode*)(*nit)->node;
            if (!visited[child->id]) q.push(child);
        }
    }
    return nodes;
}

/**
 * Re-orients each node so neighbors[0] points at the parent (rootward).
 * @param tree   Tree to re-orient.
 * @param max_id Maximum node id.
 */
static void orientTreeToRoot(PhyloTree* tree, int max_id) {
    PhyloNode* root = (PhyloNode*)tree->root;

    if ((int)g_scratch.orient_visited.size() < max_id + 1) {
        g_scratch.orient_visited.assign(max_id + 1, false);
        g_scratch.orient_visited_ids.reserve(max_id + 1);
    }
    for (int id : g_scratch.orient_visited_ids) g_scratch.orient_visited[id] = false;
    g_scratch.orient_visited_ids.clear();

    queue<pair<PhyloNode*, PhyloNode*>> q;
    q.push({root, nullptr});

    while (!q.empty()) {
        PhyloNode* node = (PhyloNode*)q.front().first; PhyloNode* parent = (PhyloNode*)q.front().second; q.pop();
        if (g_scratch.orient_visited[node->id]) continue;
        g_scratch.orient_visited[node->id] = true;
        g_scratch.orient_visited_ids.push_back(node->id);

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
            if (!g_scratch.orient_visited[child->id]) q.push({child, node});
        }
    }

    g_scratch.root = root;
}

/**
 * Applies an SPR move by rewiring neighbor pointers.
 * @param src_parent Parent of source (will be relocated).
 * @param sibling1   First non-src neighbor of src_parent.
 * @param sibling2   Second non-src neighbor of src_parent.
 * @param dst        Regraft target.
 * @param dst_parent Parent of dst.
 */
static void applySPRMove(PhyloNode* src_parent, PhyloNode* sibling1, PhyloNode* sibling2,
                         PhyloNode* dst, PhyloNode* dst_parent) {
    if (dst_parent == sibling1) {
        sibling1->updateNeighbor(src_parent, sibling2);  sibling2->updateNeighbor(src_parent, sibling1);
        src_parent->updateNeighbor(sibling1, dst);        src_parent->updateNeighbor(sibling2, dst_parent);
        dst->updateNeighbor(dst_parent, src_parent);      dst_parent->updateNeighbor(dst, src_parent);
    } else if (dst_parent == sibling2) {
        sibling2->updateNeighbor(src_parent, sibling1);   sibling1->updateNeighbor(src_parent, sibling2);
        src_parent->updateNeighbor(sibling2, dst);         src_parent->updateNeighbor(sibling1, dst_parent);
        dst->updateNeighbor(dst_parent, src_parent);       dst_parent->updateNeighbor(dst, src_parent);
    } else {
        sibling1->updateNeighbor(src_parent, sibling2);  sibling2->updateNeighbor(src_parent, sibling1);
        src_parent->updateNeighbor(sibling1, dst);       src_parent->updateNeighbor(sibling2, dst_parent);
        dst->updateNeighbor(dst_parent, src_parent);     dst_parent->updateNeighbor(dst, src_parent);
    }
}

/**
 * Snapshots every neighbor pointer on a node for later undo.
 * @param node  Node whose neighbors to snapshot.
 * @param saves Output collection.
 */
static inline void saveTopology(PhyloNode* node, vector<NeighborSave>& saves) {
    sprSaveNodeTopology(node, saves);
}

/**
 * Restores neighbor pointers from a save list.
 * @param saves Saves recorded by saveTopology.
 */
static inline void undoTopology(vector<NeighborSave>& saves) {
    sprUndoTopology(saves);
}

/**
 * Marks the 6 nodes touched by each move as dirty.
 * @param moves  Selected SPR moves.
 * @param dirty  In/out: dirty flag bitset.
 * @param max_id Upper bound on dirty index.
 */
static void markDirty(const vector<SPRCandidate>& moves, vector<bool>& dirty, int max_id) {
    for (const auto& m : moves) {
        if (m.src->id < max_id) dirty[m.src->id] = true;
        if (m.src_parent->id < max_id) dirty[m.src_parent->id] = true;
        if (m.sibling1->id < max_id) dirty[m.sibling1->id] = true;
        if (m.sibling2->id < max_id) dirty[m.sibling2->id] = true;
        if (m.dst->id < max_id) dirty[m.dst->id] = true;
        if (m.dst_parent->id < max_id) dirty[m.dst_parent->id] = true;
    }
}

/**
 * Tests whether any node within radius edges of `node` is dirty.
 * @param node    Source node to check around.
 * @param radius  Max BFS depth.
 * @param dirty   Dirty bitset (indexed by node id).
 * @param max_id  Upper bound for visited bitset size.
 * @return true if some node within radius is dirty.
 */
static bool isNeighborhoodDirty(PhyloNode* node, int radius, const vector<bool>& dirty, int max_id) {
    if ((int)g_scratch.bfs_visited.size() < max_id + 1) {
        g_scratch.bfs_visited.assign(max_id + 1, false);
        g_scratch.bfs_visited_ids.reserve(max_id + 1);
    }
    for (int id : g_scratch.bfs_visited_ids) g_scratch.bfs_visited[id] = false;
    g_scratch.bfs_visited_ids.clear();

    queue<pair<PhyloNode*, int>> q;
    q.push({node, 0});
    bool found = false;

    while (!q.empty() && !found) {
        PhyloNode* cur = q.front().first;
        int depth = q.front().second;
        q.pop();

        if (g_scratch.bfs_visited[cur->id]) continue;
        g_scratch.bfs_visited[cur->id] = true;
        g_scratch.bfs_visited_ids.push_back(cur->id);

        if (dirty[cur->id]) { found = true; break; }

        if (depth < radius) {
            FOR_NEIGHBOR_IT(cur, nullptr, nit) {
                PhyloNode* neighbor = (PhyloNode*)(*nit)->node;
                if (!g_scratch.bfs_visited[neighbor->id]) {
                    q.push({neighbor, depth + 1});
                }
            }
        }
    }

    for (int id : g_scratch.bfs_visited_ids) g_scratch.bfs_visited[id] = false;
    return found;
}

/**
 * Context passed to searchDestinations.
 */
struct DFSContext {
    const SPRSourceState* state;
    PhyloTree* tree;
    PhyloNode* src;
    PhyloNode* src_parent;
    PhyloNode* sibling1;
    PhyloNode* sibling2;
    vector<SPRCandidate>* candidates;
    vector<bool>* dfs_visited;
    vector<int>* dfs_visited_ids;
    int radius;
    int* moves_evaluated;
    int* moves_pruned;
    int best_delta;
    const nuc_one_hot* src_states;
    bool wall_active;
    std::chrono::high_resolution_clock::time_point wall_deadline;
    bool* wall_dfs_hit;
};

/**
 * Branch-and-bound estimate of the best delta achievable in a subtree.
 * @param ctx          DFS context.
 * @param node         Subtree root being bounded.
 * @param parent_delta Best delta above this subtree.
 * @return Upper bound on best delta within this subtree.
 */
static int computeSubtreeBound(const DFSContext& ctx, PhyloNode* node, int parent_delta) {
    const vector<int>* diffs = ctx.tree->fitchDiffsFor(node);
    if (!diffs || diffs->empty()) return parent_delta;

    int max_improvement = 0;
    const nuc_one_hot* node_states = ctx.tree->fitchMajorArrayFor(node);
    for (int ptn : *diffs) {
        if (ctx.src_states[ptn] & node_states[ptn])
            max_improvement += ctx.tree->fitchPatternFreq(ptn);
    }
    return parent_delta - max_improvement;
}

/**
 * Recursive DFS to discover candidate destinations for a given source.
 * @param ctx          DFS context.
 * @param node         Current node (candidate regraft target).
 * @param from         Neighbor we came from.
 * @param dist         Edge count from entry point.
 * @param parent_delta Best delta above this subtree.
 */
static void searchDestinations(const DFSContext& ctx, PhyloNode* node, PhyloNode* from, int dist, int parent_delta) {
    if (ctx.wall_active) {
        if (*ctx.wall_dfs_hit) return;
        if (std::chrono::high_resolution_clock::now() > ctx.wall_deadline) {
            *ctx.wall_dfs_hit = true;
            return;
        }
    }
    if ((*ctx.dfs_visited)[node->id]) return;
    (*ctx.dfs_visited)[node->id] = true;
    ctx.dfs_visited_ids->push_back(node->id);

    int current_delta = 0;
    bool evaluated = false;

    if (from != ctx.src_parent && from != ctx.src && node != ctx.src && node != ctx.src_parent) {
        current_delta = ctx.state->evaluate(node, from);
        (*ctx.moves_evaluated)++;
        evaluated = true;
        if (current_delta < 0)
            ctx.candidates->push_back({ctx.src, ctx.src_parent, ctx.sibling1, ctx.sibling2, node, from, current_delta});
    }

    if (ctx.radius > 0 && dist >= ctx.radius) return;

    int delta_for_bound = evaluated ? current_delta : parent_delta;

    FOR_NEIGHBOR_IT(node, nullptr, nit) {
        PhyloNode* next = (PhyloNode*)(*nit)->node;
        if ((*ctx.dfs_visited)[next->id]) continue;

        int child_bound = computeSubtreeBound(ctx, next, delta_for_bound);
        if (child_bound >= 0) {
            (*ctx.moves_pruned)++;
            continue;
        }

        searchDestinations(ctx, next, node, dist + 1, delta_for_bound);
    }
}

SPROptimizer::SPROptimizer(PhyloTree* tree) : tree(tree), current_parsimony_score(0),
    best_score_seen(0) {
    if (!tree) throw std::invalid_argument("SPROptimizer: tree cannot be NULL");
}

SPROptimizer::~SPROptimizer() {}


/**
 * Greedily picks a non-conflicting batch of moves (sorted by delta).
 * @param candidates In/out: candidates list (sorted by delta in place).
 * @param max_id     Upper bound for the used-node bitset.
 * @param deferred   Optional output: conflicting candidates skipped.
 * @return Selected non-conflicting moves.
 */
static vector<SPRCandidate> selectMoves(vector<SPRCandidate>& candidates, int max_id,
                                         vector<SPRCandidate>* deferred = nullptr) {
    sort(candidates.begin(), candidates.end(),
         [](const SPRCandidate& a, const SPRCandidate& b) { return a.delta < b.delta; });

    if ((int)g_scratch.select_used.size() < max_id + 1) {
        g_scratch.select_used.assign(max_id + 1, false);
        g_scratch.select_used_ids.reserve(max_id + 1);
    }
    for (int id : g_scratch.select_used_ids) g_scratch.select_used[id] = false;
    g_scratch.select_used_ids.clear();

    vector<SPRCandidate> selected;
    for (const auto& move : candidates) {
        if (g_scratch.select_used[move.src->id] || g_scratch.select_used[move.src_parent->id] ||
            g_scratch.select_used[move.sibling1->id] || g_scratch.select_used[move.sibling2->id] ||
            g_scratch.select_used[move.dst->id] || g_scratch.select_used[move.dst_parent->id]) {
            if (deferred) deferred->push_back(move);
            continue;
        }
        selected.push_back(move);
        g_scratch.select_used[move.src->id] = true;         g_scratch.select_used_ids.push_back(move.src->id);
        g_scratch.select_used[move.src_parent->id] = true;  g_scratch.select_used_ids.push_back(move.src_parent->id);
        g_scratch.select_used[move.sibling1->id] = true;    g_scratch.select_used_ids.push_back(move.sibling1->id);
        g_scratch.select_used[move.sibling2->id] = true;    g_scratch.select_used_ids.push_back(move.sibling2->id);
        g_scratch.select_used[move.dst->id] = true;          g_scratch.select_used_ids.push_back(move.dst->id);
        g_scratch.select_used[move.dst_parent->id] = true;   g_scratch.select_used_ids.push_back(move.dst_parent->id);
    }
    return selected;
}

/**
 * Collects the 6 directly-modified nodes per move.
 * @param moves Selected SPR moves.
 * @return Set of nodes directly altered.
 */
static set<PhyloNode*> collectDirtyNodes(const vector<SPRCandidate>& moves) {
    set<PhyloNode*> dirty;
    for (const auto& m : moves) {
        dirty.insert(m.src);
        dirty.insert(m.src_parent);
        dirty.insert(m.sibling1);
        dirty.insert(m.sibling2);
        dirty.insert(m.dst);
        dirty.insert(m.dst_parent);
    }
    return dirty;
}

/**
 * Collects dirty nodes and all ancestors up to root (for incremental diff update).
 * @param moves Selected SPR moves.
 * @return Set of dirty nodes plus their ancestors.
 */
static set<PhyloNode*> collectDirtyNodesWithAncestors(const vector<SPRCandidate>& moves) {
    set<PhyloNode*> dirty = collectDirtyNodes(moves);
    // Add ancestors of each dirty node up to root
    set<PhyloNode*> ancestors;
    for (PhyloNode* node : dirty) {
        PhyloNode* cur = getParent(node);
        while (cur) {
            if (ancestors.count(cur)) break;  // already traced this path
            ancestors.insert(cur);
            cur = getParent(cur);
        }
    }
    dirty.insert(ancestors.begin(), ancestors.end());
    return dirty;
}

int SPROptimizer::optimizeAtRadius(int radius, int known_score, double wall_seconds) {
    int cur_score = (known_score > 0) ? known_score : tree->fitchRecomputeWithDiffs();
    int initial_score = cur_score;
    int max_id = tree->fitchMaxNodeId();

    auto loop_start = high_resolution_clock::now();
    auto wall_exceeded = [&]() {
        if (wall_seconds <= 0.0) return false;
        double el = duration_cast<milliseconds>(
            high_resolution_clock::now() - loop_start).count() / 1000.0;
        return el > wall_seconds;
    };

    // Initialize bitset for pattern membership
    bitsetInit(tree->fitchNumPatterns());

    string radius_str = (radius == 0) ? "unbounded" : to_string(radius);
    cout << "=== Batch SPR (radius " << radius_str << ") ===" << endl;
    cout << "  Starting score: " << cur_score << endl;

    vector<bool> dirty_nodes(max_id + 1, false);
    bool has_dirty = false;
    int dirty_check_radius = (radius == 0) ? 32 : radius;

    vector<PhyloNode*> all_nodes = collectAllNodes(tree, max_id);

    vector<bool> dfs_visited(max_id + 1, false);
    vector<int> dfs_visited_ids;
    dfs_visited_ids.reserve(max_id + 1);

    for (int round = 0; round < MAX_ROUNDS_PER_RADIUS; round++) {
        if (wall_exceeded()) {
            cout << "  Wall budget (" << wall_seconds
                 << "s) reached after " << round << " rounds; stopping" << endl;
            break;
        }
        vector<SPRCandidate> candidates;
        int moves_evaluated = 0, src_skipped = 0, moves_pruned = 0;
        int wall_break_in_src_loop = 0;
        bool dfs_wall_hit = false;

        for (PhyloNode* src : all_nodes) {
            if (wall_exceeded()) { wall_break_in_src_loop = 1; break; }
            if (src == tree->root || src->degree() < 2) continue;

            PhyloNode* src_parent = getParent(src);
            if (!src_parent || src_parent->degree() != BINARY_NODE_DEGREE) continue;

            if (has_dirty &&
                !isNeighborhoodDirty(src, dirty_check_radius, dirty_nodes, max_id)) {
                src_skipped++;
                continue;
            }

            PhyloNode* sibling1 = NULL, *sibling2 = NULL;
            FOR_NEIGHBOR_IT(src_parent, src, it) {
                PhyloNode* neighbor = (PhyloNode*)(*it)->node;
                if (!sibling1) sibling1 = neighbor; else if (!sibling2) sibling2 = neighbor;
            }
            if (!sibling1 || !sibling2) continue;

            SPRSourceState state(tree, src, src_parent);
            if (!state.isValid()) continue;

            if (!state.hasImprovementPotential()) continue;

            for (int id : dfs_visited_ids) dfs_visited[id] = false;
            dfs_visited_ids.clear();
            dfs_visited[src->id] = true;          dfs_visited_ids.push_back(src->id);
            dfs_visited[src_parent->id] = true;   dfs_visited_ids.push_back(src_parent->id);

            state.setupBaseBitset();

            DFSContext dfs_ctx = {&state, tree, src, src_parent, sibling1, sibling2,
                                  &candidates, &dfs_visited, &dfs_visited_ids,
                                  radius, &moves_evaluated, &moves_pruned,
                                  0, tree->fitchMajorArrayFor(src),
                                  (wall_seconds > 0.0),
                                  loop_start + std::chrono::milliseconds(
                                      (long long)(wall_seconds * 1000.0)),
                                  &dfs_wall_hit};
            searchDestinations(dfs_ctx, sibling1, src_parent, 0, 0);
            searchDestinations(dfs_ctx, sibling2, src_parent, 0, 0);
            if (dfs_wall_hit) {
                wall_break_in_src_loop = 1;
                bitsetClear();
                break;
            }

            // Clear bitset fully for next source
            bitsetClear();
        }

        if (candidates.empty()) {
            cout << "  Round " << round + 1 << ": no improving moves ("
                 << moves_evaluated << " eval";
            if (moves_pruned > 0) cout << ", " << moves_pruned << " pruned";
            if (src_skipped > 0) cout << ", " << src_skipped << " skip";
            if (wall_break_in_src_loop) cout << " [WALL_BREAK]";
            cout << ")" << endl;
            break;
        }
        // If wall budget hit mid-src-loop, apply what we collected and exit.
        bool break_after_round = (wall_break_in_src_loop != 0);

        vector<SPRCandidate> deferred;
        vector<SPRCandidate> selected = selectMoves(candidates, max_id, &deferred);

        vector<NeighborSave> saves;
        for (const auto& move : selected) {
            saveTopology(move.src_parent, saves);
            saveTopology(move.sibling1, saves);
            saveTopology(move.sibling2, saves);
            saveTopology(move.dst, saves);
            saveTopology(move.dst_parent, saves);
        }

        int expected = 0;
        for (const auto& move : selected)
            { applySPRMove(move.src_parent, move.sibling1, move.sibling2, move.dst, move.dst_parent); expected += move.delta; }

        set<PhyloNode*> dirty_set = collectDirtyNodes(selected);
        int new_score = tree->fitchRecomputeScoreDirty(dirty_set);

        cout << "  Round " << round + 1 << ": " << candidates.size() << " found, "
             << selected.size() << " applied";
        if (!deferred.empty()) cout << " (" << deferred.size() << " deferred)";
        cout << " (exp=" << expected
             << " act=" << (new_score - cur_score) << ") -> " << new_score << endl;

        if (new_score > cur_score) {
            undoTopology(saves);
            tree->fitchRecomputeWithDiffs();
            orientTreeToRoot(tree, max_id);
            SPRDeltaExact::precomputeDepths(tree);
            cout << "  Score worsened - reverted" << endl;
            break;
        }

        if (new_score == cur_score) {
            undoTopology(saves);
            tree->fitchRecomputeWithDiffs();
            orientTreeToRoot(tree, max_id);
            SPRDeltaExact::precomputeDepths(tree);
            cout << "  No improvement - reverted" << endl;
            break;
        }

        set<PhyloNode*> dirty_with_ancestors = collectDirtyNodesWithAncestors(selected);
        tree->fitchUpdateDiffsDirty(dirty_with_ancestors);
        orientTreeToRoot(tree, max_id);
        SPRDeltaExact::precomputeDepths(tree);

        markDirty(selected, dirty_nodes, (int)all_nodes.size());
        has_dirty = true;

        if (break_after_round) {
            cur_score = new_score;
            cout << "  Wall budget exceeded mid-round, finishing this round and stopping" << endl;
            break;
        }

        cur_score = new_score;

        // Re-evaluate any conflict-deferred moves under the updated topology
        // and re-admit those that still strictly improve.
        if (!deferred.empty()) {
            int recycled = 0;
            for (auto& dm : deferred) {
                if (!dm.src || !dm.src_parent || !dm.dst || !dm.dst_parent) continue;
                if (dm.src_parent->degree() != BINARY_NODE_DEGREE) continue;

                PhyloNode* sibling1 = NULL, *sibling2 = NULL;
                FOR_NEIGHBOR_IT(dm.src_parent, dm.src, it) {
                    PhyloNode* neighbor = (PhyloNode*)(*it)->node;
                    if (!sibling1) sibling1 = neighbor; else if (!sibling2) sibling2 = neighbor;
                }
                if (!sibling1 || !sibling2) continue;

                SPRSourceState state(tree, dm.src, dm.src_parent);
                if (!state.isValid()) continue;

                state.setupBaseBitset();
                int new_delta = state.evaluate(dm.dst, dm.dst_parent);
                bitsetClear();
                if (new_delta < 0) {
                    dm.sibling1 = sibling1;
                    dm.sibling2 = sibling2;
                    dm.delta = new_delta;
                    candidates.clear();
                    candidates.push_back(dm);
                    recycled++;
                }
            }
            if (recycled > 0) {
                cout << "  Recycled " << recycled << " deferred moves" << endl;
            }
        }
    }

    cout << "  Result: " << initial_score << " -> " << cur_score
         << " (delta=" << (initial_score - cur_score) << ")" << endl;
    return cur_score;
}

int SPROptimizer::optimizeTree(const SPROptimizeOptions& opts) {
    int max_passes    = opts.max_passes;
    int max_radius    = opts.max_radius;
    int ratchet_iters = opts.ratchet_iters;
    int ratchet_seed  = opts.ratchet_seed;
    int ratchet_runs  = opts.ratchet_runs;
    double wall_seconds = opts.wall_seconds;

    if (ratchet_runs < 1) ratchet_runs = 1;
    setActiveSPRTree(tree);

    auto t0 = high_resolution_clock::now();
    int initial_score = tree->fitchRunForSPR();
    auto t1 = high_resolution_clock::now();

    cout << "Fitch: score=" << initial_score
         << " (max_passes=" << max_passes
         << ", max_radius=" << max_radius
         << ", init=" << duration_cast<milliseconds>(t1-t0).count() << "ms)" << endl;

    int tracked_score = initial_score;
    int max_id = tree->fitchMaxNodeId();

    orientTreeToRoot(tree, max_id);
    SPRDeltaExact::precomputeDepths(tree);

    auto optimizer_total_start = high_resolution_clock::now();
    auto optimizer_total_elapsed = [&]() {
        return duration_cast<milliseconds>(
            high_resolution_clock::now() - optimizer_total_start).count() / 1000.0;
    };
    auto wall_exceeded = [&]() {
        return wall_seconds > 0.0 && optimizer_total_elapsed() > wall_seconds;
    };

    // Main pass loop: radius-escalating SPR with strict improvement.
    for (int pass = 0; pass < max_passes; pass++) {
        int start = tracked_score;
        int consecutive_empty = 0;
        int cur = start;

        for (int r = 1; ; r *= 2) {
            r = min(r, max_radius);

            int after = optimizeAtRadius(r, cur);

            if (after >= cur) {
                consecutive_empty++;
                if (consecutive_empty >= RADIUS_STALL_LIMIT && r < max_radius) {
                    cout << "  " << consecutive_empty << " consecutive radii with no improvement, skipping higher" << endl;
                    break;
                }
            } else {
                consecutive_empty = 0;
                cur = after;
            }

            if (r == max_radius) break;
        }

        int end = cur;
        tracked_score = end;
        double improvement = (start > 0) ? (double)(start - end) / start : 0.0;
        cout << "Pass " << pass + 1 << ": " << start << " -> " << end
             << " (delta=" << (start - end)
             << ", improvement=" << fixed << setprecision(4) << (improvement * 100) << "%)" << endl;

        if (end >= start) break;
        if (improvement < CONVERGENCE_THRESHOLD) {
            cout << "  Converged (improvement < " << (CONVERGENCE_THRESHOLD * 100) << "%)" << endl;
            break;
        }
    }

    // Optional parsimony ratchet (Nixon 1999) phase.
    if (ratchet_iters > 0 && !wall_exceeded()) {
        cout << "\n=== Ratchet phase (" << ratchet_iters << " iters"
             << (ratchet_runs > 1 ? " x " + std::to_string(ratchet_runs) + " runs" : "")
             << ", seed=" << ratchet_seed << ") ===" << endl;
        int pre_ratchet = tracked_score;

        std::vector<int> orig_freq;
        tree->fitchSnapshotPatternFreq(orig_freq);
        int nptn = (int)orig_freq.size();

        vector<PhyloNode*> all_nodes_for_save = collectAllNodes(tree, max_id);
        auto saveAll = [&](vector<NeighborSave>& out) {
            out.clear();
            out.reserve(all_nodes_for_save.size() * 3);
            for (PhyloNode* n : all_nodes_for_save) saveTopology(n, out);
        };

        vector<NeighborSave> pre_ratchet_topology;
        saveAll(pre_ratchet_topology);
        std::vector<nuc_one_hot> pre_ratchet_major = tree->fitchSaveNodeMajor();
        PhyloTree::FitchAuxSnapshot pre_ratchet_aux = tree->fitchSaveAux();

        vector<NeighborSave> overall_best_topology = pre_ratchet_topology;
        std::vector<nuc_one_hot> overall_best_major = pre_ratchet_major;
        PhyloTree::FitchAuxSnapshot overall_best_aux = pre_ratchet_aux;
        int overall_best_score = tracked_score;

        for (int run = 0; run < ratchet_runs; run++) {
            if (ratchet_runs > 1 && run > 0) {
                cout << "-- Run " << (run + 1) << "/" << ratchet_runs << " --" << endl;
                undoTopology(pre_ratchet_topology);
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
                    if ((rng() & 3u) == 0u) {
                        tree->fitchScalePatternFreq(ptn, 2);
                        reweighted++;
                    }
                }

                auto remaining_budget = [&]() -> double {
                    if (wall_seconds <= 0.0) return 0.0;
                    double r = wall_seconds - optimizer_total_elapsed();
                    return (r < 0.5) ? 0.5 : r;
                };

                int weighted_score = tree->fitchRecomputeWithDiffs();
                orientTreeToRoot(tree, max_id);
                SPRDeltaExact::precomputeDepths(tree);
                bool capped = false;
                for (int r = 1; ; r *= 2) {
                    r = min(r, max_radius);
                    int after = optimizeAtRadius(r, weighted_score, remaining_budget());
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
                    int after = optimizeAtRadius(r, unweighted_score, remaining_budget());
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
                    accept = true;
                } else if (delta == 0) {
                    accept = ((int)(rng() % 100u) < RATCHET_ZERO_ACCEPT_PCT);
                } else {
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
                    undoTopology(best_saves);
                    tree->fitchRecomputeWithDiffs();
                    orientTreeToRoot(tree, max_id);
                    SPRDeltaExact::precomputeDepths(tree);
                    tracked_score = best_score;
                    tree_at_best = true;
                }
            }
            if (!tree_at_best) {
                undoTopology(best_saves);
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
            undoTopology(overall_best_topology);
            tree->fitchRestoreNodeMajor(overall_best_major);
            tree->fitchRestoreAux(overall_best_aux);
            orientTreeToRoot(tree, max_id);
            SPRDeltaExact::precomputeDepths(tree);
            tracked_score = overall_best_score;
            cout << "Ratchet best across " << ratchet_runs << " runs: "
                << pre_ratchet << " -> " << overall_best_score
                << " (delta=" << (pre_ratchet - overall_best_score) << ")" << endl;
        }
    }

    // Use tracked score — caller will do its own computeParsimony verification
    current_parsimony_score = tracked_score;
    best_score_seen = current_parsimony_score;

    cout << "\nSPR complete: " << initial_score << " -> " << current_parsimony_score
         << " (total delta=" << (initial_score - current_parsimony_score) << ")" << endl;

    setActiveSPRTree(nullptr);
    return current_parsimony_score;
}
