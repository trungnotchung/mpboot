#include "sproptimize.h"
#include "fitch.h"
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
#include <functional>

using namespace std;
using namespace std::chrono;

static const int SPR_MAX_RADIUS = 16;
static const int DEFAULT_MAX_PASSES = 1;
static const int MAX_ROUNDS_PER_RADIUS = 100;
static const int BINARY_NODE_DEGREE = 3;

static PhyloNode* s_root = nullptr;

static inline PhyloNode* getParent(PhyloNode* node) {
    if (!node || node == s_root) return nullptr;
    return (PhyloNode*)node->neighbors[0]->node;
}

static inline void collectDiffsMerged(Fitch* fitch, PhyloNode* node, vector<int>& out) {
    mergeDiffsSorted(out, fitch->getFitchDiffs(node));
}

// Precomputed per-source state for O(M) delta evaluation.
// Scratch buffers are reused across evaluate() calls to avoid heap allocation.
class SPRSourceState {
public:
    SPRSourceState(Fitch* fitch, PhyloNode* src, PhyloNode* src_parent);

    int evaluate(PhyloNode* dst, PhyloNode* dst_parent) const;
    bool isValid() const { return valid; }

private:
    Fitch* fitch;
    PhyloNode* src;
    PhyloNode* src_parent;
    PhyloNode* grandparent;
    PhyloNode* sibling;
    bool valid;

    int num_patterns;

    const nuc_one_hot* src_states;
    const nuc_one_hot* sibling_states;
    const nuc_one_hot* src_parent_states;
    const nuc_one_hot* grandparent_states;

    PhyloNode* gp_sibling;      // other child of grandparent (not src_parent)
    const nuc_one_hot* gp_sibling_states;

    // Path from grandparent to root. Also used for Case A propagation above gp
    // (starting from index 1, since index 0 is grandparent itself).
    vector<PathStep> src_to_root;
    vector<int> base_diffs;

    // Scratch buffers (mutable to allow reuse in const evaluate methods)
    mutable vector<PathStep> scratch_path;
    mutable vector<int> scratch_affected;

    int evalCaseA(PhyloNode* dst, PhyloNode* dst_parent) const;
    int evalCaseB1(PhyloNode* dst) const;
    int evalCaseB2(PhyloNode* dst, PhyloNode* dst_parent, PhyloNode* lca) const;
};

SPRSourceState::SPRSourceState(Fitch* fitch, PhyloNode* src, PhyloNode* src_parent)
    : fitch(fitch), src(src), src_parent(src_parent), valid(false) {
    num_patterns = fitch->getNumPatterns();

    grandparent = getParent(src_parent);
    sibling = findOtherChild(src_parent, grandparent, src);
    if (!sibling || !grandparent) return;
    valid = true;

    src_states = fitch->getMajorArrayForNode(src);
    sibling_states = fitch->getMajorArrayForNode(sibling);
    src_parent_states = fitch->getMajorArrayForNode(src_parent);
    grandparent_states = fitch->getMajorArrayForNode(grandparent);

    gp_sibling = findOtherChild(grandparent, getParent(grandparent), src_parent);
    gp_sibling_states = gp_sibling ? fitch->getMajorArrayForNode(gp_sibling) : nullptr;

    // Build src_to_root path: grandparent → parent(grandparent) → ... → root
    // For Case B: used directly (gp is index 0)
    // For Case A: propagation above gp starts at index 1 (gp's parent)
    {
        PhyloNode* prev = src_parent;
        PhyloNode* cur = grandparent;
        while (cur) {
            PhyloNode* cur_parent = getParent(cur);
            PhyloNode* other = findOtherChild(cur, cur_parent, prev);
            const nuc_one_hot* other_states = other
                ? fitch->getMajorArrayForNode(other)
                : fitch->getMajorArrayForNode(cur);
            src_to_root.push_back({cur, other_states, fitch->getMajorArrayForNode(cur)});
            prev = cur;
            cur = cur_parent;
        }
    }

    collectDiffsMerged(fitch, src, base_diffs);
    collectDiffsMerged(fitch, sibling, base_diffs);
    collectDiffsMerged(fitch, src_parent, base_diffs);
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
    scratch_path.clear();
    {
        PhyloNode* prev = dst;
        PhyloNode* cur = dst_parent;
        while (cur && cur != src_parent) {
            PhyloNode* cur_parent = getParent(cur);
            PhyloNode* other = findOtherChild(cur, cur_parent, prev);
            const nuc_one_hot* other_states = other
                ? fitch->getMajorArrayForNode(other)
                : fitch->getMajorArrayForNode(cur);
            scratch_path.push_back({cur, other_states, fitch->getMajorArrayForNode(cur)});
            prev = cur;
            cur = cur_parent;
        }
    }

    // Collect affected patterns using scratch buffer
    scratch_affected.clear();
    scratch_affected.insert(scratch_affected.end(), base_diffs.begin(), base_diffs.end());
    collectDiffsMerged(fitch, dst, scratch_affected);
    for (const auto& step : scratch_path) collectDiffsMerged(fitch, step.node, scratch_affected);

    const nuc_one_hot* dst_states = fitch->getMajorArrayForNode(dst);
    int score = 0;

    for (int ptn : scratch_affected) {
        int freq = fitch->getPatternFreq(ptn);
        nuc_one_hot src_fitch = src_states[ptn], sibling_fitch = sibling_states[ptn];
        nuc_one_hot dst_fitch = dst_states[ptn], sp_fitch = src_parent_states[ptn];

        // 1. src_parent penalty change
        score += (((src_fitch & dst_fitch) ? 0 : 1) - ((src_fitch & sibling_fitch) ? 0 : 1)) * freq;

        // 2. Propagate dst_parent → sibling
        nuc_one_hot old_fitch = dst_fitch;
        nuc_one_hot new_fitch = fitchMerge(src_fitch, dst_fitch);
        propagatePath(scratch_path, ptn, freq, old_fitch, new_fitch, score);

        nuc_one_hot new_sibling_fitch = (!scratch_path.empty() && old_fitch != new_fitch)
            ? new_fitch : sibling_fitch;

        // 3. At grandparent
        if (grandparent->isLeaf()) {
            score += rootEdgeDelta(grandparent_states[ptn], sp_fitch, new_sibling_fitch) * freq;
        } else {
            nuc_one_hot gp_sibling_fitch = gp_sibling_states ? gp_sibling_states[ptn] : 0xF;
            nuc_one_hot new_gp_fitch;
            score += penaltyDelta(sp_fitch, new_sibling_fitch, gp_sibling_fitch, new_gp_fitch) * freq;

            if (new_gp_fitch != grandparent_states[ptn]) {
                nuc_one_hot old_prop = grandparent_states[ptn], new_prop = new_gp_fitch;
                // Propagate above gp: start at index 1 (gp's parent), skip gp itself
                propagatePath(src_to_root, 1, src_to_root.size(), ptn, freq, old_prop, new_prop, score);
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
        ? src_parent : src_to_root[src_path_len - 1].node;
    const nuc_one_hot* src_branch_states = fitch->getMajorArrayForNode(src_branch);

    PhyloNode* lca_parent = getParent(lca);
    PhyloNode* lca_sibling = findOtherChild(lca, lca_parent, src_branch);
    const nuc_one_hot* lca_sibling_states = lca_sibling ? fitch->getMajorArrayForNode(lca_sibling) : nullptr;
    const nuc_one_hot* lca_states = fitch->getMajorArrayForNode(lca);

    // Collect affected patterns using scratch buffer
    scratch_affected.clear();
    scratch_affected.insert(scratch_affected.end(), base_diffs.begin(), base_diffs.end());
    for (size_t i = 0; i < src_path_len; i++) collectDiffsMerged(fitch, src_to_root[i].node, scratch_affected);
    collectDiffsMerged(fitch, lca, scratch_affected);

    int score = 0;

    for (int ptn : scratch_affected) {
        int freq = fitch->getPatternFreq(ptn);
        nuc_one_hot src_fitch = src_states[ptn], sibling_fitch = sibling_states[ptn];
        nuc_one_hot sp_fitch = src_parent_states[ptn];
        int old_sp_penalty = (src_fitch & sibling_fitch) ? 0 : 1;

        // 1. src-side propagation: sp_fitch → sibling_fitch
        nuc_one_hot src_old = sp_fitch, src_new = sibling_fitch;
        propagatePath(src_to_root, 0, src_path_len, ptn, freq, src_old, src_new, score);

        nuc_one_hot src_child_old = src_branch_states[ptn];
        nuc_one_hot src_child_new = (src_old != src_new) ? src_new : src_child_old;

        // 2. New src_parent penalty
        score += (((src_fitch & src_child_new) ? 0 : 1) - old_sp_penalty) * freq;

        // 3. New src_parent Fitch state
        nuc_one_hot new_sp_fitch = fitchMerge(src_fitch, src_child_new);

        // 4. At lca: child changed from src_child_old to new_sp_fitch
        nuc_one_hot lca_sibling_fitch = lca_sibling_states ? lca_sibling_states[ptn] : 0xF;
        nuc_one_hot old_lca_fitch = lca_states[ptn];
        nuc_one_hot new_lca_fitch;
        score += penaltyDelta(src_child_old, new_sp_fitch, lca_sibling_fitch, new_lca_fitch) * freq;

        // 5. Propagate above lca
        if (!lca_parent && lca->isLeaf()) {
            score += rootEdgeDelta(lca_states[ptn], src_child_old, new_sp_fitch) * freq;
        } else if (new_lca_fitch != old_lca_fitch && lca_parent) {
            nuc_one_hot old_prop = old_lca_fitch, new_prop = new_lca_fitch;
            if (lca_parent->isLeaf()) {
                const nuc_one_hot* lca_parent_states = fitch->getMajorArrayForNode(lca_parent);
                score += rootEdgeDelta(lca_parent_states[ptn], old_lca_fitch, new_lca_fitch) * freq;
            } else {
                propagatePath(src_to_root, src_path_len + 1, src_to_root.size(), ptn, freq, old_prop, new_prop, score);
            }
        }
    }

    return score;
}

// ========== CASE B2: general ==========
int SPRSourceState::evalCaseB2(PhyloNode* dst, PhyloNode* dst_parent, PhyloNode* lca) const {
    // Find src_path_len: grandparent → ... → node before lca
    size_t src_path_len = 0;
    for (size_t i = 0; i < src_to_root.size(); i++) {
        if (src_to_root[i].node == lca) { src_path_len = i; break; }
        src_path_len = i + 1;
    }

    PhyloNode* src_branch = (src_path_len == 0)
        ? src_parent : src_to_root[src_path_len - 1].node;
    const nuc_one_hot* src_branch_states = fitch->getMajorArrayForNode(src_branch);

    // Build dst path using scratch buffer (dst_parent → ... → lca, exclusive)
    scratch_path.clear();
    {
        PhyloNode* prev = dst;
        PhyloNode* cur = dst_parent;
        while (cur && cur != lca) {
            PhyloNode* cur_parent = getParent(cur);
            PhyloNode* other = findOtherChild(cur, cur_parent, prev);
            const nuc_one_hot* other_states = other
                ? fitch->getMajorArrayForNode(other)
                : fitch->getMajorArrayForNode(cur);
            scratch_path.push_back({cur, other_states, fitch->getMajorArrayForNode(cur)});
            prev = cur;
            cur = cur_parent;
        }
    }

    PhyloNode* dst_branch = scratch_path.empty() ? dst : scratch_path.back().node;
    const nuc_one_hot* dst_branch_states = fitch->getMajorArrayForNode(dst_branch);
    const nuc_one_hot* lca_states = fitch->getMajorArrayForNode(lca);
    const nuc_one_hot* dst_states = fitch->getMajorArrayForNode(dst);

    // Collect affected patterns using scratch buffer
    scratch_affected.clear();
    scratch_affected.insert(scratch_affected.end(), base_diffs.begin(), base_diffs.end());
    collectDiffsMerged(fitch, dst, scratch_affected);
    for (size_t i = 0; i < src_path_len; i++) collectDiffsMerged(fitch, src_to_root[i].node, scratch_affected);
    for (const auto& step : scratch_path) collectDiffsMerged(fitch, step.node, scratch_affected);

    int score = 0;

    for (int ptn : scratch_affected) {
        int freq = fitch->getPatternFreq(ptn);
        nuc_one_hot src_fitch = src_states[ptn], sibling_fitch = sibling_states[ptn];
        nuc_one_hot dst_fitch = dst_states[ptn], sp_fitch = src_parent_states[ptn];

        // 1. src_parent penalty change
        score += (((src_fitch & dst_fitch) ? 0 : 1) - ((src_fitch & sibling_fitch) ? 0 : 1)) * freq;

        // 2. src-side: sp_fitch → sibling_fitch
        nuc_one_hot src_old = sp_fitch, src_new = sibling_fitch;
        propagatePath(src_to_root, 0, src_path_len, ptn, freq, src_old, src_new, score);

        // 3. dst-side: dst_fitch → new_sp_fitch
        nuc_one_hot new_sp_fitch = fitchMerge(src_fitch, dst_fitch);
        nuc_one_hot dst_old = dst_fitch, dst_new = new_sp_fitch;
        propagatePath(scratch_path, ptn, freq, dst_old, dst_new, score);

        // 4. At LCA: combine changes from both sides
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

        // 5. Propagate above LCA
        if (new_lca_fitch != old_lca_fitch) {
            nuc_one_hot old_prop = old_lca_fitch, new_prop = new_lca_fitch;
            propagatePath(src_to_root, src_path_len + 1, src_to_root.size(), ptn, freq, old_prop, new_prop, score);
        }
    }

    return score;
}

struct SPRCandidate {
    PhyloNode* src;
    PhyloNode* src_parent;
    PhyloNode* sibling1;
    PhyloNode* sibling2;
    PhyloNode* dst;
    PhyloNode* dst_parent;
    int delta;
};

struct NeighborSave { Neighbor* neighbor; Node* original_node; };

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
        FOR_NEIGHBOR_IT(node, nullptr, it) {
            PhyloNode* neighbor = (PhyloNode*)(*it)->node;
            if (!visited[neighbor->id]) q.push(neighbor);
        }
    }
    return nodes;
}

static void orientTreeToRoot(PhyloTree* tree, int max_id) {
    PhyloNode* root = (PhyloNode*)tree->root;

    vector<bool> visited(max_id + 1, false);
    queue<pair<PhyloNode*, PhyloNode*>> q;
    q.push({root, nullptr});

    while (!q.empty()) {
        PhyloNode* node = (PhyloNode*)q.front().first; PhyloNode* parent = (PhyloNode*)q.front().second; q.pop();
        if (visited[node->id]) continue;
        visited[node->id] = true;

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
            if (!visited[child->id]) q.push({child, node});
        }
    }

    s_root = root;
    SPRMutationOps::setRoot(root);
}

static void applySPRMove(PhyloNode* src_parent, PhyloNode* sibling1, PhyloNode* sibling2,
                         PhyloNode* dst, PhyloNode* dst_parent) {
    if (dst_parent == sibling1) {
        sibling1->updateNeighbor(src_parent, sibling2);  sibling2->updateNeighbor(src_parent, sibling1);
        sibling1->updateNeighbor(dst, src_parent);       src_parent->updateNeighbor(sibling2, dst);
        dst->updateNeighbor(sibling1, src_parent);
    } else if (dst_parent == sibling2) {
        sibling2->updateNeighbor(src_parent, sibling1);  sibling1->updateNeighbor(src_parent, sibling2);
        sibling2->updateNeighbor(dst, src_parent);       src_parent->updateNeighbor(sibling1, dst);
        dst->updateNeighbor(sibling2, src_parent);
    } else if (dst == sibling1) {
        dst->updateNeighbor(src_parent, sibling2);       sibling2->updateNeighbor(src_parent, dst);
        dst->updateNeighbor(dst_parent, src_parent);     src_parent->updateNeighbor(sibling2, dst_parent);
        dst_parent->updateNeighbor(dst, src_parent);
    } else if (dst == sibling2) {
        dst->updateNeighbor(src_parent, sibling1);       sibling1->updateNeighbor(src_parent, dst);
        dst->updateNeighbor(dst_parent, src_parent);     src_parent->updateNeighbor(sibling1, dst_parent);
        dst_parent->updateNeighbor(dst, src_parent);
    } else {
        sibling1->updateNeighbor(src_parent, sibling2);  sibling2->updateNeighbor(src_parent, sibling1);
        src_parent->updateNeighbor(sibling1, dst);       src_parent->updateNeighbor(sibling2, dst_parent);
        dst->updateNeighbor(dst_parent, src_parent);     dst_parent->updateNeighbor(dst, src_parent);
    }
}

static void saveTopology(Node* node, vector<NeighborSave>& saves) {
    for (auto it = node->neighbors.begin(); it != node->neighbors.end(); it++)
        saves.push_back({*it, (*it)->node});
}

static void undoTopology(vector<NeighborSave>& saves) {
    for (auto& save : saves) save.neighbor->node = save.original_node;
}

static void markDirty(const vector<SPRCandidate>& moves, vector<bool>& dirty,
                      int max_nodes) {
    fill(dirty.begin(), dirty.end(), false);
    for (const auto& move : moves) {
        dirty[move.src->id] = true;
        dirty[move.src_parent->id] = true;
        dirty[move.sibling1->id] = true;
        dirty[move.sibling2->id] = true;
        dirty[move.dst->id] = true;
        dirty[move.dst_parent->id] = true;
        int walk = 0;
        for (PhyloNode* ancestor = getParent(move.src_parent);
             ancestor && walk < max_nodes; ancestor = getParent(ancestor), walk++)
            dirty[ancestor->id] = true;
        walk = 0;
        for (PhyloNode* ancestor = getParent(move.dst_parent);
             ancestor && walk < max_nodes; ancestor = getParent(ancestor), walk++)
            dirty[ancestor->id] = true;
    }
}

// Reusable scratch buffer for BFS visited tracking (avoids per-call allocation)
static vector<bool> s_bfs_visited;
static vector<int> s_bfs_visited_ids;  // track which IDs were set for fast reset

static bool isNeighborhoodDirty(PhyloNode* node, int radius, const vector<bool>& dirty) {
    // Check if any dirty flag is set (empty check)
    bool has_dirty = false;
    for (size_t i = 0; i < dirty.size() && !has_dirty; i++) has_dirty = dirty[i];
    if (!has_dirty) return true;

    if (dirty[node->id]) return true;

    // Fast BFS with reusable visited buffer
    s_bfs_visited_ids.clear();
    queue<pair<PhyloNode*, int>> q;
    q.push({node, 0});
    s_bfs_visited[node->id] = true;
    s_bfs_visited_ids.push_back(node->id);

    bool found = false;
    while (!q.empty() && !found) {
        PhyloNode* cur_node = q.front().first; int depth = q.front().second; q.pop();
        if (dirty[cur_node->id]) { found = true; break; }
        if (depth < radius) {
            FOR_NEIGHBOR_IT(cur_node, nullptr, it) {
                PhyloNode* neighbor = (PhyloNode*)(*it)->node;
                if (!s_bfs_visited[neighbor->id]) {
                    s_bfs_visited[neighbor->id] = true;
                    s_bfs_visited_ids.push_back(neighbor->id);
                    q.push({neighbor, depth + 1});
                }
            }
        }
    }

    // Reset visited (only the IDs we touched — O(visited) not O(max_id))
    for (int id : s_bfs_visited_ids) s_bfs_visited[id] = false;
    return found;
}

SPROptimizer::SPROptimizer(PhyloTree* tree) : tree(tree), current_parsimony_score(0),
    best_score_seen(0) {
    if (!tree) throw std::invalid_argument("SPROptimizer: tree cannot be NULL");
}

SPROptimizer::~SPROptimizer() {}

// Select non-conflicting moves greedily (best delta first).
static vector<SPRCandidate> selectMoves(vector<SPRCandidate>& candidates, int max_id) {
    sort(candidates.begin(), candidates.end(),
         [](const SPRCandidate& a, const SPRCandidate& b) { return a.delta < b.delta; });

    vector<bool> used(max_id + 1, false);
    vector<SPRCandidate> selected;
    for (const auto& move : candidates) {
        if (used[move.src->id] || used[move.src_parent->id] ||
            used[move.sibling1->id] || used[move.sibling2->id] ||
            used[move.dst->id] || used[move.dst_parent->id])
            continue;
        selected.push_back(move);
        used[move.src->id] = true;         used[move.src_parent->id] = true;
        used[move.sibling1->id] = true;    used[move.sibling2->id] = true;
        used[move.dst->id] = true;         used[move.dst_parent->id] = true;
    }
    return selected;
}

int SPROptimizer::optimizeAtRadius(int radius, Fitch& fitch) {
    int cur_score = fitch.recomputeWithDiffs();
    int initial_score = cur_score;
    int max_id = fitch.getMaxNodeId();

    orientTreeToRoot(tree, max_id);
    SPRDeltaExact::precomputeDepths(tree);

    string radius_str = (radius == 0) ? "unbounded" : to_string(radius);
    cout << "=== Batch SPR (radius " << radius_str << ") ===" << endl;
    cout << "  Starting score: " << cur_score << endl;

    vector<bool> dirty_nodes(max_id + 1, false);
    bool has_dirty = false;
    int dirty_check_radius = (radius == 0) ? 32 : radius;

    // Initialize reusable BFS visited buffer
    s_bfs_visited.assign(max_id + 1, false);
    s_bfs_visited_ids.reserve(max_id + 1);

    // Reusable DFS visited buffer
    vector<bool> dfs_visited(max_id + 1, false);
    vector<int> dfs_visited_ids;
    dfs_visited_ids.reserve(max_id + 1);

    for (int round = 0; round < MAX_ROUNDS_PER_RADIUS; round++) {
        auto start_time = high_resolution_clock::now();
        vector<PhyloNode*> all_nodes = collectAllNodes(tree, max_id);
        vector<SPRCandidate> candidates;
        int moves_evaluated = 0, src_skipped = 0;

        for (PhyloNode* src : all_nodes) {
            if (src == tree->root || src->degree() < 2) continue;

            PhyloNode* src_parent = getParent(src);
            if (!src_parent || src_parent->degree() != BINARY_NODE_DEGREE) continue;

            if (has_dirty &&
                !isNeighborhoodDirty(src, dirty_check_radius, dirty_nodes)) {
                src_skipped++;
                continue;
            }

            PhyloNode* sibling1 = NULL, *sibling2 = NULL;
            FOR_NEIGHBOR_IT(src_parent, src, it) {
                PhyloNode* neighbor = (PhyloNode*)(*it)->node;
                if (!sibling1) sibling1 = neighbor; else if (!sibling2) sibling2 = neighbor;
            }
            if (!sibling1 || !sibling2) continue;

            SPRSourceState state(&fitch, src, src_parent);
            if (!state.isValid()) continue;

            // Reset DFS visited for this source
            for (int id : dfs_visited_ids) dfs_visited[id] = false;
            dfs_visited_ids.clear();
            dfs_visited[src->id] = true;          dfs_visited_ids.push_back(src->id);
            dfs_visited[src_parent->id] = true;   dfs_visited_ids.push_back(src_parent->id);

            function<void(PhyloNode*, PhyloNode*, int)> search_destinations =
                [&](PhyloNode* node, PhyloNode* from, int dist) {
                if (dfs_visited[node->id]) return;
                dfs_visited[node->id] = true;
                dfs_visited_ids.push_back(node->id);

                if (from != src_parent && from != src && node != src && node != src_parent) {
                    int delta = state.evaluate(node, from);
                    moves_evaluated++;
                    if (delta < 0)
                        candidates.push_back({src, src_parent, sibling1, sibling2, node, from, delta});
                }

                if (radius > 0 && dist >= radius) return;
                FOR_NEIGHBOR_IT(node, nullptr, nit) {
                    PhyloNode* next = (PhyloNode*)(*nit)->node;
                    if (!dfs_visited[next->id]) search_destinations(next, node, dist + 1);
                }
            };

            search_destinations(sibling1, src_parent, 0);
            search_destinations(sibling2, src_parent, 0);
        }

        long long elapsed = duration_cast<milliseconds>(high_resolution_clock::now() - start_time).count();

        if (candidates.empty()) {
            cout << "  Round " << round + 1 << ": no improving moves ("
                 << moves_evaluated << " eval";
            if (src_skipped > 0) cout << ", " << src_skipped << " skip";
            cout << ", " << elapsed << "ms)" << endl << flush;
            break;
        }

        vector<SPRCandidate> selected = selectMoves(candidates, max_id);

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

        int new_score = fitch.recomputeWithDiffs();
        orientTreeToRoot(tree, max_id);
        SPRDeltaExact::precomputeDepths(tree);

        cout << "  Round " << round + 1 << ": " << candidates.size() << " found, "
             << selected.size() << " applied (exp=" << expected
             << " act=" << (new_score - cur_score) << ") -> " << new_score
             << " (" << moves_evaluated << " eval";
        if (src_skipped > 0) cout << ", " << src_skipped << " skip";
        cout << ", " << elapsed << "ms)" << endl << flush;

        if (new_score >= cur_score) {
            undoTopology(saves);
            fitch.recomputeWithDiffs();
            orientTreeToRoot(tree, max_id);
            SPRDeltaExact::precomputeDepths(tree);
            cout << "  No improvement - reverted" << endl;
            break;
        }

        markDirty(selected, dirty_nodes, (int)all_nodes.size());
        has_dirty = true;
        cur_score = new_score;
        SPRDeltaExact::setCustomFitch(&fitch, cur_score);
    }

    cout << "  Result: " << initial_score << " -> " << cur_score
         << " (delta=" << (initial_score - cur_score) << ")" << endl;
    return cur_score;
}

int SPROptimizer::optimizeTree(int max_passes) {
    Fitch fitch(tree);
    int initial_score = fitch.run();
    SPRMutationOps::setCustomFitch(&fitch);
    SPRDeltaExact::setCustomFitch(&fitch, initial_score);

    cout << "Fitch: score=" << initial_score
         << " (max_passes=" << max_passes << ")" << endl;

    for (int pass = 0; pass < max_passes; pass++) {
        int start = fitch.recomputeWithDiffs();
        SPRDeltaExact::setCustomFitch(&fitch, start);

        for (int r = 1; r <= SPR_MAX_RADIUS; r *= 2)
            optimizeAtRadius(r, fitch);
        optimizeAtRadius(0, fitch);

        int end = fitch.recomputeWithDiffs();
        cout << "Pass " << pass + 1 << ": " << start << " -> " << end
             << " (delta=" << (start - end) << ")" << endl;
        if (end >= start) break;
    }

    current_parsimony_score = fitch.recompute();
    best_score_seen = current_parsimony_score;

    cout << "\nSPR complete: " << initial_score << " -> " << current_parsimony_score
         << " (total delta=" << (initial_score - current_parsimony_score) << ")" << endl;
    return current_parsimony_score;
}
