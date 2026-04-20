#include "spr_delta_exact.h"
#include "spr_mutation_ops.h"
#include "spr_utils.h"
#include "fitch.h"
#include <algorithm>
#include <map>
#include <queue>

// Static state for binary fallback
static Fitch* s_fitch = nullptr;
static int s_current_score = 0;

static std::vector<int> s_node_depth;

// Euler tour + sparse table for O(1) LCA queries
static std::vector<PhyloNode*> s_euler;
static std::vector<int> s_euler_depth;
static std::vector<int> s_first_appearance;
static std::vector<std::vector<int>> s_sparse_table;
static std::vector<int> s_log2;
static bool s_lca_ready = false;

void SPRDeltaExact::setCustomFitch(Fitch* fitch, int current_score) {
    s_fitch = fitch;
    s_current_score = current_score;
}

static void eulerTourDFS(PhyloNode* root, int max_id) {
    int estimated_nodes = max_id + 1;
    s_euler.clear();
    s_euler.reserve(estimated_nodes * 2);
    s_euler_depth.clear();
    s_euler_depth.reserve(estimated_nodes * 2);
    s_first_appearance.assign(max_id + 1, -1);
    s_node_depth.assign(max_id + 1, -1);

    struct DFSFrame {
        PhyloNode* node;
        PhyloNode* parent;
        int child_idx;       // which child to process next
        int depth;
        int num_children;
        int children_offset; // offset into children_buf
    };

    static std::vector<PhyloNode*> children_buf;
    static std::vector<DFSFrame> stack;
    children_buf.clear();
    children_buf.reserve(estimated_nodes * 3);  // each node has ~3 neighbors
    stack.clear();
    stack.reserve(estimated_nodes);

    DFSFrame root_frame;
    root_frame.node = root;
    root_frame.parent = nullptr;
    root_frame.child_idx = 0;
    root_frame.depth = 0;

    int children_start = (int)children_buf.size();
    FOR_NEIGHBOR_IT(root, nullptr, it) {
        children_buf.push_back((PhyloNode*)(*it)->node);
    }
    root_frame.num_children = (int)children_buf.size() - children_start;
    root_frame.children_offset = children_start;

    s_euler.push_back(root);
    s_euler_depth.push_back(0);
    s_first_appearance[root->id] = 0;
    s_node_depth[root->id] = 0;

    stack.push_back(root_frame);

    while (!stack.empty()) {
        DFSFrame& frame = stack.back();

        if (frame.child_idx < frame.num_children) {
            PhyloNode* child = children_buf[frame.children_offset + frame.child_idx];
            frame.child_idx++;

            if (child == frame.parent) continue;
            if (child->id >= 0 && child->id <= max_id && s_node_depth[child->id] >= 0) continue;

            int child_depth = frame.depth + 1;

            s_node_depth[child->id] = child_depth;
            s_first_appearance[child->id] = (int)s_euler.size();
            s_euler.push_back(child);
            s_euler_depth.push_back(child_depth);

            DFSFrame child_frame;
            child_frame.node = child;
            child_frame.parent = frame.node;
            child_frame.child_idx = 0;
            child_frame.depth = child_depth;

            int cs = (int)children_buf.size();
            FOR_NEIGHBOR_IT(child, nullptr, it2) {
                children_buf.push_back((PhyloNode*)(*it2)->node);
            }
            child_frame.num_children = (int)children_buf.size() - cs;
            child_frame.children_offset = cs;

            stack.push_back(child_frame);
        } else {
            stack.pop_back();
            if (!stack.empty()) {
                DFSFrame& parent_frame = stack.back();
                s_euler.push_back(parent_frame.node);
                s_euler_depth.push_back(parent_frame.depth);
            }
        }
    }
}

static void buildSparseTable() {
    int n = (int)s_euler_depth.size();
    if (n == 0) return;

    s_log2.assign(n + 1, 0);
    for (int i = 2; i <= n; i++)
        s_log2[i] = s_log2[i / 2] + 1;

    int LOG = s_log2[n] + 1;
    s_sparse_table.assign(n, std::vector<int>(LOG, 0));

    for (int i = 0; i < n; i++)
        s_sparse_table[i][0] = i;

    for (int j = 1; j < LOG; j++) {
        int len = 1 << j;
        for (int i = 0; i + len <= n; i++) {
            int left = s_sparse_table[i][j - 1];
            int right = s_sparse_table[i + (len >> 1)][j - 1];
            s_sparse_table[i][j] = (s_euler_depth[left] <= s_euler_depth[right]) ? left : right;
        }
    }
}

static int rmq(int l, int r) {
    if (l > r) std::swap(l, r);
    int len = r - l + 1;
    int k = s_log2[len];
    int left = s_sparse_table[l][k];
    int right = s_sparse_table[r - (1 << k) + 1][k];
    return (s_euler_depth[left] <= s_euler_depth[right]) ? left : right;
}

void SPRDeltaExact::precomputeDepths(PhyloTree* tree) {
    s_lca_ready = false;
    s_node_depth.clear();
    if (!tree || !tree->root) return;

    int max_id = 0;
    std::queue<PhyloNode*> bfs_q;
    std::vector<bool> visited;
    int initial_size = std::max(1, tree->nodeNum * 2);
    visited.assign(initial_size, false);
    bfs_q.push((PhyloNode*)tree->root);

    while (!bfs_q.empty()) {
        PhyloNode* node = bfs_q.front(); bfs_q.pop();
        if (node->id < 0) continue;
        if (node->id >= (int)visited.size()) visited.resize(node->id + 1, false);
        if (visited[node->id]) continue;
        visited[node->id] = true;
        if (node->id > max_id) max_id = node->id;
        FOR_NEIGHBOR_IT(node, nullptr, it) {
            PhyloNode* nb = (PhyloNode*)(*it)->node;
            if (nb->id >= 0) {
                if (nb->id >= (int)visited.size() || !visited[nb->id])
                    bfs_q.push(nb);
            }
        }
    }

    eulerTourDFS((PhyloNode*)tree->root, max_id);
    buildSparseTable();
    s_lca_ready = true;
}

static std::vector<PathStep> buildPath(PhyloNode* prev, PhyloNode* cur,
                                        PhyloNode* stop_node, Fitch* fitch) {
    std::vector<PathStep> path;
    while (cur && cur != stop_node) {
        PhyloNode* cur_parent = SPRMutationOps::getParent(cur);
        PhyloNode* other = findOtherChild(cur, cur_parent, prev);
        const nuc_one_hot* other_states = other
            ? fitch->getMajorArrayForNode(other)
            : fitch->getMajorArrayForNode(cur);  // root leaf: use own state
        path.push_back({cur, other_states, fitch->getMajorArrayForNode(cur)});
        prev = cur;
        cur = cur_parent;
    }
    return path;
}

static inline void collectDiffs(Fitch* fitch, PhyloNode* node, std::vector<int>& out) {
    mergeDiffsSorted(out, fitch->getFitchDiffs(node));
}

static long long s_total_M = 0, s_total_P = 0;
static int s_move_count = 0;
static int s_max_M = 0;

static int sparseBinaryDelta(PhyloNode* src, PhyloNode* src_parent,
                              PhyloNode* dst, PhyloNode* dst_parent) {
    if (!s_fitch) return 0;
    Fitch* fitch = s_fitch;

    PhyloNode* grandparent = SPRMutationOps::getParent(src_parent);
    PhyloNode* sibling = findOtherChild(src_parent, grandparent, src);
    if (!sibling || !grandparent) return 0;

    PhyloNode* lca = SPRDeltaExact::findLCA(src_parent, dst, nullptr);
    if (!lca) return 0;

    const nuc_one_hot* src_states = fitch->getMajorArrayForNode(src);
    const nuc_one_hot* sibling_st = fitch->getMajorArrayForNode(sibling);
    const nuc_one_hot* sp_states = fitch->getMajorArrayForNode(src_parent);
    const nuc_one_hot* dst_states = fitch->getMajorArrayForNode(dst);

    int score = 0;
    int num_patterns = fitch->getNumPatterns();
    int local_M = 0; // for stats

    if (lca == src_parent) {
        // ========== CASE A: dst in sibling's subtree ==========
        std::vector<PathStep> path_to_sp = buildPath(dst, dst_parent, src_parent, fitch);

        PhyloNode* gp_sibling = findOtherChild(grandparent, SPRMutationOps::getParent(grandparent), src_parent);
        const nuc_one_hot* gp_sibling_states = gp_sibling ? fitch->getMajorArrayForNode(gp_sibling) : nullptr;
        const nuc_one_hot* gp_states = fitch->getMajorArrayForNode(grandparent);

        std::vector<PathStep> path_above_gp = buildPath(
            grandparent, SPRMutationOps::getParent(grandparent), nullptr, fitch);

        // Collect O(M) affected patterns
        std::vector<int> affected;
        collectDiffs(fitch, src, affected);
        collectDiffs(fitch, sibling, affected);
        collectDiffs(fitch, src_parent, affected);
        collectDiffs(fitch, dst, affected);
        for (const auto& step : path_to_sp) collectDiffs(fitch, step.node, affected);
        local_M = (int)affected.size();

        for (int ptn : affected) {
            int freq = fitch->getPatternFreq(ptn);
            nuc_one_hot src_fitch = src_states[ptn], sibling_fitch = sibling_st[ptn];
            nuc_one_hot dst_fitch = dst_states[ptn], sp_fitch = sp_states[ptn];

            // 1. src_parent penalty change
            score += (((src_fitch & dst_fitch) ? 0 : 1) - ((src_fitch & sibling_fitch) ? 0 : 1)) * freq;

            // 2. Propagate dst_parent → sibling
            nuc_one_hot old_fitch = dst_fitch;
            nuc_one_hot new_fitch = fitchMerge(src_fitch, dst_fitch);
            propagatePath(path_to_sp, ptn, freq, old_fitch, new_fitch, score);

            nuc_one_hot new_sibling_fitch = (!path_to_sp.empty() && old_fitch != new_fitch)
                ? new_fitch : sibling_fitch;

            // 3. At grandparent
            if (grandparent->isLeaf()) {
                score += rootEdgeDelta(gp_states[ptn], sp_fitch, new_sibling_fitch) * freq;
            } else {
                nuc_one_hot gp_sibling_fitch = gp_sibling_states ? gp_sibling_states[ptn] : 0xF;
                nuc_one_hot new_gp_fitch;
                score += penaltyDelta(sp_fitch, new_sibling_fitch, gp_sibling_fitch, new_gp_fitch) * freq;

                if (new_gp_fitch != gp_states[ptn]) {
                    nuc_one_hot old_prop = gp_states[ptn], new_prop = new_gp_fitch;
                    propagatePath(path_above_gp, ptn, freq, old_prop, new_prop, score);
                }
            }
        }

    } else if (dst == lca) {
        // ========== CASE B1: dst == lca ==========
        std::vector<PathStep> src_path = buildPath(src_parent, grandparent, lca, fitch);
        PhyloNode* src_branch = src_path.empty() ? src_parent : src_path.back().node;
        const nuc_one_hot* src_branch_states = fitch->getMajorArrayForNode(src_branch);

        PhyloNode* lca_parent = SPRMutationOps::getParent(lca);
        PhyloNode* lca_sibling = findOtherChild(lca, lca_parent, src_branch);
        const nuc_one_hot* lca_sibling_states = lca_sibling ? fitch->getMajorArrayForNode(lca_sibling) : nullptr;
        const nuc_one_hot* lca_states = fitch->getMajorArrayForNode(lca);

        // parent(lca) info — for propagation above lca
        PhyloNode* lca_parent_sibling = lca_parent ? findOtherChild(lca_parent,
            SPRMutationOps::getParent(lca_parent), lca) : nullptr;
        const nuc_one_hot* lca_parent_sibling_states = lca_parent_sibling
            ? fitch->getMajorArrayForNode(lca_parent_sibling) : nullptr;
        const nuc_one_hot* lca_parent_states = lca_parent
            ? fitch->getMajorArrayForNode(lca_parent) : nullptr;

        std::vector<PathStep> above_lca_parent = buildPath(
            lca_parent, lca_parent ? SPRMutationOps::getParent(lca_parent) : nullptr,
            nullptr, fitch);

        // Collect O(M) affected patterns
        std::vector<int> affected;
        collectDiffs(fitch, src, affected);
        collectDiffs(fitch, sibling, affected);
        collectDiffs(fitch, src_parent, affected);
        for (const auto& step : src_path) collectDiffs(fitch, step.node, affected);
        collectDiffs(fitch, lca, affected);
        local_M = (int)affected.size();

        for (int ptn : affected) {
            int freq = fitch->getPatternFreq(ptn);
            nuc_one_hot src_fitch = src_states[ptn], sibling_fitch = sibling_st[ptn];
            nuc_one_hot sp_fitch = sp_states[ptn];
            int old_sp_penalty = (src_fitch & sibling_fitch) ? 0 : 1;

            // 1. src-side propagation
            nuc_one_hot src_old = sp_fitch, src_new = sibling_fitch;
            propagatePath(src_path, ptn, freq, src_old, src_new, score);

            nuc_one_hot src_child_old = src_branch_states[ptn];
            nuc_one_hot src_child_new = (src_old != src_new) ? src_new : src_child_old;

            // 2. New src_parent penalty
            score += (((src_fitch & src_child_new) ? 0 : 1) - old_sp_penalty) * freq;

            // 3. New src_parent Fitch state
            nuc_one_hot new_sp_fitch = fitchMerge(src_fitch, src_child_new);

            // 4. At lca: child changed from src_child_old to new_sp_fitch
            nuc_one_hot lca_sib_fitch = lca_sibling_states ? lca_sibling_states[ptn] : 0xF;
            nuc_one_hot old_lca_fitch = lca_states[ptn];
            nuc_one_hot new_lca_fitch;
            score += penaltyDelta(src_child_old, new_sp_fitch, lca_sib_fitch, new_lca_fitch) * freq;

            // 5. Propagate above lca
            if (!lca_parent && lca->isLeaf()) {
                score += rootEdgeDelta(lca_states[ptn], src_child_old, new_sp_fitch) * freq;
            } else if (new_lca_fitch != old_lca_fitch && lca_parent) {
                if (lca_parent->isLeaf()) {
                    score += rootEdgeDelta(lca_parent_states[ptn], old_lca_fitch, new_lca_fitch) * freq;
                } else {
                    nuc_one_hot lca_parent_sib_fitch = lca_parent_sibling_states
                        ? lca_parent_sibling_states[ptn] : 0xF;
                    nuc_one_hot new_lca_parent_fitch;
                    score += penaltyDelta(old_lca_fitch, new_lca_fitch, lca_parent_sib_fitch,
                                          new_lca_parent_fitch) * freq;

                    if (new_lca_parent_fitch != (lca_parent_states ? lca_parent_states[ptn] : (nuc_one_hot)0xF)) {
                        nuc_one_hot old_prop = lca_parent_states[ptn], new_prop = new_lca_parent_fitch;
                        propagatePath(above_lca_parent, ptn, freq, old_prop, new_prop, score);
                    }
                }
            }
        }

    } else {
        // ========== CASE B2: general ==========
        std::vector<PathStep> src_path = buildPath(src_parent, grandparent, lca, fitch);
        std::vector<PathStep> dst_path = buildPath(dst, dst_parent, lca, fitch);

        PhyloNode* src_branch = src_path.empty() ? src_parent : src_path.back().node;
        PhyloNode* dst_branch = dst_path.empty() ? dst : dst_path.back().node;
        const nuc_one_hot* src_branch_states = fitch->getMajorArrayForNode(src_branch);
        const nuc_one_hot* dst_branch_states = fitch->getMajorArrayForNode(dst_branch);

        const nuc_one_hot* lca_states = fitch->getMajorArrayForNode(lca);
        std::vector<PathStep> above_lca = buildPath(
            lca, SPRMutationOps::getParent(lca), nullptr, fitch);

        // Collect O(M) affected patterns
        std::vector<int> affected;
        collectDiffs(fitch, src, affected);
        collectDiffs(fitch, sibling, affected);
        collectDiffs(fitch, src_parent, affected);
        collectDiffs(fitch, dst, affected);
        for (const auto& step : src_path) collectDiffs(fitch, step.node, affected);
        for (const auto& step : dst_path) collectDiffs(fitch, step.node, affected);
        local_M = (int)affected.size();

        for (int ptn : affected) {
            int freq = fitch->getPatternFreq(ptn);
            nuc_one_hot src_fitch = src_states[ptn], sibling_fitch = sibling_st[ptn];
            nuc_one_hot dst_fitch = dst_states[ptn], sp_fitch = sp_states[ptn];

            // 1. src_parent penalty change
            score += (((src_fitch & dst_fitch) ? 0 : 1) - ((src_fitch & sibling_fitch) ? 0 : 1)) * freq;

            // 2. src-side: sp_fitch → sibling_fitch
            nuc_one_hot src_old = sp_fitch, src_new = sibling_fitch;
            propagatePath(src_path, ptn, freq, src_old, src_new, score);

            // 3. dst-side: dst_fitch → new_sp_fitch
            nuc_one_hot new_sp_fitch = fitchMerge(src_fitch, dst_fitch);
            nuc_one_hot dst_old = dst_fitch, dst_new = new_sp_fitch;
            propagatePath(dst_path, ptn, freq, dst_old, dst_new, score);

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
                propagatePath(above_lca, ptn, freq, old_prop, new_prop, score);
            }
        }
    }

    // Track O(M) stats
    s_total_M += local_M;
    s_total_P += num_patterns;
    s_move_count++;
    if (local_M > s_max_M) s_max_M = local_M;

    return score;
}

void SPRDeltaExact::printSparseStats() {
    if (s_move_count == 0) {
        std::cout << "  [SPARSE STATS] no moves recorded" << std::endl;
        return;
    }
    std::cout << "  [SPARSE STATS] moves=" << s_move_count
              << " avg_M=" << (s_total_M / s_move_count)
              << " max_M=" << s_max_M
              << " P=" << (s_total_P / s_move_count)
              << " speedup=" << (1.0 * s_total_P / s_total_M) << "x"
              << " avg_ratio=" << (100.0 * s_total_M / s_total_P) << "%"
              << std::endl;
    // Reset for next round
    s_total_M = s_total_P = s_move_count = s_max_M = 0;
}

int SPRDeltaExact::calculateParsimonyDelta(PhyloNode* src, PhyloNode* src_parent,
                                           PhyloNode* dst, PhyloNode* dst_parent,
                                           PhyloNode* lca,
                                           const std::map<PhyloNode*, PhyloNode*>* node_parent) {
    if (!src || !src_parent || !dst || !dst_parent || !lca) {
        return 0;
    }

    int num_neighbors = 0;
    for (auto it = src_parent->neighbors.begin(); it != src_parent->neighbors.end(); ++it) {
        num_neighbors++;
    }
    if (num_neighbors <= 2) {
        return 0;
    }

    // Binary src_parent: use sparse Fitch-propagation approach.
    // O(M) via precomputed Fitch diffs — only iterates over patterns where
    // Fitch sets differ on affected edges. M << P for real data.
    if (num_neighbors == 3 && s_fitch) {
        return sparseBinaryDelta(src, src_parent, dst, dst_parent);
    }

    MutationCountChangeCollection mutations = SPRMutationOps::initMutationChange(src, src_parent);
    int parsimony_score_change = 0;
    MutationCountChangeCollection root_mutations_altered;

    int src_delta = computeSrcSideDelta(src, src_parent, lca,
                                        mutations, root_mutations_altered,
                                        node_parent);
    parsimony_score_change += src_delta;

    MutationCountChangeCollection dst_mutations_altered;
    int dst_delta = computeDstSideDelta(src, dst, dst_parent, lca,
                                        mutations, root_mutations_altered,
                                        node_parent, &dst_mutations_altered);
    parsimony_score_change += dst_delta;

    int lca_delta = computeAboveLCADelta(lca, root_mutations_altered);
    parsimony_score_change += lca_delta;

    // Process dst-side Fitch changes at LCA for positions not already
    // handled by src-side propagation
    for (const auto& dst_change : dst_mutations_altered) {
        bool overlaps = false;
        for (const auto& src_change : root_mutations_altered) {
            if (src_change.position == dst_change.position) {
                overlaps = true;
                break;
            }
        }
        if (!overlaps) {
            int dst_lca_delta = SPRMutationOps::recomputeFullDelta(
                lca, dst_change.position, dst_change.removed_alleles, dst_change.added_alleles,
                dst_change.par_state);
            parsimony_score_change += dst_lca_delta;
        }
    }

    return parsimony_score_change;
}

int SPRDeltaExact::calculateParsimonyDelta(PhyloNode* src, PhyloNode* src_parent,
                                           PhyloNode* dst, PhyloNode* dst_parent,
                                           PhyloTree* tree) {
    PhyloNode* lca = findLCA(src, dst, tree);
    if (!lca) {
        return 0;
    }

    return calculateParsimonyDelta(src, src_parent, dst, dst_parent, lca);
}


int SPRDeltaExact::computeSrcSideDelta(PhyloNode* src, PhyloNode* src_parent,
                                       PhyloNode* lca,
                                       MutationCountChangeCollection& mutations,
                                       MutationCountChangeCollection& root_mutations_altered,
                                       const std::map<PhyloNode*, PhyloNode*>* node_parent) {
    int score_change = 0;

    // Always compute the effect of removing src from src_parent
    SPRMutationOps::getParentAlteredRemove(root_mutations_altered, src, score_change);

    if (src_parent == lca) {
        // No intermediate nodes to walk, but we still needed getParentAlteredRemove above
        return score_change;
    }

    // Only merge here if we can't walk the path (no node_parent map).
    // Otherwise the loop below will handle it starting from src_parent.
    if (!node_parent || node_parent->empty()) {
        SPRMutationOps::mergeMutationSrcToLCA(src_parent, mutations);
    }

    int path_steps = 0;
    if (node_parent && !node_parent->empty()) {
        PhyloNode* current = src_parent;
        PhyloNode* last_src_branch_node = nullptr;
        PhyloNode* child_on_path = src;  // first iteration: changes came from removing src

        while (current != lca) {
            auto parent_it = node_parent->find(current);
            if (parent_it == node_parent->end() || !parent_it->second) {
                break;
            }

            PhyloNode* parent = parent_it->second;
            path_steps++;

            MutationCountChangeCollection new_alter_mutations;
            SPRMutationOps::getIntermediateNodesMutations(
                current, root_mutations_altered,
                new_alter_mutations, score_change,
                child_on_path
            );

            last_src_branch_node = current;
            SPRMutationOps::mergeMutationSrcToLCA(current, mutations);
            root_mutations_altered = std::move(new_alter_mutations);

            child_on_path = current;  // track for next iteration
            current = parent;
        }
    }

    return score_change;
}


int SPRDeltaExact::computeDstSideDelta(PhyloNode* src, PhyloNode* dst,
                                       PhyloNode* dst_parent, PhyloNode* lca,
                                       const MutationCountChangeCollection& mutations,
                                       const MutationCountChangeCollection& root_mutations_altered,
                                       const std::map<PhyloNode*, PhyloNode*>* node_parent,
                                       MutationCountChangeCollection* dst_added) {
    std::vector<PhyloNode*> dst_path;
    PhyloNode* current = dst;

    if (node_parent && !node_parent->empty()) {
        while (current != lca) {
            dst_path.push_back(current);

            auto parent_it = node_parent->find(current);
            if (parent_it == node_parent->end() || !parent_it->second) {
                break;
            }

            current = parent_it->second;
        }
    } else {
        dst_path.push_back(dst);
    }

    MutationCountChangeCollection current_mutations = mutations;

    for (int i = dst_path.size() - 1; i >= 1; --i) {
        current_mutations = SPRMutationOps::mergeMutationLCAToRank(
            dst_path[i], current_mutations);
    }

    if (dst == lca) {
        return checkMoveProfitableLCA(src, lca, current_mutations,
                                      root_mutations_altered, 0, dst_added);
    } else {
        return checkMoveProfitableDstNotLCA(src, dst, lca, current_mutations,
                                            root_mutations_altered, 0, dst_added);
    }
}


int SPRDeltaExact::computeAboveLCADelta(PhyloNode* lca,
                                        const MutationCountChangeCollection& lca_changes) {
    if (lca_changes.empty()) {
        return 0;
    }

    int score_change = 0;
    MutationCountChangeCollection parent_changes;

    SPRMutationOps::checkParsimonyScoreChangeAboveLCA(
        lca, score_change, lca_changes, parent_changes);

    return score_change;
}

int SPRDeltaExact::computeLCAAndAboveDelta(PhyloNode* lca, PhyloNode* src,
                                            const MutationCountChangeCollection& src_changes,
                                            const MutationCountChangeCollection& dst_changes) {
    if (src_changes.empty() && dst_changes.empty()) {
        return 0;
    }

    // Merge src and dst changes using cancel formula.
    MutationCountChangeCollection combined;
    combined.reserve(src_changes.size() + dst_changes.size());

    auto src_it = src_changes.begin(), src_end = src_changes.end();
    auto dst_it = dst_changes.begin(), dst_end = dst_changes.end();

    while (src_it != src_end && dst_it != dst_end) {
        if (src_it->position < dst_it->position) {
            combined.push_back(*src_it);
            src_it++;
        } else if (dst_it->position < src_it->position) {
            combined.push_back(*dst_it);
            dst_it++;
        } else {
            // Same position: merge using cancel formula
            uint8_t any_inc = (src_it->added_alleles & ~dst_it->removed_alleles) |
                             (dst_it->added_alleles & ~src_it->removed_alleles);
            uint8_t any_dec = (src_it->removed_alleles & ~dst_it->added_alleles) |
                             (dst_it->removed_alleles & ~src_it->added_alleles);

            if (any_inc || any_dec) {
                MutationCountChange merged(src_it->position, any_dec, any_inc);
                merged.par_state = src_it->par_state;
                merged.major_allele_set = src_it->major_allele_set;
                merged.boundary1_allele = src_it->boundary1_allele;
                combined.push_back(merged);
            }
            src_it++;
            dst_it++;
        }
    }

    while (src_it != src_end) {
        combined.push_back(*src_it);
        src_it++;
    }

    while (dst_it != dst_end) {
        combined.push_back(*dst_it);
        dst_it++;
    }

    if (combined.empty()) {
        return 0;
    }

    int score_change = 0;
    MutationCountChangeCollection parent_changes;

    SPRMutationOps::checkParsimonyScoreChangeAboveLCA(
        lca, score_change, combined, parent_changes);

    return score_change;
}

int SPRDeltaExact::checkMoveProfitableLCA(PhyloNode* src, PhyloNode* lca,
                                          const MutationCountChangeCollection& mutations,
                                          const MutationCountChangeCollection& root_mutations,
                                          int base_score_change,
                                          MutationCountChangeCollection* dst_added) {
    if (!src || !lca) {
        return base_score_change;
    }

    int score_change = base_score_change;

    PhyloNode* src_branch_node = nullptr;
    PhyloNode* curr = src;
    PhyloNode* curr_parent = SPRMutationOps::getParent(curr);

    while (curr_parent && curr_parent != lca) {
        curr = curr_parent;
        curr_parent = SPRMutationOps::getParent(curr);
    }

    if (curr_parent == lca) {
        src_branch_node = curr;
    }

    if (!src_branch_node) {
        return base_score_change;
    }

    std::vector<Mutation>* src_branch_mutations =
        SPRMutationOps::getMutations(src_branch_node, lca);

    if (!src_branch_mutations) {
        return base_score_change;
    }

    auto mutations_iter = mutations.begin();
    auto mutations_end = mutations.end();
    auto root_iter = root_mutations.begin();
    auto root_end = root_mutations.end();

    for (const Mutation& branch_mut : *src_branch_mutations) {
        if (!branch_mut.is_valid()) {
            continue;
        }

        int pos = branch_mut.position;

        bool src_has_mutation = false;
        uint8_t src_allele = 0;

        while (mutations_iter != mutations_end && mutations_iter->position < pos) {
            mutations_iter++;
        }

        if (mutations_iter != mutations_end && mutations_iter->position == pos) {
            src_has_mutation = true;
            src_allele = mutations_iter->added_alleles;
        } else {
            src_allele = branch_mut.get_par_one_hot();
        }

        uint8_t branch_allele = branch_mut.major_allele_set;

        while (root_iter != root_end && root_iter->position < pos) {
            root_iter++;
        }

        if (root_iter != root_end && root_iter->position == pos) {
            branch_allele = (branch_allele | root_iter->added_alleles) &
                          (~root_iter->removed_alleles);
        }

        uint8_t new_node_allele = src_allele & branch_allele;

        if (!new_node_allele) {
            score_change++;
            new_node_allele = src_allele | branch_allele;
        }

        uint8_t lca_parent_allele = branch_mut.get_par_one_hot();

        bool had_mutation = (branch_allele != lca_parent_allele) &&
                          !(branch_allele & lca_parent_allele);
        bool has_mutation = (new_node_allele != lca_parent_allele) &&
                          !(new_node_allele & lca_parent_allele);

        if (had_mutation && !has_mutation) {
            score_change--;
        } else if (!had_mutation && has_mutation) {
            score_change++;
        }

        if (dst_added && new_node_allele != branch_allele) {
            dst_added->emplace_back(
                pos,
                branch_allele & ~new_node_allele,
                new_node_allele & ~branch_allele
            );
            dst_added->back().par_state = lca_parent_allele;
            dst_added->back().major_allele_set = new_node_allele;
        }

    }

    for (const auto& mut_change : mutations) {
        int pos = mut_change.position;

        bool branch_has_mutation = false;
        for (const Mutation& bm : *src_branch_mutations) {
            if (bm.position == pos) { branch_has_mutation = true; break; }
            if (bm.position > pos) break;
        }

        if (!branch_has_mutation) {
            uint8_t src_allele = mut_change.added_alleles;
            uint8_t parent_allele = mut_change.par_state;

            if (src_allele && !(src_allele & parent_allele)) {
                if (!mut_change.from_src &&
                    (src_allele & mut_change.boundary1_allele)) {
                    // Fitch cascade cancels - net effect is 0
                } else {
                    score_change++;
                }

                if (dst_added) {
                    dst_added->emplace_back(
                        pos,
                        0,
                        src_allele & ~parent_allele
                    );
                    dst_added->back().par_state = parent_allele;
                    dst_added->back().major_allele_set = src_allele | parent_allele;
                }
            }
        }
    }

    return score_change;
}

int SPRDeltaExact::checkMoveProfitableDstNotLCA(PhyloNode* src, PhyloNode* dst,
                                                PhyloNode* lca,
                                                const MutationCountChangeCollection& mutations,
                                                const MutationCountChangeCollection& root_mutations,
                                                int base_score_change,
                                                MutationCountChangeCollection* dst_added) {
    if (!src || !dst || !lca) {
        return base_score_change;
    }

    PhyloNode* dst_parent = SPRMutationOps::getParent(dst);
    if (!dst_parent) {
        return base_score_change;
    }

    int score_change = base_score_change;

    std::vector<Mutation>* dst_mutations = SPRMutationOps::getMutations(dst, dst_parent);

    if (!dst_mutations) {
        return base_score_change;
    }

    auto src_mutations_iter = mutations.begin();
    auto src_mutations_end = mutations.end();

    for (const Mutation& dst_mut : *dst_mutations) {
        if (!dst_mut.is_valid()) {
            continue;
        }

        int pos = dst_mut.position;

        uint8_t src_allele = 0;
        bool src_has_allele = false;

        while (src_mutations_iter != src_mutations_end &&
               src_mutations_iter->position < pos) {
            src_mutations_iter++;
        }

        if (src_mutations_iter != src_mutations_end &&
            src_mutations_iter->position == pos) {
            src_allele = src_mutations_iter->added_alleles;
            src_has_allele = true;
        } else {
            src_allele = dst_mut.get_par_one_hot();
        }

        uint8_t dst_allele = dst_mut.major_allele_set;

        uint8_t new_internal_allele = src_allele & dst_allele;

        int union_score = 0;
        if (!new_internal_allele) {
            score_change++;
            union_score = 1;
            new_internal_allele = src_allele | dst_allele;
        }

        uint8_t parent_allele = dst_mut.get_par_one_hot();

        int new_mutations = 0;
        int old_mutations = 0;

        if (!(new_internal_allele & parent_allele)) {
            new_mutations++;
        }

        if (src_has_allele && !(src_allele & new_internal_allele)) {
            new_mutations++;
        }

        if (!(dst_allele & new_internal_allele)) {
            new_mutations++;
        }

        if (!(dst_allele & parent_allele)) {
            old_mutations++;
        }

        score_change += (new_mutations - old_mutations);

        if (dst_added && new_internal_allele != dst_allele) {
            dst_added->emplace_back(
                pos,
                dst_allele & ~new_internal_allele,
                new_internal_allele & ~dst_allele
            );
            dst_added->back().par_state = parent_allele;
            dst_added->back().major_allele_set = new_internal_allele;
        }
    }

    src_mutations_iter = mutations.begin();
    for (; src_mutations_iter != src_mutations_end; ++src_mutations_iter) {
        int pos = src_mutations_iter->position;

        bool dst_has_mutation = false;
        for (const Mutation& dst_mut : *dst_mutations) {
            if (dst_mut.position == pos) {
                dst_has_mutation = true;
                break;
            }
            if (dst_mut.position > pos) {
                break;
            }
        }

        if (!dst_has_mutation) {
            uint8_t src_allele = src_mutations_iter->added_alleles;
            uint8_t parent_allele = src_mutations_iter->par_state;

            if (!(src_allele & parent_allele)) {
                bool cascade_cancel = false;
                if (!src_mutations_iter->from_src &&
                    (src_allele & src_mutations_iter->boundary1_allele)) {
                    cascade_cancel = true;
                } else {
                    score_change++;
                }
                if (dst_added) {
                    dst_added->emplace_back(
                        pos,
                        0,
                        src_allele & ~parent_allele
                    );
                    dst_added->back().par_state = parent_allele;
                    dst_added->back().major_allele_set = src_allele | parent_allele;
                }
            }

        }
    }

    return score_change;
}

PhyloNode* SPRDeltaExact::findLCA(PhyloNode* node1, PhyloNode* node2, PhyloTree* tree) {
    if (!node1 || !node2) return nullptr;
    if (node1 == node2) return node1;

    if (s_lca_ready &&
        node1->id < (int)s_first_appearance.size() && s_first_appearance[node1->id] >= 0 &&
        node2->id < (int)s_first_appearance.size() && s_first_appearance[node2->id] >= 0) {
        int l = s_first_appearance[node1->id];
        int r = s_first_appearance[node2->id];
        int idx = rmq(l, r);
        return s_euler[idx];
    }

    int d1, d2;
    if (!s_node_depth.empty()) {
        d1 = (node1->id < (int)s_node_depth.size() && s_node_depth[node1->id] >= 0)
             ? s_node_depth[node1->id] : 0;
        d2 = (node2->id < (int)s_node_depth.size() && s_node_depth[node2->id] >= 0)
             ? s_node_depth[node2->id] : 0;
    } else {
        d1 = d2 = 0;
        for (PhyloNode* c = node1; c; c = SPRMutationOps::getParent(c)) d1++;
        for (PhyloNode* c = node2; c; c = SPRMutationOps::getParent(c)) d2++;
    }

    PhyloNode* a = node1;
    PhyloNode* b = node2;
    while (d1 > d2) { a = SPRMutationOps::getParent(a); d1--; }
    while (d2 > d1) { b = SPRMutationOps::getParent(b); d2--; }
    while (a != b) {
        a = SPRMutationOps::getParent(a);
        b = SPRMutationOps::getParent(b);
    }
    return a;
}
