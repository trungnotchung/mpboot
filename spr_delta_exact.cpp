#include "spr_delta_exact.h"
#include "spr_mutation_ops.h"
#include "fitch.h"
#include <algorithm>
#include <map>
#include <queue>
#include <set>

// Static state for binary fallback
static Fitch* s_cf = nullptr;
static int s_current_score = 0;

// Precomputed depth cache for O(1) depth lookup in findLCA.
// Populated by precomputeDepths(), keyed by node pointer.
static std::map<PhyloNode*, int> s_node_depth;

void SPRDeltaExact::setCustomFitch(Fitch* cf, int current_score) {
    s_cf = cf;
    s_current_score = current_score;
}

void SPRDeltaExact::precomputeDepths(PhyloTree* tree) {
    s_node_depth.clear();
    if (!tree || !tree->root) return;
    // BFS from root, assigning depth
    std::queue<std::pair<PhyloNode*, int>> q;
    q.push({(PhyloNode*)tree->root, 0});
    std::set<PhyloNode*> visited;
    while (!q.empty()) {
        PhyloNode* node = q.front().first;
        int depth = q.front().second;
        q.pop();
        if (visited.count(node)) continue;
        visited.insert(node);
        s_node_depth[node] = depth;
        FOR_NEIGHBOR_IT(node, nullptr, it) {
            PhyloNode* nb = (PhyloNode*)(*it)->node;
            if (!visited.count(nb)) q.push({nb, depth + 1});
        }
    }
}

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

struct SparsePathNode {
    PhyloNode* node;
    const nuc_one_hot* other_arr;   // sibling (unchanged) Fitch array
    const nuc_one_hot* node_arr;    // this node's current Fitch array
};

static std::vector<SparsePathNode> buildPath(PhyloNode* prev, PhyloNode* cur,
                                               PhyloNode* stop_node, Fitch* cf) {
    std::vector<SparsePathNode> path;
    while (cur && cur != stop_node) {
        PhyloNode* par = SPRMutationOps::getParent(cur);
        PhyloNode* other = findOtherChild(cur, par, prev);
        const nuc_one_hot* other_arr = other
            ? cf->getMajorArrayForNode(other)
            : cf->getMajorArrayForNode(cur);  // root leaf: use own state
        path.push_back({cur, other_arr, cf->getMajorArrayForNode(cur)});
        prev = cur;
        cur = par;
    }
    return path;
}

static void propagatePath(const std::vector<SparsePathNode>& path, int p, int freq,
                           nuc_one_hot& old_f, nuc_one_hot& new_f, int& score) {
    for (size_t i = 0; i < path.size(); i++) {
        if (old_f == new_f) break;
        const SparsePathNode& pn = path[i];
        nuc_one_hot other_f = pn.other_arr[p];
        nuc_one_hot new_node_f;
        score += penaltyDelta(old_f, new_f, other_f, new_node_f) * freq;
        old_f = pn.node_arr[p];
        new_f = new_node_f;
    }
}

static inline void addDiffs(Fitch* cf, PhyloNode* node, std::vector<int>& out) {
    const std::vector<int>* d = cf->getFitchDiffs(node);
    if (d && !d->empty()) out.insert(out.end(), d->begin(), d->end());
}

static inline void finalizeDiffs(std::vector<int>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
}

static long long s_total_M = 0, s_total_P = 0;
static int s_move_count = 0;
static int s_max_M = 0;

static int sparseBinaryDelta(PhyloNode* src, PhyloNode* src_parent,
                              PhyloNode* dst, PhyloNode* dst_parent) {
    if (!s_cf) return 0;
    Fitch* cf = s_cf;

    PhyloNode* gp  = SPRMutationOps::getParent(src_parent);
    PhyloNode* sib = findOtherChild(src_parent, gp, src);
    if (!sib || !gp) return 0;

    PhyloNode* lca = SPRDeltaExact::findLCA(src_parent, dst, nullptr);
    if (!lca) return 0;

    const nuc_one_hot* src_arr = cf->getMajorArrayForNode(src);
    const nuc_one_hot* sib_arr = cf->getMajorArrayForNode(sib);
    const nuc_one_hot* sp_arr  = cf->getMajorArrayForNode(src_parent);
    const nuc_one_hot* dst_arr = cf->getMajorArrayForNode(dst);

    int score = 0;
    int nptn = cf->getNumPatterns();
    int local_M = 0; // for stats

    if (lca == src_parent) {
        // ========== CASE A: dst in sib's subtree ==========
        // sp bypassed (sib→gp), sp reinserted at dp with children {src, dst}.
        // Chain: dp → ... → sib → gp → above.

        std::vector<SparsePathNode> dp_to_sib = buildPath(dst, dst_parent, src_parent, cf);

        PhyloNode* gp_other = findOtherChild(gp, SPRMutationOps::getParent(gp), src_parent);
        const nuc_one_hot* gp_other_arr = gp_other ? cf->getMajorArrayForNode(gp_other) : nullptr;
        const nuc_one_hot* gp_arr = cf->getMajorArrayForNode(gp);

        std::vector<SparsePathNode> above_gp = buildPath(
            gp, SPRMutationOps::getParent(gp), nullptr, cf);

        // Collect O(M) affected patterns
        std::vector<int> affected;
        addDiffs(cf, src, affected);
        addDiffs(cf, sib, affected);
        addDiffs(cf, src_parent, affected);
        addDiffs(cf, dst, affected);
        for (const auto& pn : dp_to_sib) addDiffs(cf, pn.node, affected);
        finalizeDiffs(affected);
        local_M = (int)affected.size();

        for (int p : affected) {
            int freq = cf->getPatternFreq(p);
            nuc_one_hot src_f = src_arr[p], sib_f = sib_arr[p];
            nuc_one_hot dst_f = dst_arr[p],  sp_f = sp_arr[p];

            // 1. sp penalty: old !(src∩sib) → new !(src∩dst)
            score += (((src_f & dst_f) ? 0 : 1) - ((src_f & sib_f) ? 0 : 1)) * freq;

            // 2. Propagate dp→sib: child changed from dst_f to sp_new_f
            nuc_one_hot old_f = dst_f;
            nuc_one_hot new_f = fitchMerge(src_f, dst_f);
            propagatePath(dp_to_sib, p, freq, old_f, new_f, score);

            // After propagation, new_f = sib's new Fitch (if changed)
            nuc_one_hot new_sib_f = (!dp_to_sib.empty() && old_f != new_f) ? new_f : sib_f;

            // 3. At gp: child changed from sp_f to new_sib_f
            if (gp->isLeaf()) {
                // Root edge: root—sp → root—sib
                score += rootEdgeDelta(gp_arr[p], sp_f, new_sib_f) * freq;
            } else {
                nuc_one_hot gp_other_f = gp_other_arr ? gp_other_arr[p] : 0xF;
                nuc_one_hot new_gp_f;
                score += penaltyDelta(sp_f, new_sib_f, gp_other_f, new_gp_f) * freq;

                if (new_gp_f != gp_arr[p]) {
                    nuc_one_hot ao = gp_arr[p], an = new_gp_f;
                    propagatePath(above_gp, p, freq, ao, an, score);
                }
            }
        }

    } else if (dst == lca) {
        // ========== CASE B1: sp_new inserted on lca—src_branch edge ==========
        // src-side: sp_f→sib_f propagates gp→...→src_branch_at_lca
        // sp_new sits between lca and src_branch, children = {src, updated_src_branch}
        // At lca: child changed from src_branch to sp_new

        std::vector<SparsePathNode> src_path = buildPath(src_parent, gp, lca, cf);
        PhyloNode* src_branch = src_path.empty() ? src_parent : src_path.back().node;
        const nuc_one_hot* src_branch_arr = cf->getMajorArrayForNode(src_branch);

        PhyloNode* lca_parent    = SPRMutationOps::getParent(lca);
        PhyloNode* lca_other     = findOtherChild(lca, lca_parent, src_branch);
        const nuc_one_hot* lca_other_arr = lca_other ? cf->getMajorArrayForNode(lca_other) : nullptr;
        const nuc_one_hot* lca_arr = cf->getMajorArrayForNode(lca);

        // parent(lca) info — for propagation above lca
        PhyloNode* plca_other = lca_parent ? findOtherChild(lca_parent,
            SPRMutationOps::getParent(lca_parent), lca) : nullptr;
        const nuc_one_hot* plca_other_arr = plca_other ? cf->getMajorArrayForNode(plca_other) : nullptr;
        const nuc_one_hot* plca_arr = lca_parent ? cf->getMajorArrayForNode(lca_parent) : nullptr;

        std::vector<SparsePathNode> above_lca_parent = buildPath(
            lca_parent, lca_parent ? SPRMutationOps::getParent(lca_parent) : nullptr,
            nullptr, cf);

        // Collect O(M) affected patterns
        std::vector<int> affected;
        addDiffs(cf, src, affected);
        addDiffs(cf, sib, affected);
        addDiffs(cf, src_parent, affected);
        for (const auto& pn : src_path) addDiffs(cf, pn.node, affected);
        addDiffs(cf, lca, affected);
        finalizeDiffs(affected);
        local_M = (int)affected.size();

        for (int p : affected) {
            int freq = cf->getPatternFreq(p);
            nuc_one_hot src_f = src_arr[p], sib_f = sib_arr[p], sp_f = sp_arr[p];
            int old_sp_pen = (src_f & sib_f) ? 0 : 1;

            // 1. src-side propagation: sp_f → sib_f through gp→...→src_branch
            nuc_one_hot src_old = sp_f, src_new = sib_f;
            propagatePath(src_path, p, freq, src_old, src_new, score);

            nuc_one_hot src_child_old = src_branch_arr[p];
            nuc_one_hot src_child_new = (src_old != src_new) ? src_new : src_child_old;

            // 2. sp_new penalty: !(src ∩ updated_src_branch)
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
                // lca IS root leaf: root edge = root—src_branch → root—sp_new
                score += rootEdgeDelta(lca_arr[p], src_child_old, sp_new_f) * freq;
            } else if (new_lca_f != old_lca_f && lca_parent) {
                if (lca_parent->isLeaf()) {
                    // parent(lca) is root leaf
                    score += rootEdgeDelta(plca_arr[p], old_lca_f, new_lca_f) * freq;
                } else {
                    nuc_one_hot plca_other_f = plca_other_arr ? plca_other_arr[p] : 0xF;
                    nuc_one_hot new_plca_f;
                    score += penaltyDelta(old_lca_f, new_lca_f, plca_other_f, new_plca_f) * freq;

                    if (new_plca_f != (plca_arr ? plca_arr[p] : (nuc_one_hot)0xF)) {
                        nuc_one_hot ao = plca_arr[p], an = new_plca_f;
                        propagatePath(above_lca_parent, p, freq, ao, an, score);
                    }
                }
            }
        }

    } else {
        // ========== CASE B2: general (lca above sp, dst below lca) ==========
        // Independent src-side and dst-side paths meet at lca.

        std::vector<SparsePathNode> src_path = buildPath(src_parent, gp, lca, cf);
        std::vector<SparsePathNode> dst_path = buildPath(dst, dst_parent, lca, cf);

        PhyloNode* src_branch = src_path.empty() ? src_parent : src_path.back().node;
        PhyloNode* dst_branch = dst_path.empty() ? dst : dst_path.back().node;
        const nuc_one_hot* src_branch_arr = cf->getMajorArrayForNode(src_branch);
        const nuc_one_hot* dst_branch_arr = cf->getMajorArrayForNode(dst_branch);

        const nuc_one_hot* lca_arr = cf->getMajorArrayForNode(lca);
        std::vector<SparsePathNode> above_lca = buildPath(
            lca, SPRMutationOps::getParent(lca), nullptr, cf);

        // Collect O(M) affected patterns
        std::vector<int> affected;
        addDiffs(cf, src, affected);
        addDiffs(cf, sib, affected);
        addDiffs(cf, src_parent, affected);
        addDiffs(cf, dst, affected);
        for (const auto& pn : src_path) addDiffs(cf, pn.node, affected);
        for (const auto& pn : dst_path) addDiffs(cf, pn.node, affected);
        finalizeDiffs(affected);
        local_M = (int)affected.size();

        for (int p : affected) {
            int freq = cf->getPatternFreq(p);
            nuc_one_hot src_f = src_arr[p], sib_f = sib_arr[p];
            nuc_one_hot dst_f = dst_arr[p],  sp_f = sp_arr[p];

            // 1. sp penalty: old !(src∩sib) → new !(src∩dst)
            score += (((src_f & dst_f) ? 0 : 1) - ((src_f & sib_f) ? 0 : 1)) * freq;

            // 2. src-side: sp_f → sib_f
            nuc_one_hot src_old = sp_f, src_new = sib_f;
            propagatePath(src_path, p, freq, src_old, src_new, score);

            // 3. dst-side: dst_f → sp_new_f
            nuc_one_hot sp_new_f = fitchMerge(src_f, dst_f);
            nuc_one_hot dst_old = dst_f, dst_new = sp_new_f;
            propagatePath(dst_path, p, freq, dst_old, dst_new, score);

            // 4. At LCA: combine changes from both branches
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
                propagatePath(above_lca, p, freq, ao, an, score);
            }
        }
    }

    // Track O(M) stats
    s_total_M += local_M;
    s_total_P += nptn;
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
    if (num_neighbors == 3 && s_cf) {
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
                lca, dst_change.position, dst_change.decremented, dst_change.incremented,
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
    // Key insight from dbl_inc_dec_mutations: src branch decrements alleles
    // while dst branch increments them (opposite directions), so changes
    // cancel at the same allele. Net change per allele is at most ±1.
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
            // Same position: merge using cancel formula (dbl_inc_dec_mutations)
            uint8_t any_inc = (src_it->incremented & ~dst_it->decremented) |
                             (dst_it->incremented & ~src_it->decremented);
            uint8_t any_dec = (src_it->decremented & ~dst_it->incremented) |
                             (dst_it->decremented & ~src_it->incremented);

            if (any_inc || any_dec) {
                MutationCountChange merged(src_it->position, any_dec, any_inc);
                merged.par_state = src_it->par_state;  // Both should have same par_state (LCA's parent allele)
                merged.major_allele = src_it->major_allele;
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
            src_allele = mutations_iter->incremented;
        } else {
            // When src has no mutation at this position in the collection,
            // src's allele defaults to the LCA's allele (same as branch parent).
            // This mirrors checkMoveProfitableDstNotLCA's default behavior.
            src_allele = branch_mut.get_par_one_hot();
        }

        uint8_t branch_allele = branch_mut.all_major_allele;

        while (root_iter != root_end && root_iter->position < pos) {
            root_iter++;
        }

        if (root_iter != root_end && root_iter->position == pos) {
            branch_allele = (branch_allele | root_iter->incremented) &
                          (~root_iter->decremented);
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

        // Track Fitch set changes at the mezzanine level for LCA propagation
        if (dst_added && new_node_allele != branch_allele) {
            dst_added->emplace_back(
                pos,
                branch_allele & ~new_node_allele,    // decremented
                new_node_allele & ~branch_allele      // incremented
            );
            dst_added->back().par_state = lca_parent_allele;
            dst_added->back().major_allele = new_node_allele;
        }

    }

    // Process positions in mutations NOT in src_branch_mutations.
    // These are positions where src introduces a new allele at LCA that the
    // old src_branch didn't have. The src-side counted removing these mutations,
    // so we must account for re-inserting them at the destination.
    // (Analogous to the second loop in checkMoveProfitableDstNotLCA.)
    for (const auto& mut_change : mutations) {
        int pos = mut_change.position;

        bool branch_has_mutation = false;
        for (const Mutation& bm : *src_branch_mutations) {
            if (bm.position == pos) { branch_has_mutation = true; break; }
            if (bm.position > pos) break;
        }

        if (!branch_has_mutation) {
            uint8_t src_allele = mut_change.incremented;
            uint8_t parent_allele = mut_change.par_state;

            if (src_allele && !(src_allele & parent_allele)) {
                // For path-edge entries: if src allele is in boundary1,
                // Fitch cascade at intermediate node cancels this mutation
                if (!mut_change.from_src &&
                    (src_allele & mut_change.boundary1_allele)) {
                    // Fitch cascade cancels - net effect is 0
                } else {
                    score_change++;
                }

                // Track Fitch change at mezzanine level
                if (dst_added) {
                    dst_added->emplace_back(
                        pos,
                        0,
                        src_allele & ~parent_allele
                    );
                    dst_added->back().par_state = parent_allele;
                    dst_added->back().major_allele = src_allele | parent_allele;
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
            src_allele = src_mutations_iter->incremented;
            src_has_allele = true;
        } else {
            src_allele = dst_mut.get_par_one_hot();
        }

        uint8_t dst_allele = dst_mut.all_major_allele;

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

        // Track Fitch set changes for LCA propagation
        if (dst_added && new_internal_allele != dst_allele) {
            dst_added->emplace_back(
                pos,
                dst_allele & ~new_internal_allele,    // decremented
                new_internal_allele & ~dst_allele      // incremented
            );
            dst_added->back().par_state = parent_allele;
            dst_added->back().major_allele = new_internal_allele;
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
            uint8_t src_allele = src_mutations_iter->incremented;
            uint8_t parent_allele = src_mutations_iter->par_state;

            if (!(src_allele & parent_allele)) {
                bool cascade_cancel = false;
                if (!src_mutations_iter->from_src &&
                    (src_allele & src_mutations_iter->boundary1_allele)) {
                    cascade_cancel = true;
                } else {
                    score_change++;
                }
                // Track Fitch change: new_node has {src_allele, parent_allele}
                // Old child (dst) had {parent_allele}. Added: src_allele.
                if (dst_added) {
                    dst_added->emplace_back(
                        pos,
                        0,                              // nothing decremented
                        src_allele & ~parent_allele     // src_allele added
                    );
                    dst_added->back().par_state = parent_allele;
                    dst_added->back().major_allele = src_allele | parent_allele;
                }
            }

        }
    }

    return score_change;
}

// Depth-walk LCA: O(depth_diff + lca_depth) time, O(1) space — no heap allocation.
// Uses precomputed depths for O(1) lookup when available, else walks to root.
PhyloNode* SPRDeltaExact::findLCA(PhyloNode* node1, PhyloNode* node2, PhyloTree* tree) {
    if (!node1 || !node2) return nullptr;
    if (node1 == node2) return node1;

    int d1, d2;
    if (!s_node_depth.empty()) {
        // Use precomputed depths — O(1) lookup
        auto it1 = s_node_depth.find(node1);
        auto it2 = s_node_depth.find(node2);
        d1 = (it1 != s_node_depth.end()) ? it1->second : 0;
        d2 = (it2 != s_node_depth.end()) ? it2->second : 0;
    } else {
        // Walk to root to find depths (getParent is O(1) each)
        d1 = d2 = 0;
        for (PhyloNode* c = node1; c; c = SPRMutationOps::getParent(c)) d1++;
        for (PhyloNode* c = node2; c; c = SPRMutationOps::getParent(c)) d2++;
    }

    // Bring deeper node up to same depth
    PhyloNode* a = node1;
    PhyloNode* b = node2;
    while (d1 > d2) { a = SPRMutationOps::getParent(a); d1--; }
    while (d2 > d1) { b = SPRMutationOps::getParent(b); d2--; }

    // Walk both up until they meet
    while (a != b) {
        a = SPRMutationOps::getParent(a);
        b = SPRMutationOps::getParent(b);
    }
    return a;
}
