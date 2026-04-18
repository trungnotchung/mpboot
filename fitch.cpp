#include "fitch.h"
#include "phylotree.h"
#include "phylonode.h"
#include "onehot_encoding.h"
#include <iostream>
#include <queue>
#include <set>
#include <cassert>
#include <algorithm>

using namespace std;

extern unsigned int dna_state_map[128];

Fitch::Fitch(PhyloTree* tree) : tree(tree), nptn(0), num_nodes(0), root_side_mutation_count(0) {
    assert(tree != nullptr);
    assert(tree->aln != nullptr);
}

inline int Fitch::getIdx(PhyloNode* node) const {
    auto it = node_index.find(node);
    assert(it != node_index.end());
    return it->second;
}


nuc_one_hot Fitch::chooseRepresentative(nuc_one_hot state_set) {
    if (state_set & 0x1) return 0x1;  // A
    if (state_set & 0x2) return 0x2;  // C
    if (state_set & 0x4) return 0x4;  // G
    if (state_set & 0x8) return 0x8;  // T
    return 0xF;
}

void Fitch::buildPatternPositionMap() {
    std::fill(ptn_position.begin(), ptn_position.end(), -1);
    std::set<PhyloNode*> visited;
    std::queue<PhyloNode*> q;
    q.push((PhyloNode*)tree->root);
    while (!q.empty()) {
        PhyloNode* n = q.front(); q.pop();
        if (visited.count(n)) continue;
        visited.insert(n);
        for (auto it = n->neighbors.begin(); it != n->neighbors.end(); ++it) {
            PhyloNeighbor* nei = (PhyloNeighbor*)*it;
            PhyloNode* other = (PhyloNode*)nei->node;
            for (const Mutation& m : nei->mutations) {
                if (m.is_valid() && m.compressed_position >= 0 && m.compressed_position < nptn) {
                    ptn_position[m.compressed_position] = m.position;
                }
            }
            if (!visited.count(other)) q.push(other);
        }
    }
}

void Fitch::buildPatternToSites() {
    int nsites = tree->aln->getNSite();
    pattern_to_sites.resize(nptn);
    for (int site = 0; site < nsites; site++) {
        int ptn = tree->aln->getPatternID(site);
        pattern_to_sites[ptn].push_back(site);
    }
}

static void fillFitchDiffsRecursive(PhyloNode* node, PhyloNode* parent,
                                     Fitch* cf,
                                     std::vector<std::vector<int>>& fitch_diffs) {
    auto it_idx = cf->getNodeIndex().find(node);
    auto it_par = cf->getNodeIndex().find(parent);
    if (it_idx == cf->getNodeIndex().end() || it_par == cf->getNodeIndex().end()) return;

    int node_idx = it_idx->second;
    int par_idx = it_par->second;
    int nptn = cf->getNumPatterns();
    const nuc_one_hot* node_arr = cf->getMajorArrayForNode(node);
    const nuc_one_hot* par_arr = cf->getMajorArrayForNode(parent);

    auto& diffs = fitch_diffs[node_idx];
    diffs.clear();
    for (int p = 0; p < nptn; p++) {
        if (node_arr[p] != par_arr[p]) diffs.push_back(p);
    }
    FOR_NEIGHBOR_IT(node, parent, nit) {
        fillFitchDiffsRecursive((PhyloNode*)(*nit)->node, node, cf, fitch_diffs);
    }
}

void Fitch::computeFitchDiffs() {
    fitch_diffs.resize(num_nodes);

    PhyloNode* root = (PhyloNode*)tree->root;
    PhyloNode* root_neighbor = (PhyloNode*)root->neighbors[0]->node;

    fitch_diffs[getIdx(root)].clear();

    {
        int rn_idx = getIdx(root_neighbor);
        int r_idx = getIdx(root);
        const nuc_one_hot* rn_arr = majorAt(rn_idx);
        const nuc_one_hot* r_arr = majorAt(r_idx);
        auto& diffs = fitch_diffs[rn_idx];
        diffs.clear();
        for (int p = 0; p < nptn; p++) {
            if (rn_arr[p] != r_arr[p]) diffs.push_back(p);
        }
    }

    FOR_NEIGHBOR_IT(root_neighbor, root, it) {
        fillFitchDiffsRecursive((PhyloNode*)(*it)->node, root_neighbor, this, fitch_diffs);
    }
}

int Fitch::bottomUp(PhyloNode* node, PhyloNode* parent) {
    int idx = getIdx(node);
    nuc_one_hot* my_major = majorAt(idx);

    if (node->isLeaf()) {
        for (int ptn = 0; ptn < nptn; ptn++) {
            int state = (tree->aln->at(ptn))[node->id];
            my_major[ptn] = (nuc_one_hot)(dna_state_map[state] & 0xF);
        }
        node_penalty[idx] = 0;
        subtree_score[idx] = 0;
        return 0;
    }

    PhyloNode* child0 = nullptr;
    PhyloNode* child1 = nullptr;
    FOR_NEIGHBOR_IT(node, parent, it) {
        if (!child0) child0 = (PhyloNode*)(*it)->node;
        else child1 = (PhyloNode*)(*it)->node;
    }
    assert(child0 && child1);

    int children_score = bottomUp(child0, node) + bottomUp(child1, node);

    const nuc_one_hot* c0 = majorAt(getIdx(child0));
    const nuc_one_hot* c1 = majorAt(getIdx(child1));

    int my_penalty = 0;
    for (int ptn = 0; ptn < nptn; ptn++) {
        nuc_one_hot intersect = c0[ptn] & c1[ptn];
        if (intersect != 0) {
            my_major[ptn] = intersect;
        } else {
            my_major[ptn] = c0[ptn] | c1[ptn];
            my_penalty += tree->aln->at(ptn).frequency;
        }
    }

    node_penalty[idx] = my_penalty;
    subtree_score[idx] = children_score + my_penalty;
    return children_score + my_penalty;
}

void Fitch::topDown(PhyloNode* node, PhyloNode* parent,
                          const vector<nuc_one_hot>& parent_states) {
    int idx = getIdx(node);
    const nuc_one_hot* my_major = majorAt(idx);

    vector<nuc_one_hot> my_states(nptn);
    for (int ptn = 0; ptn < nptn; ptn++) {
        nuc_one_hot major = my_major[ptn];
        if (!parent_states.empty() && (parent_states[ptn] & major)) {
            my_states[ptn] = parent_states[ptn];
        } else {
            my_states[ptn] = chooseRepresentative(major);
        }
    }

    const nuc_one_hot* ch0_major = nullptr;
    const nuc_one_hot* ch1_major = nullptr;
    if (!node->isLeaf()) {
        PhyloNode* ch0 = nullptr;
        PhyloNode* ch1 = nullptr;
        FOR_NEIGHBOR_IT(node, parent, child_it) {
            if (!ch0) ch0 = (PhyloNode*)(*child_it)->node;
            else ch1 = (PhyloNode*)(*child_it)->node;
        }
        if (ch0) ch0_major = majorAt(getIdx(ch0));
        if (ch1) ch1_major = majorAt(getIdx(ch1));
    }

    if (parent != nullptr) {
        PhyloNeighbor* edge_to_node = (PhyloNeighbor*)parent->findNeighbor(node);
        PhyloNeighbor* edge_to_parent = (PhyloNeighbor*)node->findNeighbor(parent);

        if (edge_to_node) edge_to_node->mutations.clear();
        if (edge_to_parent) edge_to_parent->mutations.clear();

        for (int ptn = 0; ptn < nptn; ptn++) {
            if (my_states[ptn] == parent_states[ptn]) continue;
            if (ptn_is_const[ptn]) continue;

            // boundary1 = alleles with count = max_count - 1
            nuc_one_hot boundary1 = 0;
            if (ch0_major && ch1_major) {
                nuc_one_hot inter = ch0_major[ptn] & ch1_major[ptn];
                if (inter != 0) {
                    boundary1 = (ch0_major[ptn] | ch1_major[ptn]) & ~inter;
                } else {
                    boundary1 = (~(ch0_major[ptn] | ch1_major[ptn])) & 0xF;
                }
            }

            for (int site : pattern_to_sites[ptn]) {
                Mutation mut;
                mut.position = site;
                mut.compressed_position = ptn;
                mut.par_one_hot = parent_states[ptn];
                mut.mut_one_hot = my_states[ptn];
                mut.par_nuc = OneHotEncoding::oneHotToChar(parent_states[ptn]);
                mut.mut_nuc = OneHotEncoding::oneHotToChar(my_states[ptn]);
                mut.ref_nuc = mut.par_nuc;
                mut.is_missing = false;
                mut.all_major_allele = my_major[ptn];
                mut.boundary1_allele = boundary1;

                if (edge_to_node) edge_to_node->mutations.push_back(mut);
                if (edge_to_parent) edge_to_parent->mutations.push_back(mut);
            }
        }

        if (edge_to_node)
            sort(edge_to_node->mutations.begin(), edge_to_node->mutations.end());
        if (edge_to_parent)
            sort(edge_to_parent->mutations.begin(), edge_to_parent->mutations.end());
    }

    FOR_NEIGHBOR_IT(node, parent, it) {
        PhyloNode* child = (PhyloNode*)(*it)->node;
        topDown(child, node, my_states);
    }
}

int Fitch::run() {
    assert(tree->root != nullptr);
    assert(tree->root->isLeaf());  // unrooted tree convention: root is a leaf

    nptn = tree->aln->size();

    {
        int idx = 0;
        queue<PhyloNode*> q;
        set<PhyloNode*> visited;
        q.push((PhyloNode*)tree->root);
        while (!q.empty()) {
            PhyloNode* n = q.front(); q.pop();
            if (visited.count(n)) continue;
            visited.insert(n);
            node_index[n] = idx++;
            FOR_NEIGHBOR_IT(n, nullptr, it)
                if (!visited.count((PhyloNode*)(*it)->node))
                    q.push((PhyloNode*)(*it)->node);
        }
        num_nodes = idx;
    }

    node_major.assign((size_t)num_nodes * nptn, 0);
    node_penalty.assign(num_nodes, 0);
    subtree_score.assign(num_nodes, 0);

    ptn_freq.resize(nptn);
    ptn_is_const.resize(nptn);
    ptn_position.resize(nptn);
    for (int ptn = 0; ptn < nptn; ptn++) {
        ptn_freq[ptn] = tree->aln->at(ptn).frequency;
        ptn_is_const[ptn] = tree->aln->at(ptn).is_const;
    }

    buildPatternToSites();

    PhyloNode* root = (PhyloNode*)tree->root;
    PhyloNode* root_neighbor = (PhyloNode*)root->neighbors[0]->node;

    int root_idx = getIdx(root);
    nuc_one_hot* root_major = majorAt(root_idx);
    for (int ptn = 0; ptn < nptn; ptn++) {
        int state = (tree->aln->at(ptn))[root->id];
        root_major[ptn] = (nuc_one_hot)(dna_state_map[state] & 0xF);
    }

    int subtree_score = bottomUp(root_neighbor, root);

    int root_edge_score = 0;
    int rn_idx = getIdx(root_neighbor);
    const nuc_one_hot* rn_major = majorAt(rn_idx);
    for (int ptn = 0; ptn < nptn; ptn++) {
        if ((root_major[ptn] & rn_major[ptn]) == 0) {
            root_edge_score += tree->aln->at(ptn).frequency;
        }
    }

    int total_score = subtree_score + root_edge_score;

    // Virtual root state: merge root leaf and root_neighbor at the root edge.
    // The tree is unrooted, so we need a virtual root to determine top-down states.
    // This prevents spurious mutations when root leaf has ambiguous states (N=0xF).

    vector<nuc_one_hot> virtual_root_states(nptn);
    for (int ptn = 0; ptn < nptn; ptn++) {
        nuc_one_hot intersect = root_major[ptn] & rn_major[ptn];
        if (intersect != 0) {
            virtual_root_states[ptn] = chooseRepresentative(intersect);
        } else {
            virtual_root_states[ptn] = chooseRepresentative(root_major[ptn] | rn_major[ptn]);
        }
    }

    topDown(root_neighbor, root, virtual_root_states);

    // Root-side mutations: counted separately (not on root edge) so SPRDeltaExact
    // sees clean mutations on root_neighbor's parent edge.
    root_side_mutation_count = 0;
    {
        for (int ptn = 0; ptn < nptn; ptn++) {
            nuc_one_hot vr = virtual_root_states[ptn];

            nuc_one_hot root_assigned;
            if (vr & root_major[ptn]) {
                root_assigned = vr;
            } else {
                root_assigned = chooseRepresentative(root_major[ptn]);
            }

            if (root_assigned == vr) continue;
            if (ptn_is_const[ptn]) continue;

            root_side_mutation_count += (int)pattern_to_sites[ptn].size();
        }
    }

    computeFitchDiffs();

    cout << "Fitch: score=" << total_score << " patterns=" << nptn
         << " nodes=" << num_nodes << endl;

    return total_score;
}

int Fitch::countMutations() const {
    // BFS counting one direction per edge to avoid double-counting.
    int total = 0;
    set<PhyloNode*> visited;
    queue<PhyloNode*> q;
    q.push((PhyloNode*)tree->root);

    while (!q.empty()) {
        PhyloNode* node = q.front(); q.pop();
        if (visited.count(node)) continue;
        visited.insert(node);

        FOR_NEIGHBOR_IT(node, nullptr, it) {
            PhyloNode* neighbor = (PhyloNode*)(*it)->node;
            if (!visited.count(neighbor)) {
                PhyloNeighbor* edge = (PhyloNeighbor*)(*it);
                total += (int)edge->mutations.size();
                q.push(neighbor);
            }
        }
    }

    total += root_side_mutation_count;

    return total;
}

int Fitch::localBottomUp(PhyloNode* node, PhyloNode* parent) {
    int idx = getIdx(node);
    nuc_one_hot* my_major = majorAt(idx);

    if (node->isLeaf()) {
        subtree_score[idx] = 0;
        node_penalty[idx] = 0;
        return 0;
    }

    PhyloNode* child0 = nullptr;
    PhyloNode* child1 = nullptr;
    FOR_NEIGHBOR_IT(node, parent, it) {
        if (!child0) child0 = (PhyloNode*)(*it)->node;
        else child1 = (PhyloNode*)(*it)->node;
    }
    assert(child0 && child1);

    int children_score = localBottomUp(child0, node) + localBottomUp(child1, node);

    const nuc_one_hot* c0 = majorAt(getIdx(child0));
    const nuc_one_hot* c1 = majorAt(getIdx(child1));

    int my_penalty = 0;
    for (int ptn = 0; ptn < nptn; ptn++) {
        nuc_one_hot intersect = c0[ptn] & c1[ptn];
        if (intersect != 0) {
            my_major[ptn] = intersect;
        } else {
            my_major[ptn] = c0[ptn] | c1[ptn];
            my_penalty += ptn_freq[ptn];
        }
    }

    node_penalty[idx] = my_penalty;
    subtree_score[idx] = children_score + my_penalty;
    return children_score + my_penalty;
}

int Fitch::recompute() {
    assert(!ptn_freq.empty());
    assert(!node_major.empty());

    PhyloNode* root = (PhyloNode*)tree->root;
    PhyloNode* root_neighbor = (PhyloNode*)root->neighbors[0]->node;

    int subtree_score = localBottomUp(root_neighbor, root);

    int root_edge_score = 0;
    int root_idx = getIdx(root);
    int rn_idx = getIdx(root_neighbor);
    const nuc_one_hot* root_major_ptr = majorAt(root_idx);
    const nuc_one_hot* rn_major_ptr = majorAt(rn_idx);
    for (int ptn = 0; ptn < nptn; ptn++) {
        if ((root_major_ptr[ptn] & rn_major_ptr[ptn]) == 0) {
            root_edge_score += ptn_freq[ptn];
        }
    }

    int total_score = subtree_score + root_edge_score;

    vector<nuc_one_hot> virtual_root_states(nptn);
    for (int ptn = 0; ptn < nptn; ptn++) {
        nuc_one_hot intersect = root_major_ptr[ptn] & rn_major_ptr[ptn];
        if (intersect != 0) {
            virtual_root_states[ptn] = chooseRepresentative(intersect);
        } else {
            virtual_root_states[ptn] = chooseRepresentative(root_major_ptr[ptn] | rn_major_ptr[ptn]);
        }
    }

    topDown(root_neighbor, root, virtual_root_states);

    root_side_mutation_count = 0;
    {
        for (int ptn = 0; ptn < nptn; ptn++) {
            nuc_one_hot vr = virtual_root_states[ptn];

            nuc_one_hot root_assigned;
            if (vr & root_major_ptr[ptn]) {
                root_assigned = vr;
            } else {
                root_assigned = chooseRepresentative(root_major_ptr[ptn]);
            }

            if (root_assigned == vr) continue;
            if (ptn_is_const[ptn]) continue;

            root_side_mutation_count += (int)pattern_to_sites[ptn].size();
        }
    }

    computeFitchDiffs();

    return total_score;
}

int Fitch::recomputeScore() {
    assert(!ptn_freq.empty());
    assert(!node_major.empty());

    PhyloNode* root = (PhyloNode*)tree->root;
    PhyloNode* root_neighbor = (PhyloNode*)root->neighbors[0]->node;

    int subtree_score = localBottomUp(root_neighbor, root);

    int root_edge_score = 0;
    int root_idx = getIdx(root);
    int rn_idx = getIdx(root_neighbor);
    const nuc_one_hot* root_major_ptr = majorAt(root_idx);
    const nuc_one_hot* rn_major_ptr = majorAt(rn_idx);
    for (int ptn = 0; ptn < nptn; ptn++) {
        if ((root_major_ptr[ptn] & rn_major_ptr[ptn]) == 0) {
            root_edge_score += ptn_freq[ptn];
        }
    }

    return subtree_score + root_edge_score;
}

int Fitch::localBottomUpDirty(PhyloNode* node, PhyloNode* parent,
                                     bool& changed,
                                     std::vector<DirtyNodeSave>* dirty_saved,
                                     const std::set<PhyloNode*>* force_dirty) {
    int idx = getIdx(node);

    if (node->isLeaf()) {
        changed = false;
        return 0;
    }

    PhyloNode* child0 = nullptr;
    PhyloNode* child1 = nullptr;
    FOR_NEIGHBOR_IT(node, parent, it) {
        if (!child0) child0 = (PhyloNode*)(*it)->node;
        else child1 = (PhyloNode*)(*it)->node;
    }
    assert(child0 && child1);

    bool c0_changed = false, c1_changed = false;
    int s0 = localBottomUpDirty(child0, node, c0_changed, dirty_saved, force_dirty);
    int s1 = localBottomUpDirty(child1, node, c1_changed, dirty_saved, force_dirty);
    int children_score = s0 + s1;

    bool must_recompute = c0_changed || c1_changed;
    if (!must_recompute && force_dirty && force_dirty->count(node)) {
        must_recompute = true;
    }

    if (!must_recompute) {
        changed = false;
        return children_score + node_penalty[idx];
    }

    nuc_one_hot* my_major = majorAt(idx);
    const nuc_one_hot* c0 = majorAt(getIdx(child0));
    const nuc_one_hot* c1 = majorAt(getIdx(child1));

    if (dirty_saved) {
        DirtyNodeSave save;
        save.idx = idx;
        save.major.assign(my_major, my_major + nptn);
        save.penalty = node_penalty[idx];
        save.sub_score = subtree_score[idx];
        dirty_saved->push_back(std::move(save));
    }

    int my_penalty = 0;
    changed = false;
    for (int ptn = 0; ptn < nptn; ptn++) {
        nuc_one_hot intersect = c0[ptn] & c1[ptn];
        nuc_one_hot new_state = intersect ? intersect : (c0[ptn] | c1[ptn]);
        if (new_state != my_major[ptn]) changed = true;
        my_major[ptn] = new_state;
        if (!intersect) my_penalty += ptn_freq[ptn];
    }

    node_penalty[idx] = my_penalty;
    subtree_score[idx] = children_score + my_penalty;
    return children_score + my_penalty;
}

int Fitch::recomputeScoreDirty(const std::set<PhyloNode*>& force_dirty) {
    assert(!ptn_freq.empty());
    assert(!node_major.empty());
    assert(!node_penalty.empty());

    PhyloNode* root = (PhyloNode*)tree->root;
    PhyloNode* root_neighbor = (PhyloNode*)root->neighbors[0]->node;

    bool changed = false;
    int sub_score = localBottomUpDirty(root_neighbor, root, changed, nullptr, &force_dirty);

    int root_edge_score = 0;
    int root_idx = getIdx(root);
    const nuc_one_hot* root_major_ptr = majorAt(root_idx);
    const nuc_one_hot* rn_major_ptr = majorAt(getIdx(root_neighbor));
    for (int ptn = 0; ptn < nptn; ptn++) {
        if ((root_major_ptr[ptn] & rn_major_ptr[ptn]) == 0) {
            root_edge_score += ptn_freq[ptn];
        }
    }

    return sub_score + root_edge_score;
}

int Fitch::recomputeScoreDirtyAndRestore(const std::set<PhyloNode*>& force_dirty) {
    assert(!ptn_freq.empty());
    assert(!node_major.empty());
    assert(!node_penalty.empty());

    PhyloNode* root = (PhyloNode*)tree->root;
    PhyloNode* root_neighbor = (PhyloNode*)root->neighbors[0]->node;

    std::vector<DirtyNodeSave> dirty_saved;

    bool changed = false;
    int sub_score = localBottomUpDirty(root_neighbor, root, changed, &dirty_saved, &force_dirty);

    int root_edge_score = 0;
    int root_idx = getIdx(root);
    const nuc_one_hot* root_major_ptr = majorAt(root_idx);
    const nuc_one_hot* rn_major_ptr = majorAt(getIdx(root_neighbor));
    for (int ptn = 0; ptn < nptn; ptn++) {
        if ((root_major_ptr[ptn] & rn_major_ptr[ptn]) == 0) {
            root_edge_score += ptn_freq[ptn];
        }
    }

    int score = sub_score + root_edge_score;

    // Restore in reverse order (top-down) to match bottom-up save order
    for (auto it = dirty_saved.rbegin(); it != dirty_saved.rend(); ++it) {
        int idx = it->idx;
        nuc_one_hot* my_major = majorAt(idx);
        std::copy(it->major.begin(), it->major.end(), my_major);
        node_penalty[idx] = it->penalty;
        subtree_score[idx] = it->sub_score;
    }

    return score;
}
