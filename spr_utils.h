#ifndef SPR_UTILS_H
#define SPR_UTILS_H

#include "mutation.h"
#include "phylonode.h"
#include <vector>
#include <algorithm>

struct SPRNeighborSave {
    Neighbor* neighbor;
    Node*     original_node;
};

inline void sprSaveNodeTopology(Node* node, std::vector<SPRNeighborSave>& saves) {
    for (auto it = node->neighbors.begin(); it != node->neighbors.end(); ++it) {
        saves.push_back({*it, (*it)->node});
    }
}

inline void sprUndoTopology(std::vector<SPRNeighborSave>& saves) {
    for (auto& s : saves) {
        s.neighbor->node = s.original_node;
    }
}

static inline nuc_one_hot fitchMerge(nuc_one_hot left, nuc_one_hot right) {
    nuc_one_hot intersect = left & right;
    return intersect ? intersect : (left | right);
}

static inline int penaltyDelta(nuc_one_hot old_child, nuc_one_hot new_child,
                                nuc_one_hot other, nuc_one_hot& new_node_out) {
    int old_penalty = (old_child & other) ? 0 : 1;
    nuc_one_hot intersect = new_child & other;
    new_node_out = intersect ? intersect : (new_child | other);
    int new_penalty = intersect ? 0 : 1;
    return new_penalty - old_penalty;
}

static inline int rootEdgeDelta(nuc_one_hot root_fitch, nuc_one_hot old_neighbor_fitch,
                                 nuc_one_hot new_neighbor_fitch) {
    int old_penalty = (root_fitch & old_neighbor_fitch) ? 0 : 1;
    int new_penalty = (root_fitch & new_neighbor_fitch) ? 0 : 1;
    return new_penalty - old_penalty;
}

static inline PhyloNode* findOtherChild(PhyloNode* node, PhyloNode* parent,
                                         PhyloNode* known_child) {
    if (!node) return nullptr;
    FOR_NEIGHBOR_IT(node, parent, nit) {
        PhyloNode* child = (PhyloNode*)(*nit)->node;
        if (child != known_child) return child;
    }
    return nullptr;
}

static inline void mergeDiffsSorted(std::vector<int>& out, const std::vector<int>* diffs) {
    if (!diffs || diffs->empty()) return;
    if (out.empty()) { out = *diffs; return; }
    size_t old_size = out.size();
    out.insert(out.end(), diffs->begin(), diffs->end());
    std::inplace_merge(out.begin(), out.begin() + old_size, out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
}

struct PathStep {
    PhyloNode* node;
    const nuc_one_hot* sibling_states;  // Fitch states of sibling at this node
    const nuc_one_hot* node_states;     // Fitch states of this node
};

// Propagate Fitch state changes along a path, accumulating score delta.
// old_fitch/new_fitch track the propagating state and are modified in place.
static inline void propagatePath(const std::vector<PathStep>& path,
                                  size_t start, size_t end,
                                  int ptn, int freq,
                                  nuc_one_hot& old_fitch, nuc_one_hot& new_fitch,
                                  int& score) {
    for (size_t i = start; i < end && i < path.size(); i++) {
        if (old_fitch == new_fitch) break;
        const auto& step = path[i];
        nuc_one_hot updated_fitch;
        score += penaltyDelta(old_fitch, new_fitch, step.sibling_states[ptn], updated_fitch) * freq;
        old_fitch = step.node_states[ptn];
        new_fitch = updated_fitch;
    }
}

// Convenience overload: propagate entire path from start=0.
static inline void propagatePath(const std::vector<PathStep>& path,
                                  int ptn, int freq,
                                  nuc_one_hot& old_fitch, nuc_one_hot& new_fitch,
                                  int& score) {
    propagatePath(path, 0, path.size(), ptn, freq, old_fitch, new_fitch, score);
}

#endif // SPR_UTILS_H
