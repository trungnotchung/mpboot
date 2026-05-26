#ifndef SPR_UTILS_H
#define SPR_UTILS_H

#include "mutation.h"
#include "phylonode.h"
#include <vector>
#include <algorithm>

/** Binary internal node degree: 1 parent + 2 children. */
static const int BINARY_NODE_DEGREE = 3;

/**
 * Snapshot of a single neighbor pointer for SPR undo.
 */
struct SPRNeighborSave {
    Neighbor* neighbor;
    Node*     original_node;
};

/**
 * Records every (Neighbor*, neighbor->node) pair on a node.
 * @param node  Node whose neighbors to snapshot.
 * @param saves Output collection (appended to).
 */
inline void sprSaveNodeTopology(Node* node, std::vector<SPRNeighborSave>& saves) {
    for (auto it = node->neighbors.begin(); it != node->neighbors.end(); ++it) {
        saves.push_back({*it, (*it)->node});
    }
}

/**
 * Restores neighbor pointers to their saved values.
 * @param saves Saves recorded by sprSaveNodeTopology.
 */
inline void sprUndoTopology(std::vector<SPRNeighborSave>& saves) {
    for (auto& s : saves) {
        s.neighbor->node = s.original_node;
    }
}

/**
 * Fitch merge: intersection if non-empty, else union.
 * @param left  Left child Fitch state.
 * @param right Right child Fitch state.
 * @return Merged Fitch state.
 */
static inline nuc_one_hot fitchMerge(nuc_one_hot left, nuc_one_hot right) {
    nuc_one_hot intersect = left & right;
    return intersect ? intersect : (left | right);
}

/**
 * Score-delta of replacing one child's Fitch state with another at a node.
 * @param old_child    Replaced child's Fitch state.
 * @param new_child    Replacement child's Fitch state.
 * @param other        Sibling's Fitch state.
 * @param new_node_out Output: merged state after the swap.
 * @return Change in mutation count (-1, 0, or +1).
 */
static inline int penaltyDelta(nuc_one_hot old_child, nuc_one_hot new_child,
                                nuc_one_hot other, nuc_one_hot& new_node_out) {
    int old_penalty = (old_child & other) ? 0 : 1;
    nuc_one_hot intersect = new_child & other;
    new_node_out = intersect ? intersect : (new_child | other);
    int new_penalty = intersect ? 0 : 1;
    return new_penalty - old_penalty;
}

/**
 * Score-delta at the unrooted-tree's virtual root edge.
 * @param root_fitch         Root leaf's Fitch state.
 * @param old_neighbor_fitch Old neighbor's Fitch state.
 * @param new_neighbor_fitch New neighbor's Fitch state.
 * @return Change in root-edge mutation count.
 */
static inline int rootEdgeDelta(nuc_one_hot root_fitch, nuc_one_hot old_neighbor_fitch,
                                 nuc_one_hot new_neighbor_fitch) {
    int old_penalty = (root_fitch & old_neighbor_fitch) ? 0 : 1;
    int new_penalty = (root_fitch & new_neighbor_fitch) ? 0 : 1;
    return new_penalty - old_penalty;
}

/**
 * Finds the neighbor of `node` that is neither parent nor known_child.
 * @param node        Node whose neighbors to search.
 * @param parent      Neighbor to exclude (parent direction).
 * @param known_child Neighbor to exclude (known child).
 * @return The third neighbor, or nullptr.
 */
static inline PhyloNode* findOtherChild(PhyloNode* node, PhyloNode* parent,
                                         PhyloNode* known_child) {
    if (!node) return nullptr;
    FOR_NEIGHBOR_IT(node, parent, nit) {
        PhyloNode* child = (PhyloNode*)(*nit)->node;
        if (child != known_child) return child;
    }
    return nullptr;
}

/**
 * Merges a sorted-unique diff list into out (also sorted-unique on return).
 * @param out   Existing sorted-unique vector; modified in place.
 * @param diffs Sorted-unique source vector to merge in (may be null).
 */
static inline void mergeDiffsSorted(std::vector<int>& out, const std::vector<int>* diffs) {
    if (!diffs || diffs->empty()) return;
    if (out.empty()) { out = *diffs; return; }
    size_t old_size = out.size();
    out.insert(out.end(), diffs->begin(), diffs->end());
    std::inplace_merge(out.begin(), out.begin() + old_size, out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
}

/**
 * One node along a propagation path (for delta evaluation).
 */
struct PathStep {
    PhyloNode* node;
    const nuc_one_hot* sibling_states;
    const nuc_one_hot* node_states;
};

/**
 * Propagates Fitch state changes upward along a path, accumulating delta.
 * Stops early when old_fitch == new_fitch (change absorbed).
 * @param path      Sequence of PathSteps.
 * @param start     Inclusive start index.
 * @param end       Exclusive end index.
 * @param ptn       Pattern index.
 * @param freq      Pattern frequency multiplier.
 * @param old_fitch In/out: previous state propagating up.
 * @param new_fitch In/out: new state propagating up.
 * @param score     In/out: accumulator for delta * freq.
 */
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

/**
 * Convenience overload: propagates the entire path from start=0.
 */
static inline void propagatePath(const std::vector<PathStep>& path,
                                  int ptn, int freq,
                                  nuc_one_hot& old_fitch, nuc_one_hot& new_fitch,
                                  int& score) {
    propagatePath(path, 0, path.size(), ptn, freq, old_fitch, new_fitch, score);
}

#endif // SPR_UTILS_H
