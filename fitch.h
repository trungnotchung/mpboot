#ifndef FITCH_H
#define FITCH_H

#include <vector>
#include <map>
#include <set>
#include <cstdint>
#include "mutation.h"

class PhyloTree;
class PhyloNode;

class Fitch {
public:
    Fitch(PhyloTree* tree);

    // Full Fitch: reads alignment, populates mutations, returns parsimony score.
    int run();

    // Recompute after topology change (no alignment access). Must call run() first.
    int recompute();

    // Score-only recompute: skips mutation generation.
    int recomputeScore();

    int countMutations() const;

    // Access: node_major[getIdx(node) * nptn + pattern]
    const std::vector<nuc_one_hot>& getNodeMajor() const { return node_major; }

    const std::vector<std::vector<int>>& getPatternToSites() const { return pattern_to_sites; }

    int getNumPatterns() const { return nptn; }

    int getPatternFreq(int ptn) const { return ptn_freq[ptn]; }

    // Returns -1 if no mutation exists at this pattern.
    int getPositionForPattern(int ptn) const { return ptn_position[ptn]; }

    // Must be called AFTER mutations are populated on the tree edges.
    void buildPatternPositionMap();

    const nuc_one_hot* getMajorArrayForNode(PhyloNode* node) const {
        auto it = node_index.find(node);
        if (it == node_index.end()) return nullptr;
        return &node_major[it->second * nptn];
    }

    // Sorted pattern indices where this node's Fitch set differs from parent's.
    const std::vector<int>* getFitchDiffs(PhyloNode* node) const {
        auto it = node_index.find(node);
        if (it == node_index.end()) return nullptr;
        return &fitch_diffs[it->second];
    }

    nuc_one_hot getMajorForNode(PhyloNode* node, int ptn) const {
        auto it = node_index.find(node);
        if (it == node_index.end()) return 0;
        return node_major[it->second * nptn + ptn];
    }

    std::vector<nuc_one_hot> saveNodeMajor() const { return node_major; }
    void restoreNodeMajor(const std::vector<nuc_one_hot>& saved) { node_major = saved; }

    // Score recompute skipping clean subtrees. Modifies node_major for dirty nodes.
    int recomputeScoreDirty(const std::set<PhyloNode*>& force_dirty);

    // Like recomputeScoreDirty() but auto-saves/restores node_major.
    int recomputeScoreDirtyAndRestore(const std::set<PhyloNode*>& force_dirty);

    const std::map<PhyloNode*, int>& getNodeIndex() const { return node_index; }

private:
    PhyloTree* tree;

    // node_major[index * nptn + pattern] = one-hot Fitch state set
    std::vector<nuc_one_hot> node_major;

    // std::map used instead of unordered_map due to ext/hash_map conflicts in this codebase
    std::map<PhyloNode*, int> node_index;

    std::vector<std::vector<int>> pattern_to_sites;

    // Cached pattern metadata (avoids alignment access in recompute())
    std::vector<int> ptn_freq;
    std::vector<bool> ptn_is_const;
    std::vector<int> ptn_position;

    int nptn;
    int num_nodes;

    // Root-side mutation count: mutations between root leaf and virtual root.
    // Stored separately (not on the root edge) so SPRDeltaExact sees clean mutations.
    // Added to countMutations() total.
    int root_side_mutation_count;

    // Per-node cached penalty and subtree score (for dirty recompute optimization)
    std::vector<int> node_penalty;
    std::vector<int> subtree_score;

    std::vector<std::vector<int>> fitch_diffs;

    inline int getIdx(PhyloNode* node) const;
    inline nuc_one_hot* majorAt(int idx) { return &node_major[idx * nptn]; }
    inline const nuc_one_hot* majorAt(int idx) const { return &node_major[idx * nptn]; }

    int bottomUp(PhyloNode* node, PhyloNode* parent);

    // Like bottomUp but uses existing leaf node_major (no alignment access).
    int localBottomUp(PhyloNode* node, PhyloNode* parent);

    struct DirtyNodeSave {
        int idx;
        std::vector<nuc_one_hot> major;
        int penalty;
        int sub_score;
    };

    // Bottom-up Fitch with dirty-path optimization (clean subtrees use cached scores).
    int localBottomUpDirty(PhyloNode* node, PhyloNode* parent,
                           bool& changed,
                           std::vector<DirtyNodeSave>* dirty_saved,
                           const std::set<PhyloNode*>* force_dirty);

    void topDown(PhyloNode* node, PhyloNode* parent,
                 const std::vector<nuc_one_hot>& parent_states);
    void buildPatternToSites();
    void computeFitchDiffs();

    // Picks lowest bit (A > C > G > T priority). Returns 0xF for gap/empty.
    static nuc_one_hot chooseRepresentative(nuc_one_hot state_set);
};

#endif // FITCH_H
