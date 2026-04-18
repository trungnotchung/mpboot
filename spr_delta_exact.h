#ifndef SPR_DELTA_EXACT_H
#define SPR_DELTA_EXACT_H

#include "spr_mutation_ops.h"
#include "phylotree.h"
#include "phylonode.h"

class Fitch;

class SPRDeltaExact {
public:
    // Set Fitch pointer and current score for binary node fallback.
    static void setCustomFitch(Fitch* cf, int current_score);

    // Compute exact parsimony delta for SPR move.
    static int calculateParsimonyDelta(PhyloNode* src, PhyloNode* src_parent,
                                       PhyloNode* dst, PhyloNode* dst_parent,
                                       PhyloNode* lca,
                                       const std::map<PhyloNode*, PhyloNode*>* node_parent = nullptr);

    // Compute delta with automatic LCA finding.
    static int calculateParsimonyDelta(PhyloNode* src, PhyloNode* src_parent,
                                       PhyloNode* dst, PhyloNode* dst_parent,
                                       PhyloTree* tree);

    // Find LCA of two nodes.
    static PhyloNode* findLCA(PhyloNode* node1, PhyloNode* node2, PhyloTree* tree);

public:
    // Source-side delta (removing src from src_parent to LCA).
    static int computeSrcSideDelta(PhyloNode* src, PhyloNode* src_parent,
                                   PhyloNode* lca,
                                   MutationCountChangeCollection& mutations_out,
                                   MutationCountChangeCollection& root_mutations_altered,
                                   const std::map<PhyloNode*, PhyloNode*>* node_parent);

    // Destination-side delta (inserting src at dst).
    static int computeDstSideDelta(PhyloNode* src, PhyloNode* dst,
                                   PhyloNode* dst_parent, PhyloNode* lca,
                                   const MutationCountChangeCollection& mutations,
                                   const MutationCountChangeCollection& root_mutations_altered,
                                   const std::map<PhyloNode*, PhyloNode*>* node_parent,
                                   MutationCountChangeCollection* dst_added = nullptr);

    // Above-LCA propagation (LCA to root).
    static int computeAboveLCADelta(PhyloNode* lca,
                                    const MutationCountChangeCollection& lca_changes);

    // Combined LCA merge and above-LCA propagation.
    static int computeLCAAndAboveDelta(PhyloNode* lca, PhyloNode* src,
                                       const MutationCountChangeCollection& src_changes,
                                       const MutationCountChangeCollection& dst_changes);

    // Check if move is profitable when dst == LCA.
    static int checkMoveProfitableLCA(PhyloNode* src, PhyloNode* lca,
                                      const MutationCountChangeCollection& mutations,
                                      const MutationCountChangeCollection& root_mutations,
                                      int base_score_change,
                                      MutationCountChangeCollection* dst_added = nullptr);

    // Check if move is profitable when dst != LCA.
    static int checkMoveProfitableDstNotLCA(PhyloNode* src, PhyloNode* dst,
                                            PhyloNode* lca,
                                            const MutationCountChangeCollection& mutations,
                                            const MutationCountChangeCollection& root_mutations,
                                            int base_score_change,
                                            MutationCountChangeCollection* dst_added = nullptr);

    // Print O(M) sparse delta statistics and reset counters.
    static void printSparseStats();

    // Precompute node depths for O(1) lookup in findLCA.
    static void precomputeDepths(PhyloTree* tree);
};

#endif // SPR_DELTA_EXACT_H
