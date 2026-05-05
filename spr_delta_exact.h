#ifndef SPR_DELTA_EXACT_H
#define SPR_DELTA_EXACT_H

#include "spr_mutation_ops.h"
#include "phylotree.h"
#include "phylonode.h"

/**
 * Exact SPR delta computation. Uses sparse Fitch path for binary nodes,
 * mutation-based fallback otherwise.
 */
class SPRDeltaExact {
public:
    /**
     * Computes parsimony delta for an SPR move with given LCA.
     * @param src         Source node.
     * @param src_parent  Source's parent.
     * @param dst         Regraft target.
     * @param dst_parent  Parent of dst.
     * @param lca         LCA of src_parent and dst.
     * @param node_parent Optional parent map (for mutation-based path).
     * @return Score delta.
     */
    static int calculateParsimonyDelta(PhyloNode* src, PhyloNode* src_parent,
                                       PhyloNode* dst, PhyloNode* dst_parent,
                                       PhyloNode* lca,
                                       const std::map<PhyloNode*, PhyloNode*>* node_parent = nullptr);

    /**
     * Computes parsimony delta with auto LCA via tree's LCA table.
     * @param src        Source node.
     * @param src_parent Source's parent.
     * @param dst        Regraft target.
     * @param dst_parent Parent of dst.
     * @param tree       Tree containing precomputed LCA table.
     * @return Score delta.
     */
    static int calculateParsimonyDelta(PhyloNode* src, PhyloNode* src_parent,
                                       PhyloNode* dst, PhyloNode* dst_parent,
                                       PhyloTree* tree);

    /**
     * Finds LCA of two nodes.
     * @param node1 First node.
     * @param node2 Second node.
     * @param tree  Tree with LCA table.
     * @return Lowest common ancestor.
     */
    static PhyloNode* findLCA(PhyloNode* node1, PhyloNode* node2, PhyloTree* tree);

public:
    /**
     * Source-side delta: removing src from src_parent up to LCA.
     * @param src                    Source node.
     * @param src_parent             Source's parent.
     * @param lca                    LCA (stop point).
     * @param mutations_out          Output: changes propagating to dst-side.
     * @param root_mutations_altered Output: root-edge effects.
     * @param node_parent            Parent map.
     * @return Score delta along src->LCA path.
     */
    static int computeSrcSideDelta(PhyloNode* src, PhyloNode* src_parent,
                                   PhyloNode* lca,
                                   MutationCountChangeCollection& mutations_out,
                                   MutationCountChangeCollection& root_mutations_altered,
                                   const std::map<PhyloNode*, PhyloNode*>* node_parent);

    /**
     * Destination-side delta: inserting src on the (dst, dst_parent) edge.
     * @param src                    Source node.
     * @param dst                    Regraft target.
     * @param dst_parent             Parent of dst.
     * @param lca                    LCA.
     * @param mutations              Changes from src side.
     * @param root_mutations_altered Root-edge effects from src side.
     * @param node_parent            Parent map.
     * @param dst_added              Optional output: changes added at dst.
     * @return Score delta along LCA->dst path.
     */
    static int computeDstSideDelta(PhyloNode* src, PhyloNode* dst,
                                   PhyloNode* dst_parent, PhyloNode* lca,
                                   const MutationCountChangeCollection& mutations,
                                   const MutationCountChangeCollection& root_mutations_altered,
                                   const std::map<PhyloNode*, PhyloNode*>* node_parent,
                                   MutationCountChangeCollection* dst_added = nullptr);

    /**
     * Propagates LCA's net allele changes from LCA upward to root.
     * @param lca         The LCA node.
     * @param lca_changes Net changes accumulated at LCA.
     * @return Score delta on LCA->root path.
     */
    static int computeAboveLCADelta(PhyloNode* lca,
                                    const MutationCountChangeCollection& lca_changes);

    /**
     * Combined LCA merge plus above-LCA propagation.
     * @param lca         LCA node.
     * @param src         Source node (for direction).
     * @param src_changes Changes arriving from src side.
     * @param dst_changes Changes arriving from dst side.
     * @return Score delta from merge plus above-LCA propagation.
     */
    static int computeLCAAndAboveDelta(PhyloNode* lca, PhyloNode* src,
                                       const MutationCountChangeCollection& src_changes,
                                       const MutationCountChangeCollection& dst_changes);

    /**
     * Profitability check when dst == LCA.
     * @param src               Source node.
     * @param lca               LCA == dst.
     * @param mutations         Src-side changes.
     * @param root_mutations    Root-edge changes.
     * @param base_score_change Score change accumulated so far.
     * @param dst_added         Optional output: changes added at dst.
     * @return Total score delta.
     */
    static int checkMoveProfitableLCA(PhyloNode* src, PhyloNode* lca,
                                      const MutationCountChangeCollection& mutations,
                                      const MutationCountChangeCollection& root_mutations,
                                      int base_score_change,
                                      MutationCountChangeCollection* dst_added = nullptr);

    /**
     * Profitability check when dst != LCA.
     * @param src               Source node.
     * @param dst               Regraft target.
     * @param lca               LCA.
     * @param mutations         Src-side changes.
     * @param root_mutations    Root-edge changes.
     * @param base_score_change Score change accumulated so far.
     * @param dst_added         Optional output: changes added at dst.
     * @return Total score delta.
     */
    static int checkMoveProfitableDstNotLCA(PhyloNode* src, PhyloNode* dst,
                                            PhyloNode* lca,
                                            const MutationCountChangeCollection& mutations,
                                            const MutationCountChangeCollection& root_mutations,
                                            int base_score_change,
                                            MutationCountChangeCollection* dst_added = nullptr);

    /**
     * Prints sparse-delta statistics and resets counters.
     */
    static void printSparseStats();

    /**
     * Builds the LCA table on tree (call after every topology change).
     * @param tree Tree to build LCA table on.
     */
    static void precomputeDepths(PhyloTree* tree);
};

#endif // SPR_DELTA_EXACT_H
