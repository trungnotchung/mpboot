#ifndef SPR_MUTATION_OPS_H
#define SPR_MUTATION_OPS_H

#include <vector>
#include <cstdint>
#include "nucleotide_utils.h"
#include "phylotree.h"
#include "phylonode.h"
#include "mutation.h"

class MutationCountChange {
public:
    int position;
    nuc_one_hot removed_alleles;
    nuc_one_hot added_alleles;
    nuc_one_hot par_state;
    nuc_one_hot major_allele_set;
    nuc_one_hot boundary1_allele;
    bool from_src;

    MutationCountChange()
        : position(-1), removed_alleles(0), added_alleles(0),
          par_state(0), major_allele_set(0), boundary1_allele(0), from_src(false) {}

    MutationCountChange(int pos, nuc_one_hot dec, nuc_one_hot inc)
        : position(pos), removed_alleles(dec), added_alleles(inc),
          par_state(0), major_allele_set(0), boundary1_allele(0), from_src(false) {}

    MutationCountChange(int pos, nuc_one_hot dec, nuc_one_hot inc,
                       nuc_one_hot par, nuc_one_hot major, nuc_one_hot boundary1)
        : position(pos), removed_alleles(dec), added_alleles(inc),
          par_state(par), major_allele_set(major), boundary1_allele(boundary1), from_src(false) {}

    // Score change for intermediate nodes (on src->LCA or LCA->dst path).
    int scoreDeltaInternal() const {
        if (added_alleles && (added_alleles & par_state)) return -1;
        if (removed_alleles && (removed_alleles & par_state)) return +1;
        return 0;
    }

    // Score change for terminal nodes (src removal or dst insertion point).
    int scoreDeltaTerminal() const {
        if (added_alleles && !(added_alleles & par_state)) return +1;
        if (removed_alleles && !(removed_alleles & par_state)) return -1;
        return 0;
    }

    // Only recompute major allele set if the change affects boundary-1 alleles.
    bool isSensitive() const {
        return (removed_alleles & boundary1_allele) || (added_alleles & boundary1_allele);
    }

    bool hasChange() const { return removed_alleles != 0 || added_alleles != 0; }
    nuc_one_hot affectedAlleles() const { return removed_alleles | added_alleles; }

    bool operator<(const MutationCountChange& other) const { return position < other.position; }
    bool operator<(int pos) const { return position < pos; }
    friend bool operator<(int pos, const MutationCountChange& mcc) { return pos < mcc.position; }
    bool operator==(const MutationCountChange& other) const { return position == other.position; }
};

using MutationCountChangeCollection = std::vector<MutationCountChange>;

namespace MutationCountChangeUtils {
    MutationCountChangeCollection merge_sorted(
        const MutationCountChangeCollection& a,
        const MutationCountChangeCollection& b);

    MutationCountChangeCollection::const_iterator find_position(
        const MutationCountChangeCollection& collection,
        int position);

    bool is_sorted(const MutationCountChangeCollection& collection);
}

class SPRMutationOps {
public:
    /// Initialize src's mutations as count changes.
    static MutationCountChangeCollection initMutationChange(PhyloNode* src, PhyloNode* src_parent);

    /// Merge src's mutations upward through ancestor (src_parent to LCA).
    static void mergeMutationSrcToLCA(PhyloNode* ancestor,
                                      MutationCountChangeCollection& mutations);

    /// Merge mutations downward (LCA to dst).
    static MutationCountChangeCollection mergeMutationLCAToRank(
        PhyloNode* child_on_path,
        const MutationCountChangeCollection& mutations);

    /// Compute effect of removing src from its parent.
    static void getParentAlteredRemove(MutationCountChangeCollection& output,
                                       PhyloNode* src,
                                       int& parsimony_score_change);

    /// Propagate changes through intermediate nodes (src_parent to LCA).
    static void getIntermediateNodesMutations(
        PhyloNode* node,
        const MutationCountChangeCollection& child_changes,
        MutationCountChangeCollection& parent_changes,
        int& score_change,
        PhyloNode* child_on_path = nullptr);

    /// Get parent node (first neighbor), or nullptr if root.
    static PhyloNode* getParent(PhyloNode* node);

    /// Get mutations on edge between node and dad (nullptr = parent edge).
    static std::vector<Mutation>* getMutations(PhyloNode* node, PhyloNode* dad);

    /// Get src's one-hot allele at a specific position (0 if no mutation).
    static uint8_t getSrcAlleleAtPosition(PhyloNode* src, int position);

    /// Get sibling's one-hot allele at a specific position.
    static uint8_t getSiblingAllele(PhyloNode* parent, PhyloNode* src, int position);

    /// Recompute major allele set for a node given allele count changes.
    static int recomputeMajorAllele(PhyloNode* node, int position,
                                    uint8_t dec_allele, uint8_t inc_allele,
                                    uint8_t& new_major, uint8_t& new_boundary1,
                                    int* children_node_delta = nullptr);

    /// Compute full parsimony delta (children-to-node + node-to-parent edge change).
    static int recomputeFullDelta(PhyloNode* node, int position,
                                  uint8_t dec_allele, uint8_t inc_allele,
                                  uint8_t change_par_state = 0);

    /// Propagate allele changes from LCA upward toward root.
    static void checkParsimonyScoreChangeAboveLCA(
        PhyloNode* start_node,
        int& parsimony_score_change,
        const MutationCountChangeCollection& allele_changes,
        MutationCountChangeCollection& parent_changes);
};

#endif // SPR_MUTATION_OPS_H
