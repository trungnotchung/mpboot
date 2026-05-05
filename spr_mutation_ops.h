#ifndef SPR_MUTATION_OPS_H
#define SPR_MUTATION_OPS_H

#include <vector>
#include <cstdint>
#include "nucleotide_utils.h"
#include "phylotree.h"
#include "phylonode.h"
#include "mutation.h"

/**
 * One pending allele-count change at a single genomic position.
 */
class MutationCountChange {
public:
    int position;                    ///< Genomic position (or pattern index).
    nuc_one_hot removed_alleles;     ///< Alleles whose count went down.
    nuc_one_hot added_alleles;       ///< Alleles whose count went up.
    nuc_one_hot par_state;           ///< Parent's allele state.
    nuc_one_hot major_allele_set;    ///< Alleles tied for max count.
    nuc_one_hot boundary1_allele;    ///< Alleles at count = max_count - 1.
    bool from_src;                   ///< True if from src removal.

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

    /** Score change for intermediate nodes (-1, 0, or +1). */
    int scoreDeltaInternal() const {
        if (added_alleles && (added_alleles & par_state)) return -1;
        if (removed_alleles && (removed_alleles & par_state)) return +1;
        return 0;
    }

    /** Score change for terminal nodes (-1, 0, or +1). */
    int scoreDeltaTerminal() const {
        if (added_alleles && !(added_alleles & par_state)) return +1;
        if (removed_alleles && !(removed_alleles & par_state)) return -1;
        return 0;
    }

    /** True if change may flip the major allele set. */
    bool isSensitive() const {
        return (removed_alleles & boundary1_allele) || (added_alleles & boundary1_allele);
    }

    /** True if any allele was added or removed. */
    bool hasChange() const { return removed_alleles != 0 || added_alleles != 0; }

    /** Union of removed and added alleles. */
    nuc_one_hot affectedAlleles() const { return removed_alleles | added_alleles; }

    bool operator<(const MutationCountChange& other) const { return position < other.position; }
    bool operator<(int pos) const { return position < pos; }
    friend bool operator<(int pos, const MutationCountChange& mcc) { return pos < mcc.position; }
    bool operator==(const MutationCountChange& other) const { return position == other.position; }
};

using MutationCountChangeCollection = std::vector<MutationCountChange>;

namespace MutationCountChangeUtils {
    /**
     * Merges two sorted-by-position collections.
     * @param a First sorted collection.
     * @param b Second sorted collection.
     * @return Merged sorted collection.
     */
    MutationCountChangeCollection merge_sorted(
        const MutationCountChangeCollection& a,
        const MutationCountChangeCollection& b);

    /**
     * Binary search for first entry with given position.
     * @param collection Sorted collection.
     * @param position   Position to find.
     * @return Iterator to entry, or end() if not found.
     */
    MutationCountChangeCollection::const_iterator find_position(
        const MutationCountChangeCollection& collection,
        int position);

    /**
     * @return true if collection is sorted by position.
     */
    bool is_sorted(const MutationCountChangeCollection& collection);
}

/**
 * Mutation-based SPR delta operations. Requires fresh PhyloNeighbor::mutations.
 */
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
