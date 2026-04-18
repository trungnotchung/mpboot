#ifndef SPR_MUTATION_OPS_H
#define SPR_MUTATION_OPS_H

#include <vector>
#include <cstdint>
#include "onehot_encoding.h"
#include "phylotree.h"
#include "phylonode.h"
#include "mutation.h"

class MutationCountChange {
public:
    int position;
    uint8_t decremented;
    uint8_t incremented;
    uint8_t par_state;
    uint8_t major_allele;
    uint8_t boundary1_allele;
    bool from_src;

    MutationCountChange()
        : position(-1), decremented(0), incremented(0),
          par_state(0), major_allele(0), boundary1_allele(0), from_src(false) {}

    MutationCountChange(int pos, uint8_t dec, uint8_t inc)
        : position(pos), decremented(dec), incremented(inc),
          par_state(0), major_allele(0), boundary1_allele(0), from_src(false) {}

    MutationCountChange(int pos, uint8_t dec, uint8_t inc,
                       uint8_t par, uint8_t major, uint8_t boundary1)
        : position(pos), decremented(dec), incremented(inc),
          par_state(par), major_allele(major), boundary1_allele(boundary1), from_src(false) {}

    // Score change for intermediate nodes (on src->LCA or LCA->dst path).
    int get_default_change_internal() const {
        if (incremented && (incremented & par_state)) return -1;
        if (decremented && (decremented & par_state)) return +1;
        return 0;
    }

    // Score change for terminal nodes (src removal or dst insertion point).
    int get_default_change_terminal() const {
        if (incremented && !(incremented & par_state)) return +1;
        if (decremented && !(decremented & par_state)) return -1;
        return 0;
    }

    // Only recompute major allele set if the change affects boundary-1 alleles.
    bool is_sensitive() const {
        return (decremented & boundary1_allele) || (incremented & boundary1_allele);
    }

    bool has_change() const { return decremented != 0 || incremented != 0; }
    uint8_t get_net_change() const { return decremented | incremented; }

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

class Fitch;

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

    /// Set the root node pointer. Must be called after orientTreeToRoot().
    static void setRoot(PhyloNode* root);

    /// Set the Fitch pointer. Must be called after Fitch::run().
    static void setCustomFitch(Fitch* cf);

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
