//
// C++ Interface: phylonode
//
// Description:
//
//
// Author: BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PHYLONODE_H
#define PHYLONODE_H

#include "node.h"
#include "mutation.h"
#include "nucleotide_utils.h"

typedef short int UBYTE;

/**
A neighbor in a phylogenetic tree

    @author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>
 */
class PhyloNeighbor : public Neighbor
{
    friend class PhyloNode;
    friend class PhyloTree;
    friend class IQTree;
    friend class PhyloSuperTree;
    friend class PlacementOptimizer;  // Needs access to partial_pars for exact scoring

public:
    friend class TinaTree;
    friend class PhyloSuperTreePlen;
    friend class ParsTree; // DTH ##############

    /**
        construct class with a node and length		@param anode the other end of the branch

        @param alength length of branch
     */
    PhyloNeighbor(Node *anode, double alength) : Neighbor(anode, alength) {
        partial_lh = NULL;
        partial_lh_computed = 0;
        lh_scale_factor = 0.0;
        partial_pars = NULL;
    }

    /**
        construct class with a node and length
        @param anode the other end of the branch
        @param alength length of branch
        @param aid branch ID
     */
    PhyloNeighbor(Node *anode, double alength, int aid) : Neighbor(anode, alength, aid) {
        partial_lh = NULL;
        partial_lh_computed = 0;
        lh_scale_factor = 0.0;
        partial_pars = NULL;
    }

    /**
        tell that the partial likelihood vector is not computed
     */
    inline void clearPartialLh() {
        partial_lh_computed = 0;
    }

    /**
     *  tell that the partial likelihood vector is computed
     */
    inline void unclearPartialLh() {
        partial_lh_computed = 1;
    }

    /**
        clear all partial likelihood recursively in forward direction
        @param dad dad of this neighbor
     */
    void clearForwardPartialLh(Node *dad);


    /**
     * All mutations on this branch
     */
    std::vector<Mutation> mutations;

    /**
     * Number of leaves in the subtree rooted at this node
     */
    int num_leaves;

    /**
     * Clear all mutations on this branch
     */
    void clearMutations();
private:
    /**
        true if the partial likelihood was computed
     */
    int partial_lh_computed;

    /**
        vector containing the partial likelihoods
     */
    double *partial_lh;

    /**
        likelihood scaling factor
     */
    double lh_scale_factor;

    /**
        vector containing number of scaling events per pattern // NEW!
     */
    UBYTE *scale_num;

    /**
        vector containing the partial parsimony scores
     */
    UINT *partial_pars;
};

/**
A node in a phylogenetic tree

    @author BUI Quang Minh, Steffen Klaere, Arndt von Haeseler <minh.bui@univie.ac.at>
 */
class PhyloNode : public Node
{
    friend class PhyloTree;

public:
    /**
        constructor
     */
    PhyloNode();

    /**
        constructor
        @param aid id of this node
     */
    PhyloNode(int aid);

    /**
        constructor
        @param aid id of this node
        @param aname name of this node
     */
    PhyloNode(int aid, int aname);

    /**
        constructor
        @param aid id of this node
        @param aname name of this node
     */
    PhyloNode(int aid, const char *aname);

    /**
        initialization
     */
    void init();

    void setMissingNode(int index);

    bool checkMissingNode();

    int getMissingIndex();

    /**
        add a neighbor
        @param node the neighbor node
        @param length branch length
        @param id branch ID
     */
    virtual void addNeighbor(Node *node, double length, int id = -1);

    /**
        tell that all partial likelihood vectors below this node are not computed
     */
    void clearAllPartialLh(PhyloNode *dad);

    /**
        tell that all partial likelihood vectors (in reverse direction) below this node are not computed
     */
    void clearReversePartialLh(PhyloNode *dad);

    PhyloNode *dad;

    int missingIndex;

    /**
     * DFS index for ordering nodes in mutation state reassignment
     * Used by backward pass (max-heap) and forward pass (min-heap)
     */
    int dfs_index;

    // Incremental Mutation Update State Storage

    /**
     * Current nucleotide state assignment at each position
     * Key: alignment position
     * Value: assigned state (one-hot encoded: 1=A, 2=C, 4=G, 8=T)
     */
    std::map<int, nuc_one_hot> assigned_states;

    /**
     * Major allele set from Fitch bottom-up pass at each position
     * Key: alignment position
     * Value: major allele set (bit-packed, multiple bits set if ambiguous)
     */
    std::map<int, nuc_one_hot> major_alleles;

    /**
     * Quick lookup: position → mutations on edges from this node
     * Used for fast mutation updates during incremental reassignment
     */
    std::map<int, std::vector<Mutation*>> mutations_by_position;

    /**
     * Helper: Get assigned state at a position
     */
    nuc_one_hot getStateAt(int position) const {
        auto it = assigned_states.find(position);
        return (it != assigned_states.end()) ? it->second : NUC_N;
    }

    /**
     * Helper: Set assigned state at a position
     */
    void setStateAt(int position, nuc_one_hot state) {
        assigned_states[position] = state;
    }

    /**
     * Helper: Get major allele at a position
     */
    nuc_one_hot getMajorAlleleAt(int position) const {
        auto it = major_alleles.find(position);
        return (it != major_alleles.end()) ? it->second : NUC_N;
    }

    /**
     * Helper: Build mutation position map from edge mutations
     * Call after mutations are created/updated
     */
    void buildMutationPositionMap();

    /**
     * Helper: Clear incremental state (for cleanup)
     */
    void clearIncrementalState() {
        assigned_states.clear();
        major_alleles.clear();
        mutations_by_position.clear();
    }
};

/**
    Node vector
 */
typedef vector<PhyloNode *> PhyloNodeVector;

class PlacementCandidateNode
{
public:
    PhyloNode *node;
    PhyloNeighbor *node_branch;
    std::vector<Mutation> *missing_sample_mutations;
    std::vector<Mutation> *excess_mutations;

    int *best_set_difference;
    size_t *best_node_num_leaves;
    PhyloNode *best_node;
    PhyloNeighbor *best_node_branch;

    PlacementCandidateNode() {}
};

#endif