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

public:
    friend class TinaTree;
    friend class PhyloSuperTreePlen;
    friend class ParsTree; // DTH ##############

    /**
        construct class with a node and length		@param anode the other end of the branch

        @param alength length of branch
     */
    PhyloNeighbor(Node *anode, double alength) : Neighbor(anode, alength)
    {
        partial_lh = NULL;
        partial_lh_computed = 0;
        lh_scale_factor = 0.0;
        partial_pars = NULL;
        mutations.clear();
        canMove = 0;
    }

    /**
        construct class with a node and length
        @param anode the other end of the branch
        @param alength length of branch
        @param aid branch ID
     */
    PhyloNeighbor(Node *anode, double alength, int aid) : Neighbor(anode, alength, aid)
    {
        partial_lh = NULL;
        partial_lh_computed = 0;
        lh_scale_factor = 0.0;
        partial_pars = NULL;
        mutations.clear();
        canMove = 0;
    }

    /**
        tell that the partial likelihood vector is not computed
     */
    inline void clearPartialLh()
    {
        partial_lh_computed = 0;
    }

    /**
     *  tell that the partial likelihood vector is computed
     */
    inline void unclearPartialLh()
    {
        partial_lh_computed = 1;
    }

    /**
        clear all partial likelihood recursively in forward direction
        @param dad dad of this neighbor
     */
    void clearForwardPartialLh(Node *dad);

    std::vector<Mutation> mutations;

    std::vector<Mutation> saved_mutations;

    int num_leaves, distance;

    void clear_mutations();

    void add_mutation(Mutation mut);

    // mark this branch can be moved.
    void mark_can_move();

    // unmark this branch can be moved.
    void unmark_can_move();

    // check if this branch can be moved.
    bool check_can_move();

    // mark this branch can do spr.
    void mark_can_do_spr();

    // unmark this branch can do spr.
    void unmark_can_do_spr();

    // check if this branch can do spr.
    bool check_can_do_spr();

    // save this branch's mutations.
    void save_mutations();

    // restore this branch's mutations.
    void restore_mutation();

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

    /**
     * check if this branch can be movedor do SPR
     */
    int canMove;
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
};

/**
    Node vector
 */
typedef vector<PhyloNode *> PhyloNodeVector;

class CandidateNode
{
public:
    PhyloNode *node;
    PhyloNeighbor *node_branch;
    std::vector<Mutation> *missing_sample_mutations;

    int *best_set_difference;
    int *set_difference;
    size_t *best_node_num_leaves;
    size_t distance;
    size_t *best_distance;
    size_t index;
    size_t *best_index;
    size_t *num_best;
    PhyloNode *best_node;
    PhyloNeighbor *best_node_branch;

    std::vector<bool> *node_has_unique;
    std::vector<size_t> *best_j_vec;

    bool *has_unique;

    std::vector<Mutation> *excess_mutations;

    CandidateNode()
    {
    }
};

#endif