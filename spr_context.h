#ifndef SPR_CONTEXT_H
#define SPR_CONTEXT_H

class PhyloTree;

/**
 * Sets the active SPR tree.
 * @param tree The tree to register, or nullptr to clear.
 */
void setActiveSPRTree(PhyloTree* tree);

/**
 * @return The currently registered SPR tree, or nullptr.
 */
PhyloTree* getActiveSPRTree();

#endif // SPR_CONTEXT_H
