#include "spr_context.h"

static PhyloTree* s_active_spr_tree = nullptr;

void setActiveSPRTree(PhyloTree* tree) { s_active_spr_tree = tree; }
PhyloTree* getActiveSPRTree() { return s_active_spr_tree; }
