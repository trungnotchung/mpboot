#include "spr_context.h"

static PhyloTree* g_active_spr_tree = nullptr;

void setActiveSPRTree(PhyloTree* tree) { g_active_spr_tree = tree; }
PhyloTree* getActiveSPRTree() { return g_active_spr_tree; }
