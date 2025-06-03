#ifndef PLACEMENT_H
#define PLACEMENT_H

#include "tools.h"
#include "fstream"
#include "timeutil.h"

/**
 * Place new samples onto existing tree
 */
void placeNewSamplesOntoExistingTree(Params &params);

/**
 * Check if origin tree doesn't change.
 */
void checkCorrectTree(char *originTreeFile, char *newTreeFile);
#endif
