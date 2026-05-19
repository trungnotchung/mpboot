#ifndef PLACEMENT_H
#define PLACEMENT_H

#include "tools.h"
#include "fstream"
#include "timeutil.h"

const int MAX_SEQUENCE = 10000;

/** Place new VCF samples onto an existing reference tree. */
void placeNewSamplesOntoExistingTree(Params &params);

#endif
