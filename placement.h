#ifndef PLACEMENT_H
#define PLACEMENT_H

#include "tools.h"
#include "fstream"
#include "timeutil.h"

const int INF = (int)1e9 + 7;

int readFile(ifstream &inFileStream, char *outFileName, int numRow);

// read VCF file per 8 lines
int readVCFFile(IQTree *tree, Alignment **alignment, Params &params);

// add more K row using mutation like usher
void addMoreRowMutation(Params &params);

// check if origin tree doesn't change.
void checkCorrectTree(char *originTreeFile, char *newTreeFile);
#endif
