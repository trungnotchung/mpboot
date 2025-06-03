#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "phylotree.h"
#include "alignment.h"
#include "iqtree.h"
#include "mutation.h"
#include "placement.h"

void checkCorectTree(char *originTreeFile, char *newTreeFile)
{
	cout << "================= Check correct tree ================\n";
	IQTree *originTree = new IQTree;
	bool originIsRooted = false;
	originTree->readTree(originTreeFile, originIsRooted);

	IQTree *newTree = new IQTree;
	bool newIsRooted = false;
	newTree->readTree(newTreeFile, newIsRooted);

	vector<string> originLeafName;
	originTree->getLeafName(originLeafName);

	newTree->assignRoot(originLeafName[0]);
	sort(originLeafName.begin(), originLeafName.end());
	newTree->initInfoNode(originLeafName);

	if (newTree->compareTree(originTree))
		cout << "Correct tree\n";
	else
		cout << "Wrong tree\n";

	delete originTree;
	delete newTree;
}

void configLeafNames(IQTree *tree, Node *node, Node *dad)
{
	if (node->isLeaf())
	{
		node->id = tree->aln->getSeqID(node->name);
	}
	FOR_NEIGHBOR_IT(node, dad, it)
	configLeafNames(tree, (*it)->node, node);
}

void initializeNewColumn(IQTree *tree, Alignment *alignment, vector<int> &rotatedPermutationColumn)
{
	int nsite = rotatedPermutationColumn.size();
	vector<int> permCol(nsite);
	vector<int> compressedPermCol(nsite);
	if (alignment->existingSampleMutations.size())
	{
		for (int site = 0; site < nsite; ++site)
		{
			int col = rotatedPermutationColumn[site];
			compressedPermCol[site] = alignment->existingSampleMutations[0][col].compressed_position;
			permCol[site] = alignment->existingSampleMutations[0][col].position;
		}
	}
	alignment->ungroupSitePattern();
	tree->add_row = true;
	tree->root_states = new UINT[(alignment->size() + 7) / 8 + 1];
	tree->initMutation(permCol, compressedPermCol);
}

int readInitialAlignment(ifstream &INT_MAXileStream, char *outFileName, int numInitialRow)
{
	ofstream outFile(outFileName);
	if (!outFile.is_open())
	{
		cout << "Cannot open outputfile :" << outFileName << '\n';
		exit(1);
	}
	string line;
	int currentRow = 0;
	while (getline(INT_MAXileStream, line))
	{
		if (line == "")
		{
			continue;
		}
		outFile << line << '\n';
		++currentRow;
		if (currentRow >= numInitialRow)
		{
			break;
		}
	}
	outFile.close();
	return currentRow;
}

int readVCFFile(IQTree *tree, Alignment **alignment, Params &params)
{
	char *alnFile = params.aln_file;
	ifstream in;
	in.exceptions(ios::failbit | ios::badbit);
	in.open(alnFile);
	string line;
	in.exceptions(ios::badbit);

	// Read first 12 lines and create tree alignment
	int totalColumn = readInitialAlignment(in, "temp.vcf", 12) - 1; // Read first 12 lines and write to temp.vcf
	*alignment = new Alignment("temp.vcf", params.sequence_type, params.intype, params.num_existing_sample);
	(*alignment)->ungroupSitePattern();
	std::remove("temp.vcf");
	tree->setAlignment(*alignment);
	tree->aln = *alignment;

	vector<int> rotatedColumnPermutation = (*alignment)->findRotatedColumnPermutation();
	initializeNewColumn(tree, *alignment, rotatedColumnPermutation);

	while (true)
	{
		int numProcessedColumn = (*alignment)->readPartialVCF(in, params.sequence_type, rotatedColumnPermutation, params.num_existing_sample, totalColumn, 8);
		if (numProcessedColumn == 0)
		{
			// Process all columns
			break;
		}
		tree->clearAllPartialLH();
		totalColumn += numProcessedColumn;
		initializeNewColumn(tree, *alignment, rotatedColumnPermutation);
	}

	in.close();
	return totalColumn;
}

void placeNewSamplesOntoExistingTree(Params &params)
{
	cout << "\n========== Start initial data structure ==========\n";

	Alignment *alignment;
	IQTree *tree;
	tree = new IQTree;
	char *fileName = params.mutation_tree_file;
	bool isRooted = false;

	tree->readTree(fileName, isRooted);
	int numColumn = readVCFFile(tree, &alignment, params) + 1;
	// Init new tree's memory
	tree->allocateMutationMemory(numColumn);
	// free memory
	delete[] tree->root_states;
	tree->add_row = false;
	cout << "Tree parsimony after init mutations: " << tree->computeParsimonyScoreMutation() << '\n';

	cout << "\n========== Starting placement core ==========\n";
	int numSample = min((int)alignment->missingSampleMutations.size(), params.num_missing_sample);

	auto startTime = getCPUTime();
	for (int i = 0; i < numSample; ++i)
	{
		vector<pair<PhyloNode *, PhyloNeighbor *>> bfs = tree->breadth_first_expansion();
		int totalNodes = (int)bfs.size();

		CandidateNode inp;
		int bestSetDifference = INT_MAX;
		size_t bestNodeNumLeaves = INT_MAX;
		size_t bestDistance = INT_MAX;
		std::vector<Mutation> excessMutations;
		std::vector<bool> nodeHasUnique(totalNodes, false);
		bool bestNodeHasUnique = false;
		size_t bestIndex = 0;

		inp.best_set_difference = &bestSetDifference;
		inp.best_node_num_leaves = &bestNodeNumLeaves;
		inp.best_distance = &bestDistance;
		inp.node = (PhyloNode *)tree->root->neighbors[0]->node;
		inp.node_branch = (PhyloNeighbor *)inp.node->findNeighbor(tree->root);
		inp.missing_sample_mutations = &alignment->missingSampleMutations[i];
		inp.excess_mutations = &excessMutations;
		inp.has_unique = &bestNodeHasUnique;
		inp.node_has_unique = &(nodeHasUnique);
		inp.best_index = &bestIndex;

		tree->initDataCalculatePlacementMutation(inp);
		tree->optimizedCalculatePlacementMutation(inp, 0, true);

		for (int j = 0; j < totalNodes; ++j)
		{
			if (inp.best_node == bfs[j].first)
			{
				bestIndex = j;
			}
		}
		*inp.best_set_difference = INT_MAX;
		inp.index = bestIndex;
		inp.node = bfs[bestIndex].first;
		inp.node_branch = bfs[bestIndex].second;
		tree->calculatePlacementMutation(inp, false, true);
		tree->addNewSample(bfs[bestIndex].first, bfs[bestIndex].second, excessMutations, i, alignment->missingSampleNames[i]);
	}

	cout << "\n========== Finished placement core ==========\n";
	cout << "New tree's parsimony score: " << tree->computeParsimonyScoreMutation() << '\n';
	cout << "Time: " << fixed << setprecision(3) << (double)(getCPUTime() - startTime) << " seconds\n";
	cout << "Memory: " << getMemory() << " KB\n";

	delete alignment;
	alignment = NULL;
	delete tree;
}