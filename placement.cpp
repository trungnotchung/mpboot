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

string ppRunOriginalTbr(Alignment *alignment, Params &params, string newickTree = "")
{
	cout << "\n========== Start TBR core ==========\n";
	IQTree *tree = new IQTree(alignment);
	tree->params = &params;
	if (newickTree == "")
	{
		bool isRooted = params.is_rooted;
		tree->readTree(params.mutation_tree_file, isRooted);
	}
	else
	{
		tree->readTreeString(newickTree);
	}

	ofstream fout1("tree_tbr_1.txt");
	tree->drawTree(fout1, WT_SORT_TAXA | WT_NEWLINE);

	// becasue the tree is read from file and file contains its name, not its id
	configLeafNames(tree, tree->root, NULL);

	tree->initializeAllPartialPars();
	tree->clearAllPartialLH();
	tree->curScore = tree->computeParsimony();
	cout << "Tree's score before running tbr: " << tree->curScore << '\n';

	double start_time = getCPUTime();
	string newTreeString = tree->ppRunOriginalTbr();
	ofstream fout2("tree_tbr_2.txt");
	tree->drawTree(fout2, WT_SORT_TAXA | WT_NEWLINE);
	double end_time = getCPUTime();
	cout << "Time running TBR: " << fixed << setprecision(3) << (double)(end_time - start_time) << " seconds\n";
	cout << "Memory: " << getMemory() << " KB\n";

	if (params.pp_test_optimize)
	{
		ofstream fout("newTree.txt");
		tree->printTree(fout, WT_SORT_TAXA | WT_NEWLINE);
		fout.close();
		checkCorectTree(params.original_tree_file, "newTree.txt");
	}

	return newTreeString;
}

string ppRunOriginalSpr(Alignment *alignment, Params &params, string newickTree = "")
{
	cout << "\n========== Start spr core ==========\n";
	IQTree *tree = new IQTree(alignment);
	tree->params = &params;
	if (newickTree == "")
	{
		bool isRooted = params.is_rooted;
		tree->readTree(params.mutation_tree_file, isRooted);
	}
	else
	{
		tree->readTreeString(newickTree);
	}

	ofstream fout1("tree1.txt");
	tree->drawTree(fout1, WT_SORT_TAXA | WT_NEWLINE);
	configLeafNames(tree, tree->root, NULL);

	tree->initializeAllPartialPars();
	tree->clearAllPartialLH();
	tree->curScore = tree->computeParsimony();
	cout << "Tree's score before running spr: " << tree->curScore << '\n';

	double startTime = getCPUTime();
	string newTreeString = tree->ppRunOriginalSpr();
	ofstream fout2("tree2.txt");
	tree->drawTree(fout2, WT_SORT_TAXA | WT_NEWLINE);
	double endTime = getCPUTime();
	cout << "Time running SPR: " << fixed << setprecision(3) << (double)(endTime - startTime) << " seconds\n";
	cout << "Memory: " << getMemory() << " KB\n";

	if (params.pp_test_optimize)
	{
		ofstream fout("newTree.txt");
		tree->printTree(fout, WT_SORT_TAXA | WT_NEWLINE);
		fout.close();
		checkCorectTree(params.original_tree_file, "newTree.txt");
	}
	return newTreeString;
}

void initialize(IQTree *tree, Alignment *alignment, vector<int> &savePermCol, vector<int> &permCol, vector<int> &compressedPermCol)
{
	permCol.resize(savePermCol.size());
	compressedPermCol.resize(savePermCol.size());
	if (alignment->existingSamples.size())
	{
		for (int j = 0; j < savePermCol.size(); ++j)
		{
			int p = savePermCol[j];
			compressedPermCol[j] = alignment->existingSamples[0][p].compressed_position;
			permCol[j] = alignment->existingSamples[0][p].position;
		}
	}
	alignment->ungroupSitePattern();
	tree->add_row = true;
	tree->save_branch_states_dad = new UINT[(alignment->size() + 7) / 8 + 1];
	tree->computeParsimony();
	tree->initMutation(permCol, compressedPermCol);
}

int readFile(ifstream &inFileStream, char *outFileName, int numRow)
{
	ofstream outFile(outFileName);
	if (!outFile.is_open())
	{
		cout << "Cannot open file " << outFileName << '\n';
		return 0;
	}
	string line;
	int curRow = 0;
	while (getline(inFileStream, line))
	{
		if (line == "")
		{
			continue;
		}
		outFile << line << '\n';
		++curRow;
		if (curRow >= numRow)
		{
			break;
		}
	}
	outFile.close();
	return curRow;
}

int readVCFFile(IQTree *tree, Alignment **alignment, Params &params)
{
	char *alnFile = params.aln_file;
	ifstream in;
	in.exceptions(ios::failbit | ios::badbit);
	in.open(alnFile);
	string line;
	in.exceptions(ios::badbit);

	int totalColumn = readFile(in, "temp.vcf", 12) - 1;
	*alignment = new Alignment("temp.vcf", params.sequence_type, params.intype, params.numStartRow);
	(*alignment)->ungroupSitePattern();
	std::remove("temp.vcf");

	tree->setAlignment(*alignment);
	tree->aln = *alignment;

	vector<int> permCol = (*alignment)->findPermCol();
	vector<int> savePermCol = permCol;
	vector<int> compressedPermCol = permCol;
	initialize(tree, *alignment, savePermCol, permCol, compressedPermCol);

	auto startTime = getCPUTime();
	while (true)
	{
		startTime = getCPUTime();
		int numColumn = (*alignment)->readPartialVCF(in, params.sequence_type, savePermCol, params.numStartRow, totalColumn, 8);
		if (numColumn == 0)
		{
			break;
		}

		tree->clearAllPartialLH();
		totalColumn += numColumn;
		startTime = getCPUTime();
		initialize(tree, *alignment, savePermCol, permCol, compressedPermCol);
	}

	in.close();
	return totalColumn;
}

void addMoreRowMutation(Params &params)
{
	Alignment *alignment;

	IQTree *tree;
	tree = new IQTree;

	char *fileName = params.mutation_tree_file;
	bool isRooted = false;

	if (params.tree_zip_file != NULL)
	{
		tree->readTree(params.tree_zip_file, fileName, isRooted);
	}
	else
	{
		tree->readTree(fileName, isRooted);
	}

	int vecSize = readVCFFile(tree, &alignment, params) + 1;
	// Init new tree's memory
	tree->cur_missing_sample_mutations.resize(vecSize);
	tree->cur_ancestral_mutations.resize(vecSize);
	tree->visited_missing_sample_mutations.resize(vecSize);
	tree->visited_ancestral_mutations.resize(vecSize);
	tree->cur_excess_mutations.resize(vecSize);
	tree->visited_excess_mutations.resize(vecSize);

	cout << "\n========== Start placement core ==========\n";

	// free memory
	delete[] tree->save_branch_states_dad;
	tree->add_row = false;

	cout << "Tree parsimony after init mutations: " << tree->computeParsimonyScoreMutation() << '\n';
	int numSample = (int)alignment->missingSamples.size();
	vector<MutationNode> missingSamples(numSample);
	for (int i = 0; i < (int)alignment->missingSamples.size(); ++i)
	{
		missingSamples[i].mutations = alignment->missingSamples[i];
		missingSamples[i].name = alignment->remainName[i];
	}
	numSample = min(numSample, params.numAddRow);

	auto startTime = getCPUTime();

	for (int i = 0; i < numSample; ++i)
	{
		vector<pair<PhyloNode *, PhyloNeighbor *>> bfs = tree->breadth_first_expansion();

		int totalNodes = (int)bfs.size();
		CandidateNode inp;
		int bestSetDifference = INF;
		size_t bestNodeNumLeaves = INF;
		size_t bestDistance = INF;
		std::vector<Mutation> excessMutations;
		std::vector<bool> nodeHasUnique(totalNodes, false);
		bool bestNodeHasUnique = false;
		size_t bestJ = 0;

		inp.best_set_difference = &bestSetDifference;
		inp.best_node_num_leaves = &bestNodeNumLeaves;
		inp.best_distance = &bestDistance;
		inp.node = (PhyloNode *)tree->root->neighbors[0]->node;
		inp.node_branch = (PhyloNeighbor *)inp.node->findNeighbor(tree->root);
		inp.missing_sample_mutations = &missingSamples[i].mutations;
		inp.excess_mutations = &excessMutations;
		inp.has_unique = &bestNodeHasUnique;
		inp.node_has_unique = &(nodeHasUnique);
		inp.best_j = &bestJ;

		tree->initDataCalculatePlacementMutation(inp);
		tree->optimizedCalculatePlacementMutation(inp, 0, true);

		for (int j = 0; j < totalNodes; ++j)
		{
			if (inp.best_node == bfs[j].first)
			{
				bestJ = j;
				break;
			}
		}
		*inp.best_set_difference = INF;
		inp.j = bestJ;
		inp.node = bfs[bestJ].first;
		inp.node_branch = bfs[bestJ].second;
		tree->calculatePlacementMutation(inp, false, true);
		tree->addNewSample(bfs[bestJ].first, bfs[bestJ].second, excessMutations, i, missingSamples[i].name);
	}
	cout << "New tree's parsimony score: " << tree->computeParsimonyScoreMutation() << '\n';
	cout << "Time: " << fixed << setprecision(3) << (double)(getCPUTime() - startTime) << " seconds\n";
	cout << "Memory: " << getMemory() << " KB\n";

	// free memory
	tree->cur_missing_sample_mutations.clear();
	tree->cur_ancestral_mutations.clear();
	tree->visited_missing_sample_mutations.clear();
	tree->visited_ancestral_mutations.clear();
	tree->cur_excess_mutations.clear();
	tree->visited_excess_mutations.clear();

	delete alignment;
	alignment = NULL;
	delete tree;
}