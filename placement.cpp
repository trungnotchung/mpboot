#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "phylotree.h"
#include "alignment.h"
#include "iqtree.h"
#include "mutation.h"
#include "placement.h"

void initAlignment(IQTree *tree, Alignment *alignment, vector<int> &rotated_permutation_column) {
	int nsite = rotated_permutation_column.size();
	vector<int> perm_col(nsite);
	vector<int> compressed_perm_col(nsite);
	if (alignment->existing_sample_mutations.size()) {
		for (int site = 0; site < nsite; ++site) {
			int col = rotated_permutation_column[site];
			compressed_perm_col[site] = alignment->existing_sample_mutations[0][col].compressed_position;
			perm_col[site] = alignment->existing_sample_mutations[0][col].position;
		}
	}
	alignment->ungroupSitePattern();
	tree->add_row = true;
	tree->root_states = new UINT[(alignment->size() + 7) / 8 + 1];
	tree->initMutation(perm_col, compressed_perm_col);
}

int readInitialAlignment(ifstream &in_file_stream, char *out_file_name, int num_initial_rows) {
	ofstream out_file(out_file_name);
	if (!out_file.is_open()) {
		cout << "Cannot open outputfile :" << out_file_name << '\n';
		exit(1);
	}
	string line;
	int num_processed_rows = 0;
	while (getline(in_file_stream, line)) {
		if (line == "") {
			continue;
		}
		out_file << line << '\n';
		++num_processed_rows;
		if (num_processed_rows >= num_initial_rows) {
			break;
		}
	}
	out_file.close();
	return num_processed_rows;
}

int readVCFFile(IQTree *tree, Alignment*& alignment, Params &params) {
	if (params.num_existing_sequences + params.num_missing_sequences <= MAX_SEQUENCE) {
		alignment = new Alignment(params.aln_file, params.sequence_type, params.intype, params.num_existing_sequences);
		tree->setAlignment(alignment);
		tree->aln = alignment;
		vector<int> rotatedColumnPermutation = alignment->findRotatedColumnPermutation();
		initAlignment(tree, alignment, rotatedColumnPermutation);
		return (alignment)->getNSite();
	}

	ifstream in;
	in.exceptions(ios::failbit | ios::badbit);
	in.open(params.aln_file);
	string line;
	in.exceptions(ios::badbit);

	// Read first 12 lines and create tree alignment
	int totalColumn = readInitialAlignment(in, "temp.vcf", 12) - 1; // Read first 12 lines and write to temp.vcf
	alignment = new Alignment("temp.vcf", params.sequence_type, params.intype, params.num_existing_sequences);
	alignment->ungroupSitePattern();
	std::remove("temp.vcf");
	tree->setAlignment(alignment);
	tree->aln = alignment;

	vector<int> rotatedColumnPermutation = alignment->findRotatedColumnPermutation();
	initAlignment(tree, alignment, rotatedColumnPermutation);

	while (true) {
		int numProcessedColumn = (alignment)->readPartialVCF(in, params.sequence_type, rotatedColumnPermutation, params.num_existing_sequences, totalColumn, 8);
		if (numProcessedColumn == 0)
			break;
		tree->clearAllPartialLH();
		totalColumn += numProcessedColumn;
		initAlignment(tree, alignment, rotatedColumnPermutation);
	}

	in.close();
	return totalColumn;
}

void placeNewSamplesOntoExistingTree(Params &params) {
	cout << "\n========== Start initial data structure ==========\n";

	Alignment *alignment;
	IQTree *tree = new IQTree;
	bool is_rooted = false;

	tree->readTree(params.mutation_tree_file, is_rooted);
	int sequence_length = readVCFFile(tree, alignment, params) + 1;
	// Init new tree's memory
	tree->allocateMutationMemory(sequence_length);
	// free memory
	delete[] tree->root_states;
	tree->add_row = false;
	cout << "Tree parsimony after init mutations: " << tree->computeParsimonyScoreMutation() << '\n';

	cout << "\n========== Starting placement core ==========\n";
	int num_sequences = min((int)alignment->missing_sample_mutations.size(), params.num_missing_sequences);

	auto start_time = getCPUTime();
	for (int i = 0; i < num_sequences; ++i) {
		vector<pair<PhyloNode *, PhyloNeighbor *>> bfs = tree->breadth_first_expansion();
		int total_nodes = (int)bfs.size();

		CandidateNode inp;
		int best_set_difference = INT_MAX;
		size_t best_node_num_leaves = INT_MAX;
		size_t best_distance = INT_MAX;
		std::vector<Mutation> excess_mutations;
		std::vector<bool> node_has_unique(total_nodes, false);
		bool best_node_has_unique = false;
		size_t best_index = 0;

		inp.best_set_difference = &best_set_difference;
		inp.best_node_num_leaves = &best_node_num_leaves;
		inp.best_distance = &best_distance;
		inp.node = (PhyloNode *)tree->root->neighbors[0]->node;
		inp.node_branch = (PhyloNeighbor *)inp.node->findNeighbor(tree->root);
		inp.missing_sample_mutations = &alignment->missing_sample_mutations[i];
		inp.excess_mutations = &excess_mutations;
		inp.has_unique = &best_node_has_unique;
		inp.node_has_unique = &(node_has_unique);
		inp.best_index = &best_index;

		tree->initDataCalculatePlacementMutation(inp);
		tree->optimizedCalculatePlacementMutation(inp, 0, true);

		for (int j = 0; j < total_nodes; ++j) {
			if (inp.best_node == bfs[j].first) {
				best_index = j;
			}
		}
		*inp.best_set_difference = INT_MAX;
		inp.index = best_index;
		inp.node = bfs[best_index].first;
		inp.node_branch = bfs[best_index].second;
		tree->calculatePlacementMutation(inp, false, true);
		tree->addNewSample(bfs[best_index].first, bfs[best_index].second, excess_mutations, i, alignment->missing_seq_names[i]);
	}

	alignment->addToAlignmentNewSequences(alignment->missing_seq_names, alignment->missing_sequences);

	cout << "\n========== Finished placement core ==========\n";
	cout << "Time: " << fixed << setprecision(3) << (double)(getCPUTime() - start_time) << " seconds\n";
	cout << "Memory: " << getMemory() << " KB\n";
	
	cout << "New tree's parsimony score computed by mutation: " << tree->computeParsimonyScoreMutation() << '\n';
	tree->deleteAllPartialLh();
	cout << "New tree's parsimony score computed by fitch: " << tree->computeParsimony() << '\n';

	delete alignment;
	alignment = NULL;
	delete tree;
}

void checkCorectTree(char *origin_tree_file, char *new_tree_file) {
	cout << "================= Start checking correct tree ================\n";
	IQTree *origin_tree = new IQTree;
	bool origin_tree_is_rooted = false;
	origin_tree->readTree(origin_tree_file, origin_tree_is_rooted);

	IQTree *new_tree = new IQTree;
	bool new_tree_is_rooted = false;
	new_tree->readTree(new_tree_file, new_tree_is_rooted);

	vector<string> origin_tree_leaves_name;
	origin_tree->getLeafName(origin_tree_leaves_name);

	new_tree->assignRoot(origin_tree_leaves_name[0]);
	sort(origin_tree_leaves_name.begin(), origin_tree_leaves_name.end());
	new_tree->initInfoNode(origin_tree_leaves_name);

	if (new_tree->compareTree(origin_tree)) {
		cout << "Finish checking correct tree: Correct tree detected\n";
	}
	else {
		cout << "Finish checking correct tree: Wrong tree detected\n";
	}

	delete origin_tree;
	delete new_tree;
}

void configLeafNames(IQTree *tree, Node *node, Node *dad) {
	if (node->isLeaf()) {
		node->id = tree->aln->getSeqID(node->name);
	}
	FOR_NEIGHBOR_IT(node, dad, it)
	configLeafNames(tree, (*it)->node, node);
}