#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "phylotree.h"
#include "alignment.h"
#include "iqtree.h"
#include "mutation.h"
#include "placement.h"

const int VCF_HEADER_LINES = 12;  // Number of header lines in VCF file
const int BATCH_SIZE = 8;         // Number of columns to process in each batch
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
	int totalColumn = readInitialAlignment(in, "temp.vcf", VCF_HEADER_LINES) - 1; // Read header lines and write to temp.vcf
	alignment = new Alignment("temp.vcf", params.sequence_type, params.intype, params.num_existing_sequences);
	alignment->ungroupSitePattern();
	std::remove("temp.vcf");
	tree->setAlignment(alignment);
	tree->aln = alignment;

	vector<int> rotatedColumnPermutation = alignment->findRotatedColumnPermutation();
	initAlignment(tree, alignment, rotatedColumnPermutation);

	while (true) {
		int numProcessedColumn = (alignment)->readPartialVCF(in, params.sequence_type, rotatedColumnPermutation, 
			params.num_existing_sequences, totalColumn, BATCH_SIZE);
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
		tree->initNodeDataPlaceNewSample();
		PlacementCandidateNode input;
		int best_set_difference = INT_MAX;
		size_t best_node_num_leaves = INT_MAX;
		std::vector<Mutation> excess_mutations;

		input.best_set_difference = &best_set_difference;
		input.best_node_num_leaves = &best_node_num_leaves;
		input.node = (PhyloNode *)tree->root->neighbors[0]->node;
		input.node_branch = (PhyloNeighbor *)input.node->findNeighbor(tree->root);
		input.missing_sample_mutations = &alignment->missing_sample_mutations[i];
		input.excess_mutations = &excess_mutations;

		tree->initNewSampleMutations(input);
		tree->optimizedFindPositionPlaceNewSample(input, 0);
		input.node = input.best_node;
		input.node_branch = input.best_node_branch;
		tree->computeExcessMutations(input);
		tree->addNewSample(input.best_node, input.best_node_branch, excess_mutations, i, alignment->missing_seq_names[i]);
	}
	cout << "Time: " << fixed << setprecision(3) << (double)(getCPUTime() - start_time) << " seconds\n";
	cout << "Memory: " << getMemory() << " KB\n";
	cout << "New tree's parsimony score computed by mutation: " << tree->computeParsimonyScoreMutation() << '\n';
	cout << "\n========== Finished placement core ==========\n";

	if (params.pp_spr) {
		tree->params = &params;
		tree->sprTransformationWithoutBreakingOriginalTree();

		if (params.pp_verify_preserved_tree) {
			ofstream fout("new_tree.treefile");
			tree->printTree(fout, WT_SORT_TAXA | WT_NEWLINE);
			fout.close();
			checkCorrectTree(params.mutation_tree_file, "new_tree.treefile", params.num_existing_sequences);
			std::remove("new_tree.treefile");
		}
	}

	delete alignment;
	alignment = NULL;
	delete tree;
}

void checkCorrectTree(char *origin_tree_file, char *new_tree_file, int n_original) {
	cout << "\n================= Start checking correct tree ================\n";
	IQTree *origin_tree = new IQTree;
	bool origin_tree_is_rooted = false;
	origin_tree->readTree(origin_tree_file, origin_tree_is_rooted);

	IQTree *new_tree = new IQTree;
	bool new_tree_is_rooted = false;
	new_tree->readTree(new_tree_file, new_tree_is_rooted);

	if (new_tree->compareTreeByHash(origin_tree, n_original)) {
		cout << "Correct tree detected\n";
	}
	else {
		cout << "Wrong tree detected\n";
	}
	cout << "\n================= Finished checking correct tree ================\n";

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