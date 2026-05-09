#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "phylotree.h"
#include "alignment.h"
#include "iqtree.h"
#include "mutation.h"
#include "placement.h"
#include "sproptimize.h"
extern int runSPRUnitTests(PhyloTree* tree);
extern int runSPRDeltaTests(PhyloTree* tree);
extern int runTBRUnitTests(PhyloTree* tree);
#include <queue>
#include <set>

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
	char* vcf_file = params.aln_file ? params.aln_file : params.user_file;

	if (params.num_existing_sequences + params.num_missing_sequences <= MAX_SEQUENCE) {
		alignment = new Alignment(vcf_file, params.sequence_type, params.intype, params.num_existing_sequences);
		tree->setAlignment(alignment);
		tree->aln = alignment;
		vector<int> rotatedColumnPermutation = alignment->findRotatedColumnPermutation();
		initAlignment(tree, alignment, rotatedColumnPermutation);
		return (alignment)->getNSite();
	}

	ifstream in;
	in.exceptions(ios::failbit | ios::badbit);
	in.open(vcf_file);
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

	tree->allocateMutationMemory(sequence_length);
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
	cout << "\n========== Finished placement core ==========\n";
	cout << "Time: " << fixed << setprecision(3) << (double)(getCPUTime() - start_time) << " seconds\n";
	cout << "Memory: " << getMemory() << " KB\n";
	cout << "New tree's parsimony score computed by mutation: " << tree->computeParsimonyScoreMutation() << '\n';

	alignment->addToAlignmentNewSequences(alignment->missing_seq_names, alignment->missing_sequences);
	tree->deleteAllPartialLh();

	int placement_score = tree->computeParsimony();
	cout << "Placement parsimony score (Fitch): " << placement_score << "\n";

	if (params.test_delta) {
		cout << "\n========== Running SPRDeltaExact Debug Tests ==========\n";
		tree->initializeAllPartialPars();
		int failures = runSPRDeltaTests(tree);
		if (failures > 0) {
			cerr << "\n[ERROR] " << failures << " delta test(s) FAILED" << endl;
		}
		cout << "========== SPRDeltaExact Debug Tests Complete ==========\n\n";

		cout << "\n========== Running TBR Unit Tests ==========\n";
		int tbr_failures = runTBRUnitTests(tree);
		if (tbr_failures > 0) {
			cerr << "\n[ERROR] " << tbr_failures << " TBR test(s) FAILED" << endl;
		}
		cout << "========== TBR Unit Tests Complete ==========\n\n";
	} else if (params.spr_test) {
		cout << "\n========== Running SPR Unit Tests ==========\n";
		tree->initializeAllPartialPars();
		int failures = runSPRUnitTests(tree);
		if (failures > 0) {
			cerr << "\n[ERROR] " << failures << " test(s) FAILED" << endl;
		}
		cout << "========== SPR Unit Tests Complete ==========\n\n";

		cout << "\n========== Running TBR Unit Tests ==========\n";
		int tbr_failures = runTBRUnitTests(tree);
		if (tbr_failures > 0) {
			cerr << "\n[ERROR] " << tbr_failures << " TBR test(s) FAILED" << endl;
		}
		cout << "========== TBR Unit Tests Complete ==========\n\n";
	} else if (params.spr_optimize) {
		cout << "\n========== Starting post-placement SPR optimization ==========\n";
		auto spr_start_time = getCPUTime();
		auto spr_wall_start = std::chrono::high_resolution_clock::now();

		SPROptimizer optimizer(tree);
		SPROptimizeOptions opts;
		opts.max_passes    = params.spr_max_passes;
		opts.max_radius    = params.spr_max_radius;
		opts.ratchet_iters = params.spr_ratchet_iterations;
		opts.ratchet_seed  = params.spr_ratchet_seed;
		opts.ratchet_runs  = params.spr_ratchet_runs;
		opts.wall_seconds  = params.spr_wall_seconds;
		opts.tbr_iters     = params.spr_tbr_iters;
		opts.tbr_radius_a  = params.spr_tbr_radius_a;
		opts.tbr_radius_b  = params.spr_tbr_radius_b;
		int best_score = optimizer.optimizeTree(opts);

		double wall_secs = std::chrono::duration<double>(
			std::chrono::high_resolution_clock::now() - spr_wall_start).count();
		cout << "SPR optimization time: " << fixed << setprecision(3)
		     << wall_secs << " seconds (wall), "
		     << (double)(getCPUTime() - spr_start_time) << " seconds (cpu)\n";
		cout << "Final parsimony score after SPR: " << best_score << '\n';
		tree->deleteAllPartialLh();
		cout << "Final parsimony score computed by fitch: " << tree->computeParsimony() << '\n';
		cout << "========== Finished SPR optimization ==========\n\n";
	}

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
	origin_tree->getLeavesName(origin_tree_leaves_name);

	new_tree->assignRoot(origin_tree_leaves_name[0]);
	sort(origin_tree_leaves_name.begin(), origin_tree_leaves_name.end());
	new_tree->initNodeData(origin_tree_leaves_name);

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