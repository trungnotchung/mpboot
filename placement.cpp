#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "phylotree.h"
#include "alignment.h"
#include "iqtree.h"
#include "mutation.h"
#include "placement.h"
#include "optimizer.h"
#include "benchmark_stats.h"
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

	if (params.pp_num_existing + params.pp_num_missing <= MAX_SEQUENCE) {
		alignment = new Alignment(vcf_file, params.sequence_type, params.intype, params.pp_num_existing);
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
	alignment = new Alignment("temp.vcf", params.sequence_type, params.intype, params.pp_num_existing);
	alignment->ungroupSitePattern();
	std::remove("temp.vcf");
	tree->setAlignment(alignment);
	tree->aln = alignment;

	vector<int> rotatedColumnPermutation = alignment->findRotatedColumnPermutation();
	initAlignment(tree, alignment, rotatedColumnPermutation);

	while (true) {
		int numProcessedColumn = (alignment)->readPartialVCF(in, params.sequence_type, rotatedColumnPermutation, 
			params.pp_num_existing, totalColumn, BATCH_SIZE);
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
	BenchmarkStats bench_stats;
	double pipeline_start = getRealTime();

	cout << "\n========== Start initial data structure ==========\n";

	Alignment *alignment;
	IQTree *tree = new IQTree;
	bool is_rooted = false;

	tree->readTree(params.pp_tree_file, is_rooted);

	int sequence_length = readVCFFile(tree, alignment, params) + 1;

	tree->allocateMutationMemory(sequence_length);
	delete[] tree->root_states;
	tree->add_row = false;

	cout << "Tree parsimony after init mutations: " << tree->computeParsimonyScoreMutation() << '\n';

	cout << "\n========== Starting placement core ==========\n";
	int num_sequences = min((int)alignment->missing_sample_mutations.size(), params.pp_num_missing);

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
	double placement_end = getRealTime();
	bench_stats.placement_time = placement_end - pipeline_start;
	cout << "Time: " << fixed << setprecision(3) << (double)(getCPUTime() - start_time) << " seconds\n";
	cout << "Memory: " << getMemory() << " KB\n";
	cout << "New tree's parsimony score computed by mutation: " << tree->computeParsimonyScoreMutation() << '\n';

	alignment->addToAlignmentNewSequences(alignment->missing_seq_names, alignment->missing_sequences);
	tree->deleteAllPartialLh();

	int placement_score = tree->computeParsimony();
	bench_stats.initial_parsimony = placement_score;
	cout << "Placement parsimony score (Fitch): " << placement_score << "\n";

	if (params.pp_optimize) {
		cout << "\n========== Starting post-placement optimization ==========\n";
		auto spr_start_time = getCPUTime();
		auto spr_wall_start = std::chrono::high_resolution_clock::now();

		PlacementOptimizer optimizer(tree);
		PlacementOptimizeOptions opts;
		opts.max_passes     = params.pp_max_passes;
		opts.max_radius     = params.pp_max_radius;
		opts.wall_seconds   = params.pp_wall_seconds;
		opts.ratchet_iters  = params.pp_ratchet_iters;
		opts.ratchet_seed   = params.pp_ratchet_seed;
		opts.ratchet_runs   = params.pp_ratchet_runs;
		opts.tbr_iters      = params.pp_tbr_iters;
		opts.tbr_max_radius = params.pp_tbr_max_radius;
		int best_score = optimizer.optimizeTree(opts, &bench_stats);

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

	std::string tree_file = params.out_prefix;
	tree_file += ".treefile";
	tree->printTree(tree_file.c_str(), WT_TAXON_ID | WT_SORT_TAXA);
	cout << "Final tree written to: " << tree_file << '\n';

	double pipeline_end = getRealTime();
	bench_stats.total_time = pipeline_end - pipeline_start;
	bench_stats.peak_memory_mb = BenchmarkStats::getCurrentMemoryMB();
	bench_stats.final_parsimony = tree->computeParsimony();

	std::string benchmark_file = params.out_prefix;
	benchmark_file += ".benchmark.json";

	bench_stats.printSummary();
	bench_stats.writeToJSON(benchmark_file);

	delete alignment;
	alignment = NULL;
	delete tree;
}

