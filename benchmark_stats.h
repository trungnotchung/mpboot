#ifndef BENCHMARK_STATS_H
#define BENCHMARK_STATS_H

#include <string>

/**
 * Tracks timing and memory metrics for placement pipeline benchmarking.
 * Designed for always-on benchmarking in placement mode (-pp_on).
 */
struct BenchmarkStats {
    double placement_time;
    double spr_phase_a_time;
    double ratchet_phase_b_time;
    double tbr_phase_c_time;
    double total_time;

    double peak_memory_mb;

    int initial_parsimony;
    int final_parsimony;

    BenchmarkStats();

    /**
     * Print human-readable summary to stdout.
     */
    void printSummary() const;

    /**
     * Write benchmark results to JSON file.
     * @param output_path Full path to output JSON file
     */
    void writeToJSON(const std::string& output_path) const;

    /**
     * Get current RSS memory usage in MB.
     * @return Memory usage in megabytes
     */
    static double getCurrentMemoryMB();
};

#endif // BENCHMARK_STATS_H
