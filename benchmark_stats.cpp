#include "benchmark_stats.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <iomanip>

#ifdef __APPLE__
#include <mach/mach.h>
#elif __linux__
#include <sys/resource.h>
#endif

BenchmarkStats::BenchmarkStats()
    : placement_time(0.0)
    , spr_phase_a_time(0.0)
    , ratchet_phase_b_time(0.0)
    , tbr_phase_c_time(0.0)
    , total_time(0.0)
    , peak_memory_mb(0.0)
    , initial_parsimony(0)
    , final_parsimony(0)
{
}

double BenchmarkStats::getCurrentMemoryMB() {
#ifdef __APPLE__
    struct mach_task_basic_info info;
    mach_msg_type_number_t size = MACH_TASK_BASIC_INFO_COUNT;
    kern_return_t kerr = task_info(mach_task_self(),
                                   MACH_TASK_BASIC_INFO,
                                   (task_info_t)&info,
                                   &size);
    if (kerr == KERN_SUCCESS) {
        return info.resident_size / (1024.0 * 1024.0);
    }
    return 0.0;
#elif __linux__
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss / 1024.0;
#else
    return 0.0;
#endif
}

void BenchmarkStats::printSummary() const {
    std::cout << "\n=== MPBoot Placement Benchmark ===" << std::endl;
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Placement:        " << std::setw(8) << placement_time << "s" << std::endl;
    std::cout << "SPR Phase A:      " << std::setw(8) << spr_phase_a_time << "s" << std::endl;
    std::cout << "Ratchet Phase B:  " << std::setw(8) << ratchet_phase_b_time << "s" << std::endl;
    std::cout << "TBR Phase C:      " << std::setw(8) << tbr_phase_c_time << "s" << std::endl;
    std::cout << "Total:            " << std::setw(8) << total_time << "s" << std::endl;
    std::cout << "Peak Memory:      " << std::setw(8) << peak_memory_mb << " MB" << std::endl;
    std::cout << std::endl;
    std::cout << "Initial Parsimony: " << initial_parsimony << std::endl;
    std::cout << "Final Parsimony:   " << final_parsimony << std::endl;
    std::cout << "===================================" << std::endl;
}

void BenchmarkStats::writeToJSON(const std::string& output_path) const {
    std::time_t now = std::time(nullptr);
    char timestamp[32];
    std::strftime(timestamp, sizeof(timestamp), "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&now));

    std::ostringstream json;
    json << std::fixed << std::setprecision(2);
    json << "{\n";
    json << "  \"benchmark_version\": \"1.0\",\n";
    json << "  \"timestamp\": \"" << timestamp << "\",\n";
    json << "  \"timing\": {\n";
    json << "    \"placement_seconds\": " << placement_time << ",\n";
    json << "    \"spr_phase_a_seconds\": " << spr_phase_a_time << ",\n";
    json << "    \"ratchet_phase_b_seconds\": " << ratchet_phase_b_time << ",\n";
    json << "    \"tbr_phase_c_seconds\": " << tbr_phase_c_time << ",\n";
    json << "    \"total_seconds\": " << total_time << "\n";
    json << "  },\n";
    json << "  \"memory\": {\n";
    json << "    \"peak_memory_mb\": " << peak_memory_mb << "\n";
    json << "  },\n";
    json << "  \"parsimony\": {\n";
    json << "    \"initial_score\": " << initial_parsimony << ",\n";
    json << "    \"final_score\": " << final_parsimony << "\n";
    json << "  }\n";
    json << "}\n";

    std::ofstream out(output_path.c_str());
    if (out.is_open()) {
        out << json.str();
        out.close();
        std::cout << "\nBenchmark saved to: " << output_path << std::endl;
    } else {
        std::cerr << "Warning: Could not write benchmark file: " << output_path << std::endl;
    }
}
