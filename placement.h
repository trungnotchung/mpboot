#ifndef PLACEMENT_H
#define PLACEMENT_H

#include "tools.h"
#include "fstream"
#include "timeutil.h"
#include <chrono>

const int MAX_SEQUENCE = 20000;

class Timer
{
   public:
    explicit Timer(std::string name)
        : name_(std::move(name)), start_(getCPUTime()) {}

    ~Timer()
    {
        std::cout << std::fixed << std::setprecision(3);
        std::cout << this->name_ << ": " << getCPUTime() - this->start_ << " seconds\n";
    }

   private:
    std::string name_;
    double start_;
};

class MemoryTracker
{
   public:
    MemoryTracker(std::string name)
        : name_(std::move(name)), start_(getMemory()) {}

    ~MemoryTracker()
    {
        std::cout << this->name_ << ": " << getMemory() - this->start_ << " KB\n";
    }

   private:
    std::string name_;
    uint64_t start_;
};

/**
 * Place new samples onto existing tree
 */
void placeNewSamplesOntoExistingTree(Params &params);

/**
 * Check if origin tree doesn't change.
 */
void checkCorrectTree(char *originTreeFile, char *newTreeFile);
#endif
