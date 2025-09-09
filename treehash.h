//
// C++ Interface: treehash
//
// Description: Tree hashing utilities for fast topology comparison
//
//
// Author: Generated with Claude Code, (C) 2025
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TREEHASH_H
#define TREEHASH_H

#include <stdint.h>
#include <vector>

/**
 * 128-bit hash for tree topology comparison supporting 1B+ leaves
 * Using dual 64-bit values for collision resistance
 */
struct TreeHash128 {
    uint64_t high;
    uint64_t low;
    
    TreeHash128() : high(0), low(0) {}
    TreeHash128(uint64_t h, uint64_t l) : high(h), low(l) {}
    
    bool operator==(const TreeHash128& other) const {
        return high == other.high && low == other.low;
    }
    
    bool operator!=(const TreeHash128& other) const {
        return !(*this == other);
    }
    
    bool operator<(const TreeHash128& other) const {
        if (high != other.high) return high < other.high;
        return low < other.low;
    }
};

/**
 * Tree hashing utility functions (implemented in treehash.cpp)
 */

/**
 * Compute hash for a leaf node based on its index
 * @param leaf_index The index of the leaf node
 * @return 128-bit hash value for the leaf
 */
TreeHash128 computeLeafHash(int leaf_index);

/**
 * Combine two hash values using XOR and addition
 * @param a First hash value
 * @param b Second hash value  
 * @return Combined hash value
 */
TreeHash128 combineHashes(const TreeHash128& a, const TreeHash128& b);

/**
 * Compute hash for an internal node from its children hashes
 * @param child_hashes Vector of child hash values (should be sorted)
 * @return 128-bit hash value for the internal node
 */
TreeHash128 computeInternalNodeHash(const std::vector<TreeHash128>& child_hashes);

#endif // TREEHASH_H