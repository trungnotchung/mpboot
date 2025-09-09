//
// C++ Implementation: treehash
//
// Description: Tree hashing utilities for fast topology comparison
//
//
// Author: Generated with Claude Code, (C) 2025
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "treehash.h"

/**
 * Compute hash for a leaf node based on its index
 * Simple but effective hash function for leaf nodes
 */
TreeHash128 computeLeafHash(int leaf_index) {
    // Use a simple but effective hash function
    uint64_t hash = (uint64_t)leaf_index;
    hash = hash * 1000000007ULL + 123456789ULL;
    return TreeHash128(hash, hash * 31ULL);
}

/**
 * Combine two hash values using XOR and addition
 * This provides good mixing properties for hash combination
 */
TreeHash128 combineHashes(const TreeHash128& a, const TreeHash128& b) {
    return TreeHash128(a.high + b.high, a.low ^ b.low);
}

/**
 * Compute hash for an internal node from its children hashes
 * Combines all child hashes into a single hash value
 */
TreeHash128 computeInternalNodeHash(const std::vector<TreeHash128>& child_hashes) {
    if (child_hashes.empty()) {
        return TreeHash128(0, 0);
    }
    
    TreeHash128 result = child_hashes[0];
    for (size_t i = 1; i < child_hashes.size(); ++i) {
        result = combineHashes(result, child_hashes[i]);
    }
    
    return result;
}