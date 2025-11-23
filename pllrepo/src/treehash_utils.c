/**
 * PLL (version 1.0.0) a software library for phylogenetic inference
 * Copyright (C) 2013 Tomas Flouri and Alexandros Stamatakis
 *
 * Tree hashing utilities for fast topology comparison
 * Supporting 1B+ leaves with 128-bit hash values
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @file treehash_utils.c
 */

#include "treehash_utils.h"
#include <stdlib.h>

/**
 * Compute hash for a leaf node based on its index
 * Simple but effective hash function for leaf nodes
 */
pllTreeHash128 pllTreeHashComputeLeaf(int leaf_index) {
    /* Use a simple but effective hash function */
    uint64_t hash = (uint64_t)leaf_index;
    hash = hash * 1000000007ULL + 123456789ULL;

    pllTreeHash128 result;
    result.high = hash;
    result.low = hash * 31ULL;
    return result;
}

/**
 * Combine two hash values using XOR and addition
 * This provides good mixing properties for hash combination
 */
pllTreeHash128 pllTreeHashCombine(const pllTreeHash128 *a, const pllTreeHash128 *b) {
    pllTreeHash128 result;
    result.high = a->high + b->high;
    result.low = a->low ^ b->low;
    return result;
}

/**
 * Compute hash for an internal node from its children hashes
 * Combines all child hashes into a single hash value
 */
pllTreeHash128 pllTreeHashComputeInternal(const pllTreeHash128 *child_hashes, int num_children) {
    if (num_children == 0) {
        return pllTreeHashInitZero();
    }

    pllTreeHash128 result = child_hashes[0];
    int i;
    for (i = 1; i < num_children; ++i) {
        result = pllTreeHashCombine(&result, &child_hashes[i]);
    }

    return result;
}

/**
 * Comparison function for sorting TreeHash128 values
 * Used by qsort for ordering hash values
 */
int pllTreeHashCompare(const void *a, const void *b) {
    const pllTreeHash128 *hash_a = (const pllTreeHash128 *)a;
    const pllTreeHash128 *hash_b = (const pllTreeHash128 *)b;

    if (hash_a->high < hash_b->high) return -1;
    if (hash_a->high > hash_b->high) return 1;

    if (hash_a->low < hash_b->low) return -1;
    if (hash_a->low > hash_b->low) return 1;

    return 0;
}

/**
 * Sort an array of TreeHash128 values in place
 * Uses standard library qsort function
 */
void pllTreeHashSort(pllTreeHash128 *hashes, int count) {
    if (count > 1) {
        qsort(hashes, count, sizeof(pllTreeHash128), pllTreeHashCompare);
    }
}