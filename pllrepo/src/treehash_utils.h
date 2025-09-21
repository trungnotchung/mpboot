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
 * @file treehash_utils.h
 */

#ifndef __pll_TREEHASH_UTILS__
#define __pll_TREEHASH_UTILS__

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * 128-bit hash structure for tree topology comparison supporting 1B+ leaves
 * Using dual 64-bit values for collision resistance
 */
typedef struct {
    uint64_t high;
    uint64_t low;
} pllTreeHash128;

/**
 * Initialize a TreeHash128 to zero
 */
static inline pllTreeHash128 pllTreeHashInitZero(void) {
    pllTreeHash128 hash;
    hash.high = 0;
    hash.low = 0;
    return hash;
}

/**
 * Initialize a TreeHash128 with high and low values
 */
static inline pllTreeHash128 pllTreeHashInit(uint64_t high, uint64_t low) {
    pllTreeHash128 hash;
    hash.high = high;
    hash.low = low;
    return hash;
}

/**
 * Compare two TreeHash128 values for equality
 */
static inline int pllTreeHashEqual(const pllTreeHash128 *a, const pllTreeHash128 *b) {
    return (a->high == b->high) && (a->low == b->low);
}

/**
 * Compare two TreeHash128 values for ordering (less than)
 */
static inline int pllTreeHashLess(const pllTreeHash128 *a, const pllTreeHash128 *b) {
    if (a->high != b->high) return a->high < b->high;
    return a->low < b->low;
}

/**
 * Compute hash for a leaf node based on its index
 * @param leaf_index The index of the leaf node
 * @return 128-bit hash value for the leaf
 */
pllTreeHash128 pllTreeHashComputeLeaf(int leaf_index);

/**
 * Combine two hash values using XOR and addition
 * @param a First hash value
 * @param b Second hash value
 * @return Combined hash value
 */
pllTreeHash128 pllTreeHashCombine(const pllTreeHash128 *a, const pllTreeHash128 *b);

/**
 * Compute hash for an internal node from its children hashes
 * @param child_hashes Array of child hash values (should be sorted)
 * @param num_children Number of children
 * @return 128-bit hash value for the internal node
 */
pllTreeHash128 pllTreeHashComputeInternal(const pllTreeHash128 *child_hashes, int num_children);

/**
 * Sort an array of TreeHash128 values in place
 * @param hashes Array of hash values to sort
 * @param count Number of hash values in the array
 */
void pllTreeHashSort(pllTreeHash128 *hashes, int count);

/**
 * Comparison function for sorting TreeHash128 values
 * @param a First hash value
 * @param b Second hash value
 * @return -1 if a < b, 0 if a == b, 1 if a > b
 */
int pllTreeHashCompare(const void *a, const void *b);

#ifdef __cplusplus
}
#endif

#endif /* __pll_TREEHASH_UTILS__ */