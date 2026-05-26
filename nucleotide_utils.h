#ifndef NUCLEOTIDE_UTILS_H
#define NUCLEOTIDE_UTILS_H

#include <cstdint>
#include <cassert>

// One-hot nucleotide encoding using bitwise operations
typedef uint8_t nuc_one_hot;

// Individual nucleotide constants
const nuc_one_hot NUC_A = 0b0001;  // 1
const nuc_one_hot NUC_C = 0b0010;  // 2
const nuc_one_hot NUC_G = 0b0100;  // 4
const nuc_one_hot NUC_T = 0b1000;  // 8
const nuc_one_hot NUC_N = 0b1111;  // 15 (ambiguous/any)

/**
 * Convert character nucleotide to one-hot encoding
 * @param c: Nucleotide character ('A', 'C', 'G', 'T', 'N', '-')
 * @return: One-hot encoded nucleotide
 */
inline nuc_one_hot char_to_one_hot(char c) {
    switch (c) {
        case 'A': case 'a': return NUC_A;
        case 'C': case 'c': return NUC_C;
        case 'G': case 'g': return NUC_G;
        case 'T': case 't': return NUC_T;

        // IUPAC ambiguity codes (proper two-nucleotide combinations)
        case 'R': case 'r': return NUC_A | NUC_G;  // puRine (A or G)
        case 'Y': case 'y': return NUC_C | NUC_T;  // pYrimidine (C or T)
        case 'K': case 'k': return NUC_G | NUC_T;  // Keto (G or T)
        case 'M': case 'm': return NUC_A | NUC_C;  // aMino (A or C)
        case 'S': case 's': return NUC_G | NUC_C;  // Strong (G or C)
        case 'W': case 'w': return NUC_A | NUC_T;  // Weak (A or T)

        // IUPAC three-nucleotide combinations
        case 'B': case 'b': return NUC_C | NUC_G | NUC_T;  // not A
        case 'D': case 'd': return NUC_A | NUC_G | NUC_T;  // not C
        case 'H': case 'h': return NUC_A | NUC_C | NUC_T;  // not G
        case 'V': case 'v': return NUC_A | NUC_C | NUC_G;  // not T

        // Any nucleotide or gap
        case 'N': case 'n': case '-': case '?': return NUC_N;

        default: return NUC_N;  // Treat unknown as ambiguous
    }
}

/**
 * Convert one-hot encoding to character (choose first set bit)
 * @param state: One-hot encoded nucleotide
 * @return: Character representation
 */
inline char one_hot_to_char(nuc_one_hot state) {
    if (state & NUC_A) return 'A';
    if (state & NUC_C) return 'C';
    if (state & NUC_G) return 'G';
    if (state & NUC_T) return 'T';
    return 'N';
}

/**
 * Count number of bits set (number of nucleotides in set)
 * @param state: One-hot encoded nucleotide set
 * @return: Number of nucleotides in the set
 */
inline int popcount(nuc_one_hot state) {
    return __builtin_popcount(state);
}

/**
 * Choose representative nucleotide from a set
 * Prefer lexicographic order: A > C > G > T
 * @param major_allele: One-hot encoded set of nucleotides
 * @return: Single nucleotide from the set
 */
inline nuc_one_hot choose_representative(nuc_one_hot major_allele) {
    if (major_allele & NUC_A) return NUC_A;
    if (major_allele & NUC_C) return NUC_C;
    if (major_allele & NUC_G) return NUC_G;
    if (major_allele & NUC_T) return NUC_T;
    return NUC_N;
}

/**
 * Check if a nucleotide is in a set
 * @param nuc: Single nucleotide (one-hot)
 * @param set: Set of nucleotides (one-hot)
 * @return: true if nuc is in set
 */
inline bool is_in_set(nuc_one_hot nuc, nuc_one_hot set) {
    return (nuc & set) != 0;
}

#endif // NUCLEOTIDE_UTILS_H
