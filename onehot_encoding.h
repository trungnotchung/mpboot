#ifndef ONEHOT_ENCODING_H
#define ONEHOT_ENCODING_H

#include <cstdint>
#include <string>

class OneHotEncoding {
public:
    // Nucleotide bit encodings
    static const uint8_t NUC_A   = 0b0001;  // 1
    static const uint8_t NUC_C   = 0b0010;  // 2
    static const uint8_t NUC_G   = 0b0100;  // 4
    static const uint8_t NUC_T   = 0b1000;  // 8
    static const uint8_t NUC_N   = 0b1111;  // 15 (all possible)
    static const uint8_t NUC_GAP = 0b0000;  // 0 (gap/missing)
    static const uint8_t NUC_INVALID = 0xFF;

    static uint8_t charToOneHot(char c);

    // Returns first nucleotide by priority A > C > G > T if ambiguous.
    static char oneHotToChar(uint8_t hot);

    // Returns string of all possible nucleotides (e.g., "AC" for A|C).
    static std::string oneHotToString(uint8_t hot);

    static inline uint8_t intersection(uint8_t a, uint8_t b) {
        return a & b;
    }

    static inline uint8_t union_op(uint8_t a, uint8_t b) {
        return a | b;
    }

    static inline bool agree(uint8_t a, uint8_t b) {
        return (a & b) != 0;
    }

    static inline int popcount(uint8_t hot) {
        return __builtin_popcount(hot);
    }

    static inline bool isUnambiguous(uint8_t hot) {
        return popcount(hot) == 1;
    }

    static inline bool isValid(uint8_t hot) {
        return hot != NUC_GAP && hot != NUC_INVALID;
    }
};

#endif // ONEHOT_ENCODING_H
