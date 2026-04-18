#include "onehot_encoding.h"
#include <cctype>

uint8_t OneHotEncoding::charToOneHot(char c) {
    c = toupper(c);

    switch (c) {
        case 'A': return NUC_A;
        case 'C': return NUC_C;
        case 'G': return NUC_G;
        case 'T':
        case 'U': return NUC_T;

        case 'N':
        case 'X': return NUC_N;

        case 'R': return NUC_A | NUC_G;
        case 'Y': return NUC_C | NUC_T;
        case 'S': return NUC_G | NUC_C;
        case 'W': return NUC_A | NUC_T;
        case 'K': return NUC_G | NUC_T;
        case 'M': return NUC_A | NUC_C;

        case 'B': return NUC_C | NUC_G | NUC_T;
        case 'D': return NUC_A | NUC_G | NUC_T;
        case 'H': return NUC_A | NUC_C | NUC_T;
        case 'V': return NUC_A | NUC_C | NUC_G;

        case '-':
        case '.':
        case '?': return NUC_GAP;

        default: return NUC_INVALID;
    }
}

char OneHotEncoding::oneHotToChar(uint8_t hot) {
    if (hot == NUC_A) return 'A';
    if (hot == NUC_C) return 'C';
    if (hot == NUC_G) return 'G';
    if (hot == NUC_T) return 'T';

    if (hot == NUC_N) return 'N';
    if (hot == NUC_GAP) return '-';
    if (hot == NUC_INVALID) return '?';

    if (hot & NUC_A) return 'A';
    if (hot & NUC_C) return 'C';
    if (hot & NUC_G) return 'G';
    if (hot & NUC_T) return 'T';

    return '?';
}

std::string OneHotEncoding::oneHotToString(uint8_t hot) {
    std::string result;

    if (hot == NUC_GAP) return "-";
    if (hot == NUC_INVALID) return "?";
    if (hot == NUC_N) return "N";

    if (hot & NUC_A) result += 'A';
    if (hot & NUC_C) result += 'C';
    if (hot & NUC_G) result += 'G';
    if (hot & NUC_T) result += 'T';

    if (result.empty()) return "?";
    return result;
}
