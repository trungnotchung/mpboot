#ifndef _MUTATION
#define _MUTATION 
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <cassert>
char get_nuc(int8_t nuc_id);

struct Mutation {
    int position;
    int compressed_position;
    char ref_nuc;
    char par_nuc;
    char mut_nuc;
    bool is_missing;

    Mutation() {
        is_missing = false;
    }

    inline bool operator < (const Mutation &m) const {
        return ((*this).position < m.position);
    }

    inline Mutation copy() const {
        Mutation m;
        m.position = position;
        m.ref_nuc = ref_nuc;
        m.par_nuc = par_nuc;
        m.mut_nuc = mut_nuc;
        m.is_missing = is_missing;
        m.compressed_position = compressed_position;
        return m;
    }

    inline bool is_masked() const {
        return (position < 0);
    }

    inline std::string get_string() const {
        if (is_masked()) {
            return "MASKED";
        }
        else {
            return get_nuc(par_nuc) + std::to_string(position) + get_nuc(mut_nuc);
        }
    }

    inline bool has_nuc(int nuc) {
        return ((1 << nuc) & mut_nuc) != 0;
    }
};
#endif