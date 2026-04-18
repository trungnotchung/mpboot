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

typedef uint8_t nuc_one_hot;

char get_nuc(int8_t nuc_id);

struct Mutation {
    int position;
    int compressed_position;
    char ref_nuc;
    char par_nuc;
    char mut_nuc;
    bool is_missing;
    unsigned char boundary1_allele;

    nuc_one_hot par_one_hot;
    nuc_one_hot mut_one_hot;
    nuc_one_hot all_major_allele;

    Mutation() {
        is_missing = false;
        boundary1_allele = 0;
        par_one_hot = 0;
        mut_one_hot = 0;
        all_major_allele = 0;
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
        m.boundary1_allele = boundary1_allele;
        m.par_one_hot = par_one_hot;
        m.mut_one_hot = mut_one_hot;
        m.all_major_allele = all_major_allele;
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

    // UShER-compatible methods for exact SPR delta calculation
    inline bool is_valid() const {
        return par_one_hot != mut_one_hot;
    }

    inline nuc_one_hot get_par_one_hot() const {
        return par_one_hot;
    }

    inline nuc_one_hot get_mut_one_hot() const {
        return mut_one_hot;
    }

    inline nuc_one_hot get_boundary1_one_hot() const {
        return boundary1_allele;
    }
};
#endif