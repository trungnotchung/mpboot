#include "spr_mutation_ops.h"
#include "fitch.h"
#include "phylonode.h"
#include <algorithm>

static void computeFitchChangeFromStored(
    uint8_t old_major, uint8_t old_bnd1,
    uint8_t dec_allele, uint8_t inc_allele,
    uint8_t& new_major, uint8_t& new_bnd1, int& children_delta) {

    new_major = old_major;
    new_bnd1 = old_bnd1;
    children_delta = 0;

    if (dec_allele) {
        uint8_t dec_in_major = dec_allele & old_major;
        if (dec_in_major) {
            uint8_t remaining_major = old_major & ~dec_allele;
            if (remaining_major == 0) {
                uint8_t remaining_bnd1 = old_bnd1 & ~dec_allele;
                new_major = old_major | remaining_bnd1;
                new_bnd1 = 0; // don't know level M-2
                children_delta = 1;
            } else {
                new_major = remaining_major;
                new_bnd1 = dec_in_major | (old_bnd1 & ~dec_allele);
                children_delta = 0;
            }
        } else if (dec_allele & old_bnd1) {
            new_bnd1 = old_bnd1 & ~dec_allele;
        }
    }

    if (inc_allele) {
        uint8_t inc_in_new_major = inc_allele & new_major;
        uint8_t inc_in_new_bnd1 = inc_allele & new_bnd1;

        if (inc_in_new_bnd1 && !(inc_in_new_major)) {
            new_major = new_major | inc_in_new_bnd1;
            new_bnd1 = new_bnd1 & ~inc_in_new_bnd1;
        } else if (!inc_in_new_major && !inc_in_new_bnd1) {
            // Approximate: actual level depends on count levels
        } else if (inc_in_new_major) {
            new_bnd1 = (new_major & ~inc_allele) | new_bnd1;
            new_major = inc_in_new_major;
            children_delta -= 1; // max increased by 1
        }
    }
}


std::vector<Mutation>* SPRMutationOps::getMutations(PhyloNode* node, PhyloNode* dad) {
    if (!node) {
        return nullptr;
    }

    if (!dad) {
        if (node->neighbors.empty()) {
            return nullptr;
        }
        PhyloNeighbor* nei = (PhyloNeighbor*)node->neighbors[0];
        return &nei->mutations;
    }

    for (auto it = node->neighbors.begin(); it != node->neighbors.end(); ++it) {
        PhyloNeighbor* nei = (PhyloNeighbor*)*it;
        if (nei->node == dad) {
            return &nei->mutations;
        }
    }

    return nullptr;
}

static PhyloNode* s_tree_root = nullptr;
static Fitch* s_custom_fitch = nullptr;

void SPRMutationOps::setRoot(PhyloNode* root) {
    s_tree_root = root;
}

void SPRMutationOps::setCustomFitch(Fitch* fitch) {
    s_custom_fitch = fitch;
}

PhyloNode* SPRMutationOps::getParent(PhyloNode* node) {
    if (!node || node->neighbors.empty()) {
        return nullptr;
    }
    if (node == s_tree_root) {
        return nullptr;
    }
    // After orientTreeToRoot(), neighbors[0] points toward parent
    return (PhyloNode*)node->neighbors[0]->node;
}


MutationCountChangeCollection SPRMutationOps::initMutationChange(PhyloNode* src, PhyloNode* src_parent) {
    MutationCountChangeCollection mutations;

    if (!src || !src_parent) {
        return mutations;
    }

    std::vector<Mutation>* src_muts = getMutations(src, src_parent);
    if (!src_muts || src_muts->empty()) {
        return mutations;
    }

    for (const Mutation& m : *src_muts) {
        if (m.is_valid() || m.major_allele_set != m.mut_one_hot) {
            mutations.emplace_back(
                m.position,
                0,  // decremented
                m.major_allele_set  // incremented
            );

            mutations.back().par_state = m.par_one_hot;
            mutations.back().major_allele_set = m.major_allele_set;
            mutations.back().boundary1_allele = m.boundary1_allele;
            mutations.back().from_src = true;
        }
    }

    return mutations;
}

void SPRMutationOps::mergeMutationSrcToLCA(PhyloNode* ancestor,
                                           MutationCountChangeCollection& mutations) {
    if (!ancestor) {
        return;
    }

    std::vector<Mutation>* ancestor_muts = getMutations(ancestor, nullptr);
    if (!ancestor_muts) {
        return;
    }

    MutationCountChangeCollection merged_mutations;
    merged_mutations.reserve(mutations.size() + ancestor_muts->size());

    auto iter = mutations.begin();
    auto end = mutations.end();

    for (const auto& m : *ancestor_muts) {
        if (!m.is_valid()) {
            continue;
        }

        while (iter != end && iter->position < m.position) {
            merged_mutations.push_back(*iter);
            iter++;
        }

        if (iter != end && iter->position == m.position) {
            nuc_one_hot new_par_nuc = m.get_par_one_hot();

            if (new_par_nuc != iter->added_alleles) {
                merged_mutations.push_back(*iter);
                merged_mutations.back().par_state = m.get_par_one_hot();
            }

            iter++;
        } else {
            merged_mutations.emplace_back(m.position, 0, m.get_mut_one_hot());
            merged_mutations.back().par_state = m.get_par_one_hot();
            merged_mutations.back().major_allele_set = m.major_allele_set;
            merged_mutations.back().boundary1_allele = m.boundary1_allele;
        }
    }

    while (iter != end) {
        merged_mutations.push_back(*iter);
        iter++;
    }

    mutations = std::move(merged_mutations);
}

MutationCountChangeCollection SPRMutationOps::mergeMutationLCAToRank(
    PhyloNode* child_on_path,
    const MutationCountChangeCollection& mutations) {
    MutationCountChangeCollection child_mutations;

    if (!child_on_path) {
        return mutations;
    }

    std::vector<Mutation>* child_edge_muts = getMutations(child_on_path, nullptr);
    if (!child_edge_muts) {
        return mutations;
    }

    child_mutations.reserve(mutations.size() + child_edge_muts->size());

    auto child_mutation_iter = child_edge_muts->begin();
    auto child_mutation_end = child_edge_muts->end();

    for (const auto& m : mutations) {
        while (child_mutation_iter != child_mutation_end &&
               !child_mutation_iter->is_valid()) {
            child_mutation_iter++;
        }

        while (child_mutation_iter != child_mutation_end &&
               child_mutation_iter->position < m.position) {

            child_mutations.emplace_back(
                child_mutation_iter->position,
                0,
                child_mutation_iter->get_par_one_hot()
            );

            child_mutations.back().par_state = child_mutation_iter->get_mut_one_hot();
            child_mutations.back().major_allele_set = child_mutation_iter->major_allele_set;
            child_mutations.back().boundary1_allele = child_mutation_iter->boundary1_allele;

            child_mutation_iter++;
        }

        if (child_mutation_iter != child_mutation_end &&
            child_mutation_iter->position == m.position) {

            nuc_one_hot new_par_allele = child_mutation_iter->get_mut_one_hot();

            if (m.added_alleles != new_par_allele) {
                child_mutations.push_back(m);
                child_mutations.back().par_state = new_par_allele;
            }

            child_mutation_iter++;
        } else {
            child_mutations.push_back(m);
        }
    }

    while (child_mutation_iter != child_mutation_end) {
        child_mutations.emplace_back(
            child_mutation_iter->position,
            0,  // decremented
            child_mutation_iter->get_par_one_hot()  // incremented
        );

        child_mutations.back().par_state = child_mutation_iter->get_mut_one_hot();
        child_mutations.back().major_allele_set = child_mutation_iter->major_allele_set;
        child_mutations.back().boundary1_allele = child_mutation_iter->boundary1_allele;

        child_mutation_iter++;
    }

    return child_mutations;
}

void SPRMutationOps::getParentAlteredRemove(MutationCountChangeCollection& output,
                                            PhyloNode* src,
                                            int& parsimony_score_change) {
    if (!src) {
        return;
    }

    PhyloNode* parent = getParent(src);
    if (!parent) {
        return;
    }

    std::vector<Mutation>* parent_mutations = getMutations(parent, nullptr);
    if (!parent_mutations) {
        return;
    }

    std::vector<Mutation>* src_mutations = getMutations(src, parent);
    if (!src_mutations) {
        return;
    }

    bool is_binary = false;
    PhyloNode* sibling = nullptr;
    if (s_custom_fitch) {
        int num_children = 0;
        PhyloNode* gp = getParent(parent);
        FOR_NEIGHBOR_IT(parent, gp, nit) {
            PhyloNode* child = (PhyloNode*)(*nit)->node;
            num_children++;
            if (child != src) sibling = child;
        }
        is_binary = (num_children == 2 && sibling != nullptr);
    }

    for (const Mutation& parent_mut : *parent_mutations) {
        if (!parent_mut.is_valid()) {
            continue;
        }

        int position = parent_mut.position;
        uint8_t old_major = parent_mut.major_allele_set;
        uint8_t old_bnd1 = parent_mut.boundary1_allele;

        uint8_t src_allele = getSrcAlleleAtPosition(src, position);
        uint8_t effective_src_allele = src_allele;
        if (effective_src_allele == 0) {
            // When src has no edge mutation, its assigned state equals parent's.
            // But src's Fitch set may be broader (e.g., {C,T} vs parent's {C}).
            // For binary nodes, use src's actual Fitch set from Fitch;
            // for non-binary, old_major is a reasonable approximation.
            if (is_binary && s_custom_fitch) {
                int ptn = parent_mut.compressed_position;
                effective_src_allele = s_custom_fitch->getMajorForNode(src, ptn);
            } else {
                effective_src_allele = old_major;
            }
        }

        uint8_t computed_major, computed_bnd1;
        int children_delta;
        computeFitchChangeFromStored(old_major, old_bnd1,
                                      effective_src_allele, 0,
                                      computed_major, computed_bnd1, children_delta);

        parsimony_score_change += children_delta - 1;

        uint8_t new_major, new_bnd1;
        if (is_binary) {
            int ptn = parent_mut.compressed_position;
            new_major = s_custom_fitch->getMajorForNode(sibling, ptn);
            new_bnd1 = 0;  // only one child remains, no boundary
        } else {
            new_major = computed_major;
            new_bnd1 = computed_bnd1;
        }

        if (new_major != old_major) {
            output.emplace_back(
                position,
                old_major & ~new_major,
                new_major & ~old_major
            );
            output.back().par_state = parent_mut.get_par_one_hot();
            output.back().major_allele_set = new_major;
            output.back().boundary1_allele = new_bnd1;
        }
    }

    for (const Mutation& src_mut : *src_mutations) {
        if (!src_mut.is_valid()) {
            continue;
        }

        bool parent_has_mutation = false;
        for (const Mutation& pm : *parent_mutations) {
            if (pm.position == src_mut.position) {
                parent_has_mutation = true;
                break;
            }
        }

        if (!parent_has_mutation) {
            parsimony_score_change--;

            // For binary nodes: removing src also changes parent's Fitch set
            // at positions where src has a mutation but parent doesn't.
            // Parent's old major at this position includes src's contribution;
            // after removal, it becomes just sibling's Fitch set.
            // We must emit propagation entries so intermediate/dst-side see them.
            if (is_binary && s_custom_fitch) {
                int ptn = src_mut.compressed_position;
                uint8_t old_par_major = s_custom_fitch->getMajorForNode(parent, ptn);
                uint8_t new_par_major = s_custom_fitch->getMajorForNode(sibling, ptn);
                if (new_par_major != old_par_major) {
                    output.emplace_back(
                        src_mut.position,
                        old_par_major & ~new_par_major,
                        new_par_major & ~old_par_major
                    );
                    output.back().par_state = src_mut.get_par_one_hot();
                    output.back().major_allele_set = new_par_major;
                    output.back().boundary1_allele = 0;
                }
            }
        }
    }

}

void SPRMutationOps::getIntermediateNodesMutations(
    PhyloNode* node,
    const MutationCountChangeCollection& child_changes,
    MutationCountChangeCollection& parent_changes,
    int& score_change,
    PhyloNode* child_on_path) {
    if (!node) {
        parent_changes = child_changes;
        return;
    }

    parent_changes.clear();
    parent_changes.reserve(child_changes.size());

    PhyloNode* parent = getParent(node);
    std::vector<Mutation>* node_mutations = getMutations(node, parent);

    if (!node_mutations) {
        parent_changes = child_changes;
        return;
    }

    bool is_binary = false;
    PhyloNode* sibling = nullptr;
    if (child_on_path && s_custom_fitch) {
        int num_children = 0;
        FOR_NEIGHBOR_IT(node, parent, nit) {
            PhyloNode* child = (PhyloNode*)(*nit)->node;
            num_children++;
            if (child != child_on_path) sibling = child;
        }
        is_binary = (num_children == 2 && sibling != nullptr);
    }

    for (const MutationCountChange& change : child_changes) {
        const Mutation* node_mut = nullptr;
        for (const Mutation& m : *node_mutations) {
            if (m.position == change.position) {
                node_mut = &m;
                break;
            }
        }

        if (node_mut && node_mut->is_valid()) {
            uint8_t all_major = node_mut->major_allele_set;
            uint8_t bnd1 = node_mut->boundary1_allele;
            uint8_t dec = change.removed_alleles;
            uint8_t inc = change.added_alleles;
            uint8_t par = node_mut->get_par_one_hot();

            size_t entries_before = parent_changes.size();

            if (inc && !dec) {
                uint8_t inc_in_major = all_major & inc;
                if (inc_in_major) {
                    score_change--;
                    uint8_t not_incremented = all_major & ~inc;
                    if (not_incremented) {
                        parent_changes.emplace_back(
                            change.position, not_incremented, 0);
                        parent_changes.back().par_state = par;
                    }
                } else if (bnd1 & inc) {
                    parent_changes.emplace_back(
                        change.position, 0, bnd1 & inc);
                    parent_changes.back().par_state = par;
                }
            } else if (dec && !inc) {
                uint8_t dec_in_major = all_major & dec;
                if (dec_in_major) {
                    uint8_t not_decremented = all_major & ~dec;
                    if (not_decremented) {
                        parent_changes.emplace_back(
                            change.position, dec_in_major, 0);
                        parent_changes.back().par_state = par;
                    } else {
                        score_change++;
                        uint8_t bnd1_not_dec = bnd1 & ~dec;
                        if (bnd1_not_dec) {
                            parent_changes.emplace_back(
                                change.position, 0, bnd1_not_dec);
                            parent_changes.back().par_state = par;
                        }
                    }
                }
            } else if (dec && inc) {
                uint8_t major_inc = all_major & inc;
                uint8_t major_dec = all_major & dec;
                uint8_t major_not_dec = all_major & ~dec;

                if (major_inc) {
                    score_change--;
                    uint8_t not_incremented = all_major & ~inc;
                    if (not_incremented) {
                        parent_changes.emplace_back(
                            change.position, not_incremented, 0);
                        parent_changes.back().par_state = par;
                    }
                } else if (major_dec) {
                    if ((bnd1 & inc) || major_not_dec) {
                        parent_changes.emplace_back(
                            change.position, major_dec, bnd1 & inc);
                        parent_changes.back().par_state = par;
                    } else {
                        uint8_t not_decremented = all_major & ~dec;
                        if (not_decremented) {
                            parent_changes.emplace_back(
                                change.position, major_dec, 0);
                            parent_changes.back().par_state = par;
                        } else {
                            score_change++;
                            uint8_t bnd1_not_dec = bnd1 & ~dec;
                            if (bnd1_not_dec) {
                                parent_changes.emplace_back(
                                    change.position, 0, bnd1_not_dec);
                                parent_changes.back().par_state = par;
                            }
                        }
                    }
                } else {
                    if (inc) {
                        uint8_t inc_in_major2 = all_major & inc;
                        if (inc_in_major2) {
                            score_change--;
                            uint8_t not_inc = all_major & ~inc;
                            if (not_inc) {
                                parent_changes.emplace_back(
                                    change.position, not_inc, 0);
                                parent_changes.back().par_state = par;
                            }
                        } else if (bnd1 & inc) {
                            parent_changes.emplace_back(
                                change.position, 0, bnd1 & inc);
                            parent_changes.back().par_state = par;
                        }
                    }
                }
            }

            // Binary node fix: score_change is correct from above, but propagation
            // entries may use bnd1 which conflates both children for binary nodes.
            // Recompute new_N_major from CF and replace propagation entries.
            if (is_binary) {
                parent_changes.resize(entries_before);

                int ptn = node_mut->compressed_position;
                uint8_t A_old = s_custom_fitch->getMajorForNode(child_on_path, ptn);
                uint8_t B = s_custom_fitch->getMajorForNode(sibling, ptn);
                uint8_t A_new = (A_old & ~dec) | inc;
                uint8_t new_intersect = A_new & B;
                uint8_t new_N_major = new_intersect ? new_intersect : (A_new | B);

                if (new_N_major != all_major) {
                    parent_changes.emplace_back(
                        change.position,
                        all_major & ~new_N_major,
                        new_N_major & ~all_major);
                    parent_changes.back().par_state = par;
                    parent_changes.back().major_allele_set = new_N_major;
                    parent_changes.back().boundary1_allele = 0;
                }
            }
        } else {
            score_change += change.scoreDeltaInternal();
        }
    }
}

uint8_t SPRMutationOps::getSrcAlleleAtPosition(PhyloNode* src, int position) {
    if (!src) {
        return 0;
    }

    PhyloNode* parent = getParent(src);
    if (!parent) {
        return 0;
    }

    std::vector<Mutation>* src_mutations = getMutations(src, parent);
    if (!src_mutations) {
        return 0;
    }

    for (const Mutation& m : *src_mutations) {
        if (m.position == position) {
            return m.major_allele_set;
        }
        if (m.position > position) {
            break;
        }
    }

    return 0;
}

uint8_t SPRMutationOps::getSiblingAllele(PhyloNode* parent, PhyloNode* src, int position) {
    if (!parent) {
        return 0;
    }

    PhyloNode* sibling = nullptr;
    PhyloNode* grandparent = getParent(parent);

    for (auto it = parent->neighbors.begin(); it != parent->neighbors.end(); ++it) {
        PhyloNode* neighbor = (PhyloNode*)((*it)->node);
        if (neighbor != src && neighbor != grandparent) {
            sibling = neighbor;
            break;
        }
    }

    if (!sibling) {
        return 0;
    }

    return getSrcAlleleAtPosition(sibling, position);
}

int SPRMutationOps::recomputeMajorAllele(PhyloNode* node, int position,
                                         uint8_t dec_allele, uint8_t inc_allele,
                                         uint8_t& new_major, uint8_t& new_boundary1,
                                         int* children_node_delta) {
    if (!node) {
        new_major = 0;
        new_boundary1 = 0;
        return 0;
    }

    int allele_counts[4] = {0, 0, 0, 0};
    PhyloNode* parent_of_node = getParent(node);

    // Used as default allele when child has no mutation at this position
    uint8_t node_major_at_pos = 0;
    {
        std::vector<Mutation>* node_muts = getMutations(node, parent_of_node);
        if (node_muts) {
            for (const Mutation& m : *node_muts) {
                if (m.position == position) {
                    node_major_at_pos = m.major_allele_set;
                    break;
                }
                if (m.position > position) break;
            }
        }
    }

    // 1: Count alleles from all children WITHOUT modification
    for (auto it = node->neighbors.begin(); it != node->neighbors.end(); ++it) {
        PhyloNode* child = (PhyloNode*)((*it)->node);

        if (child == parent_of_node) {
            continue;
        }

        uint8_t child_allele = getSrcAlleleAtPosition(child, position);
        if (child_allele == 0 && node_major_at_pos != 0) {
            child_allele = node_major_at_pos;
        }

        if (child_allele & 0x1) allele_counts[0]++;
        if (child_allele & 0x2) allele_counts[1]++;
        if (child_allele & 0x4) allele_counts[2]++;
        if (child_allele & 0x8) allele_counts[3]++;
    }

    int old_max_count = 0;
    for (int i = 0; i < 4; i++) {
        if (allele_counts[i] > old_max_count) old_max_count = allele_counts[i];
    }
    uint8_t old_major = 0;
    if (allele_counts[0] == old_max_count) old_major |= 0x1;
    if (allele_counts[1] == old_max_count) old_major |= 0x2;
    if (allele_counts[2] == old_max_count) old_major |= 0x4;
    if (allele_counts[3] == old_max_count) old_major |= 0x8;

    // 2: Apply dec/inc to counts
    if (dec_allele & 0x1) allele_counts[0]--;
    if (dec_allele & 0x2) allele_counts[1]--;
    if (dec_allele & 0x4) allele_counts[2]--;
    if (dec_allele & 0x8) allele_counts[3]--;
    if (inc_allele & 0x1) allele_counts[0]++;
    if (inc_allele & 0x2) allele_counts[1]++;
    if (inc_allele & 0x4) allele_counts[2]++;
    if (inc_allele & 0x8) allele_counts[3]++;

    int max_count = 0;
    for (int i = 0; i < 4; i++) {
        if (allele_counts[i] > max_count) {
            max_count = allele_counts[i];
        }
    }

    new_major = 0;
    if (allele_counts[0] == max_count) new_major |= 0x1;
    if (allele_counts[1] == max_count) new_major |= 0x2;
    if (allele_counts[2] == max_count) new_major |= 0x4;
    if (allele_counts[3] == max_count) new_major |= 0x8;

    int boundary1_count = max_count - 1;
    new_boundary1 = 0;
    if (boundary1_count > 0) {
        if (allele_counts[0] == boundary1_count) new_boundary1 |= 0x1;
        if (allele_counts[1] == boundary1_count) new_boundary1 |= 0x2;
        if (allele_counts[2] == boundary1_count) new_boundary1 |= 0x4;
        if (allele_counts[3] == boundary1_count) new_boundary1 |= 0x8;
    }

    // old_max - new_max, WITHOUT -1 for child removal
    if (children_node_delta) {
        *children_node_delta = old_max_count - max_count;
    }

    PhyloNode* parent = getParent(node);
    if (!parent) {
        return 0;
    }

    std::vector<Mutation>* node_mutations = getMutations(node, parent);
    uint8_t par_state = 0;
    if (node_mutations) {
        for (const Mutation& m : *node_mutations) {
            if (m.position == position) {
                par_state = m.get_par_one_hot();
                break;
            }
        }
    }

    // 3: Node-to-parent edge change only
    int score_change = 0;
    if (par_state) {
        bool old_follows_parent = (old_major & par_state) != 0;
        bool new_follows_parent = (new_major & par_state) != 0;

        if (!old_follows_parent && new_follows_parent) {
            score_change = -1;
        } else if (old_follows_parent && !new_follows_parent) {
            score_change = +1;
        }
    }

    return score_change;
}

int SPRMutationOps::recomputeFullDelta(PhyloNode* node, int position,
                                        uint8_t dec_allele, uint8_t inc_allele,
                                        uint8_t change_par_state) {
    if (!node) return 0;

    PhyloNode* parent_of_node = getParent(node);
    std::vector<Mutation>* node_muts = getMutations(node, parent_of_node);
    const Mutation* node_mut = nullptr;

    if (node_muts) {
        for (const Mutation& m : *node_muts) {
            if (m.position == position) {
                node_mut = &m;
                break;
            }
            if (m.position > position) break;
        }
    }

    if (!node_mut || !node_mut->is_valid()) {
        return 0;
    }

    uint8_t all_major = node_mut->major_allele_set;
    uint8_t bnd1 = node_mut->boundary1_allele;
    uint8_t par_state = node_mut->get_par_one_hot();
    int score_change = 0;

    if (inc_allele && !dec_allele) {
        uint8_t inc_in_major = all_major & inc_allele;
        if (inc_in_major) {
            score_change--;
        }
    } else if (dec_allele && !inc_allele) {
        uint8_t dec_in_major = all_major & dec_allele;
        if (dec_in_major) {
            uint8_t not_decremented = all_major & ~dec_allele;
            if (!not_decremented) {
                score_change++;
            }
        }
    } else if (dec_allele && inc_allele) {
        uint8_t major_inc = all_major & inc_allele;
        uint8_t major_dec = all_major & dec_allele;
        if (major_inc) {
            score_change--;
        } else if (major_dec) {
            uint8_t not_dec = all_major & ~dec_allele;
            if (!(bnd1 & inc_allele) && !not_dec) {
                score_change++;
            }
        }
    }

    uint8_t new_major, new_bnd1;
    int children_delta;
    computeFitchChangeFromStored(all_major, bnd1, dec_allele, inc_allele,
                                  new_major, new_bnd1, children_delta);
    if (par_state) {
        bool old_follows = (all_major & par_state) != 0;
        bool new_follows = (new_major & par_state) != 0;
        if (!old_follows && new_follows) {
            score_change -= 1;
        } else if (old_follows && !new_follows) {
            score_change += 1;
        }
    }

    return score_change;
}

void SPRMutationOps::checkParsimonyScoreChangeAboveLCA(
    PhyloNode* start_node,
    int& parsimony_score_change,
    const MutationCountChangeCollection& allele_changes,
    MutationCountChangeCollection& parent_changes) {
    if (!start_node) {
        parent_changes.clear();
        return;
    }

    MutationCountChangeCollection current_changes = allele_changes;
    PhyloNode* node = start_node;
    PhyloNode* parent = getParent(node);
    PhyloNode* child_for_path = nullptr;  // unknown for first iteration

    while (parent && !current_changes.empty()) {
        MutationCountChangeCollection next_changes;

        getIntermediateNodesMutations(node, current_changes,
                                     next_changes, parsimony_score_change,
                                     child_for_path);

        child_for_path = node;  // track for subsequent iterations
        current_changes = next_changes;
        node = parent;
        parent = getParent(node);
    }

    // Hit root: remaining changes that never matched a node mutation
    // Use scoreDeltaInternal for each
    if (!parent && !current_changes.empty()) {
        for (const auto& c : current_changes) {
            parsimony_score_change += c.scoreDeltaInternal();
        }
        current_changes.clear();
    }

    parent_changes = current_changes;
}

// --- MutationCountChangeUtils ---

namespace MutationCountChangeUtils {

MutationCountChangeCollection merge_sorted(
    const MutationCountChangeCollection& a,
    const MutationCountChangeCollection& b) {

    MutationCountChangeCollection result;
    result.reserve(a.size() + b.size());

    auto it_a = a.begin();
    auto it_b = b.begin();

    while (it_a != a.end() && it_b != b.end()) {
        if (it_a->position < it_b->position) {
            result.push_back(*it_a); ++it_a;
        } else if (it_b->position < it_a->position) {
            result.push_back(*it_b); ++it_b;
        } else {
            MutationCountChange combined = *it_a;
            combined.removed_alleles |= it_b->removed_alleles;
            combined.added_alleles |= it_b->added_alleles;
            result.push_back(combined);
            ++it_a; ++it_b;
        }
    }

    while (it_a != a.end()) { result.push_back(*it_a); ++it_a; }
    while (it_b != b.end()) { result.push_back(*it_b); ++it_b; }

    return result;
}

MutationCountChangeCollection::const_iterator find_position(
    const MutationCountChangeCollection& collection,
    int position) {
    return std::lower_bound(collection.begin(), collection.end(), position,
        [](const MutationCountChange& mcc, int pos) {
            return mcc.position < pos;
        });
}

bool is_sorted(const MutationCountChangeCollection& collection) {
    for (size_t i = 1; i < collection.size(); ++i) {
        if (collection[i-1].position >= collection[i].position) return false;
    }
    return true;
}

} // namespace MutationCountChangeUtils
