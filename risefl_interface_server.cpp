//
// Created by yizheng on 15/3/23.
//

#include <unordered_map>
#include <random>
#include "immintrin.h"
#include "include/risefl_interface_server.h"
#include "include/risefl_interface_common.h"
#include "include/shamir.h"
#include "include/exception.h"
#include "include/ristretto_vector.h"
#include "include/random_avx/mersenne-twister-avx2.h"

void ServerInterface::generate_hh_comb_and_bb_and_aabb() {
    bb_bias = c_scal_zero;

    rand_init(bb);
    auto check_type = predicate.check_param.check_type;

    auto seeds_at_a_m = generate_multiple_seeds(predicate.aa_seed, predicate.num_samples + 1);

//    std::mt19937_64 generator;
//    generator.seed(seeds_at_a_m[0]);
    RistScalVec top_row(dim);
    for (int l = 0; l < dim; l++) {
//        std::uniform_int_distribution<unsigned long> unif_unsigned_long;
//        RistHashbytes rand_hash;
//        for (int j = 0; j < rand_hash.hashbytes.size(); j += sizeof(unsigned long)) {
//            auto u = unif_unsigned_long(generator);
//            memcpy(rand_hash.hashbytes.data() + j, &u, sizeof(unsigned long));
//        }
//        scalar_reduce(top_row[l], rand_hash);

        top_row[l] = rist_scalar_from_hash_from_bytes(predicate.aa_seed.hashbytes, bytes_from_int(l));
        aa_bb[l] = top_row[l] * bb[0];
        if (check_type == CHECK_TYPE::COSINE_SIM) {
            aa_bb[l] += predicate.check_param.cosine_param.pivot[l] * bb[bb.size() - 1];
        }
    }
    linear_comb_block(predicate.hh_comb.elems[0], top_row, predicate.hh);
    if (check_type == CHECK_TYPE::COSINE_SIM) {
        linear_comb_block(predicate.hh_comb.elems[predicate.hh_comb.elems.size() - 1],
                          predicate.check_param.cosine_param.pivot,
                          predicate.hh);
    }



    std::vector<std::vector<int>> temp_A(NUM_THREADS, std::vector<int>(dim));
//    std::vector<std::vector<float>> temp_A_float(NUM_THREADS, std::vector<float>(dim));

    std::vector<IntCombRistScal> aa_bb_ttmath(dim);
    rist_to_ttmath(aa_bb_ttmath, aa_bb);
//    MpzVec bb_mpz(predicate.num_samples + 1);
//    mpz_vec_from_rist(bb_mpz, bb);


    for (int m = 1; m <= predicate.num_samples; m += NUM_THREADS) {
        // TODO: non-AVX2
        assert(__AVX2__);

        std::for_each(
                std::execution::par_unseq,
                temp_A.begin(),
                temp_A.end(),
                [this, &m, &temp_A, &seeds_at_a_m](auto &&temp_A_row) {
                    auto j = &temp_A_row - temp_A.data();
                    if (m + j <= predicate.num_samples) {
                        mt_avx2 rand_gen(seeds_at_a_m[m + j]);
//                        rand_gen.rand_normal_avx2(temp_A_float[j].data(), dim);
//                        for (int l = 0; l < dim; l++) {
//                            temp_A[j][l] = std::lround(temp_A_float[j][l] * (1 << predicate.random_normal_bit_shifter));
//                        }

                        rand_gen.rand_normal_avx2_shifted(temp_A[j].data(), dim, predicate.random_normal_bit_shifter);
                        linear_comb_block(predicate.hh_comb.elems[m + j],
                                          temp_A[j],
                                          predicate.random_normal_bit_shifter + 3,
                                          predicate.hh);
                    }
                }
        );


        std::for_each(std::execution::par_unseq, aa_bb_ttmath.begin(), aa_bb_ttmath.end(),
                      [this, &m, &aa_bb_ttmath, &temp_A, &check_type](auto &&aa_bb_ttmath_l) {
                          auto l = &aa_bb_ttmath_l - aa_bb_ttmath.data();
                          for (int j = 0; (j < NUM_THREADS) && (m + j <= predicate.num_samples); j++) {
                              if (temp_A[j][l] >= 0)
                                  aa_bb_ttmath_l += multiply_nonnegint(temp_A[j][l], bb[m + j]);
                              else
                                  aa_bb_ttmath_l -= multiply_nonnegint(-temp_A[j][l], bb[m + j]);
                              if (check_type == CHECK_TYPE::SPHERE) {
                                  bb_bias += RistScal(-temp_A[j][l]) * bb[m + j] * predicate.check_param.sphere_param.center[l];
                              }
                          }
                      });
    }

    predicate.hh_comb.fill_bytes();
    ttmath_to_rist(aa_bb, aa_bb_ttmath);
}

void ServerInterface::generate_dispute_table() {
    for (int i = 1; i <= num_clients; i++) {
        // if more than max_malicious_clients clients flag i, flag i
        int count = 0;
        for (int j = 1; j <= num_clients; j++) {
            if (flags_collection[i][j])
                count++;
        }
        if (count > max_malicious_clients) {
            server_flags[i] = 1;
            continue;
        }


        // if i flags more than max_malicious_clients clients, flag i
        count = 0;
        for (int j = 1; j <= num_clients; j++) {
            if (flags_collection[j][i])
                count++;
        }
        if (count > max_malicious_clients) {
            server_flags[i] = 1;
            continue;
        }
    }


    // create dispute table
    for (int i = 1; i <= num_clients; i++) {
        if (server_flags[i])
            continue;
        for (int j = 1; j <= num_clients; j++) {
            if (server_flags[j])
                continue;
            if (flags_collection[j][i]) {
                dispute_table[i][j] = 1;
            }
        }
    }
}

void ServerInterface::check_disputes() {
    for (int i = 1; i <= num_clients; i++) {
        if (server_flags[i])
            continue;
        for (int j = 1; j <= num_clients; j++) {
            if (dispute_table[i][j] && !b_shamir_verify(dispute_shares_collection[i][j],
                                                        shamir_check_string_collection[i].elems,
                                                        max_malicious_clients + 1,
                                                        j)) {
                server_flags[i] = 1;
                break;
            }
        }
    }
}

std::size_t RistElemHasher::operator()(const RistElem &h) const {
    std::size_t ret;
    memcpy(&ret, h.element.data(), sizeof(ret));
    return ret;
}


// compute discrete log of y with base g using baby-step giant-step
// from small table
long discrete_log(const RistElem &y,
                  const std::unordered_map<RistElem, long, RistElemHasher> &small_table,
                  long per_side_step_count) {
    long step = small_table.size();
    auto curr_pos = y;
    auto curr_neg = y;
    auto val = small_table.find(curr_pos);
    if (val != small_table.end()) {
        return val->second;
    }
    RistElemP3 gc_p3;
    pedersen_zero_commit(gc_p3, RistScal(step));
    RistElemCached gc_cached;
    p3_to_cached(gc_cached, gc_p3);
//    RistElem gc = scalar_mult_base(RistScal(step));
    RistElemP3 curr_pos_p3, curr_neg_p3;
    p3_from_bytes(curr_pos_p3, curr_pos.element);
    p3_from_bytes(curr_neg_p3, curr_neg.element);
    for (long i = 1; i <= per_side_step_count; i++) {
        p3_sub(curr_pos_p3, curr_pos_p3, gc_cached);
        bytes_from_p3(curr_pos.element, curr_pos_p3);
        val = small_table.find(curr_pos);
        if (val != small_table.end()) {
            return i * step + val->second;
        }
        p3_add(curr_neg_p3, curr_neg_p3, gc_cached);
        bytes_from_p3(curr_neg.element, curr_neg_p3);
        val = small_table.find(curr_neg);
        if (val != small_table.end()) {
            return -i * step + val->second;
        }
    }
    throw Abort(); // too many malicious clients, should abort, throw exception to python
}

class CommittedUpdatesIterator : public std::iterator<std::input_iterator_tag, RistElemP3Vec> {
public:
    const RistElemP3Mat &committed_updates_collection;
    const std::vector<int> &client_indices;
    int index;

    explicit CommittedUpdatesIterator(const RistElemP3Mat &committed_updates_collection,
                                      const std::vector<int> &client_indices,
                                      int index = 0) :
            committed_updates_collection(committed_updates_collection),
            client_indices(client_indices),
            index(index) {}

    CommittedUpdatesIterator &operator++() {
        index++;
        return *this;
    }

    bool operator==(CommittedUpdatesIterator other) const { return index == other.index; }

    bool operator!=(CommittedUpdatesIterator other) const { return !(*this == other); }

    CommittedUpdatesIterator begin() {
        return CommittedUpdatesIterator(committed_updates_collection, client_indices, 0);
    };

    CommittedUpdatesIterator end() {
        return CommittedUpdatesIterator(committed_updates_collection, client_indices, client_indices.size());
    };

    reference
    operator*() const { return const_cast<RistElemP3Vec &>(committed_updates_collection[client_indices[index]]); }
};

void ServerInterface::compute_final_update(bool parallel_on_clients) {
    auto check_type = predicate.check_param.check_type;

    if (protocol_type == PROTOCOL_TYPE::NON_PRIV_INT) {

        RistScalVec sum_inner_prod_sq_collection(num_clients + 1, c_scal_zero);

        auto seeds_at_a_m = generate_multiple_seeds(predicate.aa_seed, predicate.num_samples + 1);
        std::vector<int> temp(dim);

        if (check_type == CHECK_TYPE::L2NORM) {
            for (int m = 1; m <= predicate.num_samples; m++) {
                // TODO: non-AVX2
                assert(__AVX2__);
                mt_avx2 rand_gen(seeds_at_a_m[m]);
                rand_gen.rand_normal_avx2_shifted(temp.data(), dim, predicate.random_normal_bit_shifter);
                for (int i = 1; i <= num_clients; i++) {
                    RistScal ip_temp(inner_prod(temp, updates_int_collection[i]));
                    sum_inner_prod_sq_collection[i] += ip_temp * ip_temp;
                }
            }

            for (int i = 1; i <= num_clients; i++) {
                if (get_power_two_bound(predicate.check_param.l2_param.bound_sq - sum_inner_prod_sq_collection[i]) <=
                    predicate.max_bound_sq_bits) {
                    // i is valid
                    for (int j = 0; j < dim; j++) {
                        final_update_int[j] += updates_int_collection[i][j];
                    }
                } else {
                    // i is invalid
                    server_flags[i] = 1;
                }
            }
        }

        if (check_type == CHECK_TYPE::SPHERE) {
            for (int m = 1; m <= predicate.num_samples; m++) {
                // TODO: non-AVX2
                assert(__AVX2__);
                mt_avx2 rand_gen(seeds_at_a_m[m]);
                rand_gen.rand_normal_avx2_shifted(temp.data(), dim, predicate.random_normal_bit_shifter);
                for (int i = 1; i <= num_clients; i++) {
                    RistScal ip_temp(inner_prod(temp, rist_scal_vec_from_long_vec(updates_int_collection[i]) - predicate.check_param.sphere_param.center));
                    sum_inner_prod_sq_collection[i] += ip_temp * ip_temp;
                }
            }

            for (int i = 1; i <= num_clients; i++) {
                if (get_power_two_bound(predicate.check_param.sphere_param.bound_sq - sum_inner_prod_sq_collection[i]) <=
                    predicate.max_bound_sq_bits) {
                    // i is valid
                    for (int j = 0; j < dim; j++) {
                        final_update_int[j] += updates_int_collection[i][j];
                    }
                } else {
                    // i is invalid
                    server_flags[i] = 1;
                }
            }
        }

        if (check_type == CHECK_TYPE::COSINE_SIM) {
            RistScalVec pivot_inner_prod_collection(num_clients + 1, c_scal_zero);
            auto pivot = predicate.check_param.cosine_param.pivot;
            auto inner_prod_sq_multiplier = predicate.check_param.cosine_param.inner_prod_sq_multiplier;
            auto l2_sq_multiplier = predicate.check_param.cosine_param.l2_sq_multiplier;
            auto pivot_norm_sq = inner_prod(pivot, pivot);

            int weight_bit_shifter = weight_bits - 1;

            for (int m = 1; m <= predicate.num_samples; m++) {
                // TODO: non-AVX2
                assert(__AVX2__);
                mt_avx2 rand_gen(seeds_at_a_m[m]);
                rand_gen.rand_normal_avx2_shifted(temp.data(), dim, predicate.random_normal_bit_shifter);
                for (int i = 1; i <= num_clients; i++) {
                    RistScal ip_temp(inner_prod(temp, updates_int_collection[i]));
                    sum_inner_prod_sq_collection[i] += ip_temp * ip_temp;
                    pivot_inner_prod_collection[i] = inner_prod(updates_int_collection[i], pivot);
                }
            }

            for (int i = 1; i <= num_clients; i++) {
                if (get_power_two_bound(predicate.check_param.cosine_param.bound_sq - sum_inner_prod_sq_collection[i]) <=
                    predicate.max_bound_sq_bits
                    &&
                    get_power_two_bound(pivot_inner_prod_collection[i]) <= predicate.max_bound_sq_bits
                    &&
                    get_power_two_bound(inner_prod_sq_multiplier * pivot_inner_prod_collection[i] * pivot_inner_prod_collection[i]
                        - l2_sq_multiplier * sum_inner_prod_sq_collection[i] * pivot_norm_sq) <= predicate.max_bound_sq_bits
                    ) {
                    // i is valid
                    for (int j = 0; j < dim; j++) {
                        final_update_int[j] += updates_int_collection[i][j];
                    }
                } else {
                    // i is invalid
                    server_flags[i] = 1;
                }
            }
        }

        for (int i = 0; i < dim; i++) {
            final_update_float[i] =
                    static_cast<double>(final_update_int[i]) / (1L << (weight_bits - 1)) * normalizing_factor;
            final_update_float_avg[i] = final_update_float[i] / valid_client_count();
        }
        return;
    }

    if (protocol_type == PROTOCOL_TYPE::NON_PRIV_FLOAT) {
        std::vector<double> sum_inner_prod_sq_collection(num_clients + 1, 0);

        auto seeds_at_a_m = generate_multiple_seeds(predicate.aa_seed, predicate.num_samples + 1);
        std::vector<float> temp(dim);

        if (check_type == CHECK_TYPE::L2NORM) {
            for (int m = 1; m <= predicate.num_samples; m++) {
                // TODO: non-AVX2
                assert(__AVX2__);
                mt_avx2 rand_gen(seeds_at_a_m[m]);
                rand_gen.rand_normal_avx2(temp.data(), dim);
                for (int i = 1; i <= num_clients; i++) {
                    double ip_temp(inner_prod(temp, updates_float_collection[i]));
                    sum_inner_prod_sq_collection[i] += ip_temp * ip_temp;
                }
            }

            for (int i = 1; i <= num_clients; i++) {
                if (sum_inner_prod_sq_collection[i] <= float_l2_sq_multiplier * check_param_float.l2_param.bound * check_param_float.l2_param.bound) {
                    // i is valid
                    for (int j = 0; j < dim; j++) {
                        final_update_float[j] += updates_float_collection[i][j];
                    }
                } else {
                    // i is invalid
                    server_flags[i] = 1;
                }
            }
        }

        if (check_type == CHECK_TYPE::SPHERE) {
            for (int m = 1; m <= predicate.num_samples; m++) {
                // TODO: non-AVX2
                assert(__AVX2__);
                mt_avx2 rand_gen(seeds_at_a_m[m]);
                rand_gen.rand_normal_avx2(temp.data(), dim);
                for (int i = 1; i <= num_clients; i++) {
                    auto diff = updates_float_collection[i];
                    for (int l = 0; l < dim; l++) {
                        diff[l] -= check_param_float.sphere_param.center[l];
                    }
                    double ip_temp(inner_prod(temp, diff));
                    sum_inner_prod_sq_collection[i] += ip_temp * ip_temp;
                }
            }

            for (int i = 1; i <= num_clients; i++) {
                if (sum_inner_prod_sq_collection[i] <= float_l2_sq_multiplier *
                    check_param_float.sphere_param.bound * check_param_float.sphere_param.bound) {
                    // i is valid
                    for (int j = 0; j < dim; j++) {
                        final_update_float[j] += updates_float_collection[i][j];
                    }
                } else {
                    // i is invalid
                    server_flags[i] = 1;
                }
            }
        }

        if (check_type == CHECK_TYPE::COSINE_SIM) {
            std::vector<double> pivot_inner_prod_collection(num_clients + 1, 0);
            auto pivot = check_param_float.cosine_param.pivot;
            auto cosine_bound = check_param_float.cosine_param.cosine_bound;
            double cosine_bound_sq = cosine_bound * cosine_bound;
            double pivot_norm_sq = inner_prod(pivot, pivot);

            for (int m = 1; m <= predicate.num_samples; m++) {
                // TODO: non-AVX2
                assert(__AVX2__);
                mt_avx2 rand_gen(seeds_at_a_m[m]);
                rand_gen.rand_normal_avx2(temp.data(), dim);
                for (int i = 1; i <= num_clients; i++) {
                    double ip_temp(inner_prod(temp, updates_float_collection[i]));
                    sum_inner_prod_sq_collection[i] += ip_temp * ip_temp;
                    pivot_inner_prod_collection[i] = inner_prod(updates_float_collection[i], pivot);
                }
            }

            for (int i = 1; i <= num_clients; i++) {
                if (sum_inner_prod_sq_collection[i] <= float_l2_sq_multiplier &&
                    pivot_inner_prod_collection[i] > 0 &&
                    cosine_bound_sq * sum_inner_prod_sq_collection[i] * pivot_norm_sq <= float_l2_sq_multiplier * pivot_inner_prod_collection[i] * pivot_inner_prod_collection[i]) {
                    // i is valid
                    for (int j = 0; j < dim; j++) {
                        final_update_float[j] += updates_float_collection[i][j];
                    }
                } else {
                    // i is invalid
                    server_flags[i] = 1;
                }
            }
        }

        for (int j = 0; j < dim; j++) {
//            final_update_float[j] *= normalizing_factor; // normalizing_factor is 1.0 by default with NON_PRIV_FLOAT
            final_update_float_avg[j] = final_update_float[j] / valid_client_count();
        }
        return;
    }


    std::vector<int> good_clients;
    RistScalMat good_aggs(predicate.num_blinds_per_group_element);

    RistElemP3Mat valid_check_strings(predicate.num_blinds_per_group_element,
                                      RistElemP3Vec(max_malicious_clients + 1));
    for (int k = 0; k < predicate.num_blinds_per_group_element; k++) {
        std::vector<RistElemP3Vec> valid_check_strings_to_sum(max_malicious_clients + 1);
        for (int i = 1; i <= num_clients; i++) {
            if (server_flags[i])
                continue;
            for (int j = 0; j <= max_malicious_clients; j++) {
                valid_check_strings_to_sum[j].emplace_back(shamir_check_string_collection[i].elems[k][j]);
            }
        }
        for (int j = 0; j <= max_malicious_clients; j++) {
            sum_single_thread(valid_check_strings[k][j], valid_check_strings_to_sum[j]);
        }
    }

    for (int i = 1; i <= num_clients && good_clients.size() <= max_malicious_clients; i++) {
        if (server_flags[i])
            continue;
        if (!b_shamir_verify(aggregates_collection[i], valid_check_strings, max_malicious_clients + 1, i)) {
            // i submits an invalid aggregate, but its shares are included in other clients' updates, so don't flag i
            continue;
        }
        good_clients.emplace_back(i);
        for (int k = 0; k < predicate.num_blinds_per_group_element; k++)
            good_aggs[k].emplace_back(aggregates_collection[i][k]);
    }
    if (good_clients.size() < max_malicious_clients + 1)
        throw Abort(); // too many malicious clients, should abort, throw exception to python

    RistScalVec recovered_share_sum = shamir_recover(good_aggs, good_clients, max_malicious_clients + 1);

    // new algo
    if (parallel_on_clients) {
        std::vector<int> indices_to_agg;
        for (int i = 1; i <= num_clients; i++) {
            if (!server_flags[i])
                indices_to_agg.emplace_back(i);
        }
        CommittedUpdatesIterator it(committed_updates_collection, indices_to_agg);

        RistElemP3Vec powers_init(dim);
        std::for_each(std::execution::par_unseq, powers_init.begin(), powers_init.end(),
                      [this, &powers_init, &recovered_share_sum](auto &&power) {
                          auto l = &power - powers_init.data();
                          power = -recovered_share_sum[l % predicate.num_blinds_per_group_element] *
                                  predicate.hh[l / predicate.num_blinds_per_group_element];
                      });
        RistElemP3Vec powers_p3 = std::reduce(std::execution::par_unseq, it.begin(), it.end(), powers_init,
                                              [](const RistElemP3Vec &vv1, const RistElemP3Vec &vv2) {
                                                  RistElemP3Vec vv(vv1.size());
                                                  for (int i = 0; i < vv1.size(); i++) {
                                                      vv[i] = vv1[i] + vv2[i];
                                                  }
                                                  return vv;
                                              });

        std::for_each(
                std::execution::par_unseq,
                final_update_int.begin(),
                final_update_int.end(),
                [this, &recovered_share_sum, &powers_p3](auto &&u_sum) {
                    auto l = &u_sum - this->final_update_int.data();
                    RistElem power;
                    bytes_from_p3(power.element, powers_p3[l]);

                    u_sum = discrete_log(power,
                                         this->small_mult_base_table,
                                         this->num_clients * (1 << (weight_bits - small_mult_base_table_bit_size)));
                }
        );
    } else {
        std::for_each(
                std::execution::par_unseq,
                final_update_int.begin(),
                final_update_int.end(),
                [this, &recovered_share_sum](auto &&u_sum) {
                    auto l = &u_sum - this->final_update_int.data();
                    RistElemP3 u_sum_p3;
                    p3_0(u_sum_p3);
                    for (int i = 1; i <= num_clients; i++) {
                        if (!server_flags[i])
                            u_sum_p3 += committed_updates_collection[i][l];
                    }
                    RistElemP3 power_p3 = u_sum_p3 -
                                          recovered_share_sum[l % predicate.num_blinds_per_group_element] *
                                          predicate.hh[l / predicate.num_blinds_per_group_element];
                    RistElem power;

                    bytes_from_p3(power.element, power_p3);

                    u_sum = discrete_log(power,
                                         this->small_mult_base_table,
                                         this->num_clients * (1 << (weight_bits - small_mult_base_table_bit_size)));
                }
        );
    }

    for (int i = 0; i < dim; i++) {
        final_update_float[i] =
                static_cast<double>(final_update_int[i]) / (1L << (weight_bits - 1)) * normalizing_factor;
        final_update_float_avg[i] = final_update_float[i] / valid_client_count();
    }
}


void ServerInterface::check_sq_bound_proof(int i) {
    if (!b_verify_pedersen_batch_with_zero(predicate.hh_comb, shamir_check_string_collection[i].slice_col_0(),
                                           proof_collection[i].linear_comb_batch_commitments, predicate.square_key,
                                           proof_collection[i].linear_comb_single_commitments,
                                           proof_collection[i].proof_well_formed)) {
        server_flags[i] = 1;
        return;
    }
    if (!b_verify_abs_val_range_power_two_agg(predicate.square_key,
                                              proof_collection[i].linear_comb_single_commitments.no_head(),
                                              predicate.bound_elem_keys_1,
                                              predicate.bound_elem_keys_2,
                                              predicate.inner_prod_bound_bits,
                                              predicate.num_samples,
                                              proof_collection[i].proof_linear_comb_bound)) {
        server_flags[i] = 1;
        return;
    }

    if (!b_verify_batch_pedersen_with_square(predicate.square_key,
                                             proof_collection[i].linear_comb_single_commitments.no_head(),
                                             proof_collection[i].square_commitments,
                                             proof_collection[i].proof_squares)) {
        server_flags[i] = 1;
        return;
    }

    RistP3AndBytes sum_of_square_commitments;
    sum_single_thread(sum_of_square_commitments.elem, proof_collection[i].square_commitments.elems);
    RistP3AndBytes bound_minus_sum_of_square_commitments;
    bound_minus_sum_of_square_commitments.elem =
            pedersen_zero_commit_p3(predicate.check_param.l2_param.bound_sq) - sum_of_square_commitments.elem;
    bound_minus_sum_of_square_commitments.fill_bytes();
    if (!b_verify_range_power_two(predicate.square_key,
                                  bound_minus_sum_of_square_commitments,
                                  predicate.bound_elem_keys_1,
                                  predicate.bound_elem_keys_2,
                                  predicate.max_bound_sq_bits,
                                  proof_collection[i].proof_sum_range))
        server_flags[i] = 1;
}


void ServerInterface::check_linear_comb_batch_commitments_with_bb(int i) {
    LinearCombCalculator lcc;
    lcc.add(bb, proof_collection[i].linear_comb_batch_commitments.elems);
    lcc.add_base(bb_bias);
    lcc.sub(aa_bb, committed_updates_collection[i]);
//    RistElemP3 left = linear_comb(bb, proof_collection[i].linear_comb_batch_commitments.elems);
//    RistElemP3 right = linear_comb(aa_bb, committed_updates_collection[i]);
//    if (left != right) {
    if (!lcc.is_0()) {
        server_flags[i] = 1;
    }
}

void ServerInterface::check_proof(int i) {
    check_sq_bound_proof(i);
    check_linear_comb_batch_commitments_with_bb(i);
}

void ServerInterface::check_proofs() {
    for (int i = 1; i <= num_clients; i++)
        check_proof(i);
}

void ServerInterface::receive_1_bytes(const std::vector<unsigned char> &bytes, int i) {
    if (protocol_type == PROTOCOL_TYPE::NON_PRIV_INT) {
        for (int j = 0; j < dim; j++) {
            memcpy(&updates_int_collection[i][j], bytes.data() + j * sizeof(long), sizeof(long));
        }
        return;
    }
    if (protocol_type == PROTOCOL_TYPE::NON_PRIV_FLOAT) {
        for (int j = 0; j < dim; j++) {
            memcpy(&updates_float_collection[i][j], bytes.data() + j * sizeof(float), sizeof(float));
        }
        return;
    }

    if (bytes.size() != SIGNBYTES + crypto_box_PUBLICKEYBYTES) {
        server_flags[i] = 1;
        return;
    }
    std::copy(bytes.begin(), bytes.end(), bul_signed_pub_key_collection[i].begin());
}

void ServerInterface::receive_1(const std::string &bytes_str, int i) {
    receive_1_bytes(base64_decode(bytes_str), i);
}

std::vector<unsigned char> ServerInterface::send_2_bytes() {
    std::vector<unsigned char> ret(num_clients *
    (SIGNBYTES + crypto_box_PUBLICKEYBYTES));
    auto it_ret = ret.begin();
    for (int i = 1; i <= num_clients; i++) {
        it_ret = std::copy(bul_signed_pub_key_collection[i].begin(),
                           bul_signed_pub_key_collection[i].end(),
                           it_ret);
//        dh_public_key_collection[i].export_to_bytestream(it_ret);
    }
    assert(it_ret == ret.end());
    return ret;
}

std::string ServerInterface::send_2() {
    return base64_encode(send_2_bytes());
}

void ServerInterface::receive_2_bytes(const std::vector<unsigned char> &bytes, int i) {
    int size_committed_updates = dim * RISTBYTES;
    int size_encrypted_shares = num_clients * predicate.num_blinds_per_group_element * ENCRISTSCALBYTES;
    int size_check_strings = predicate.num_blinds_per_group_element * (max_malicious_clients + 1) * RISTBYTES;

    if (bytes.size() != size_committed_updates + size_encrypted_shares + size_check_strings) {
        server_flags[i] = 1;
        return;
    }

    auto it_bytes = bytes.cbegin();
    p3_vec_from_bytestream(committed_updates_collection[i], it_bytes);

    for (int j = 1; j <= num_clients; j++) {
        for (int k = 0; k < predicate.num_blinds_per_group_element; k++)
            encrypted_shamir_shares_collection[i][j][k].import_from_bytestream(it_bytes);
    }

    shamir_check_string_collection[i].import_from_bytestream(it_bytes);

    assert(it_bytes == bytes.end());
}


void ServerInterface::receive_2(const std::string &bytes_str, int i) {
    receive_2_bytes(base64_decode(bytes_str), i);
}

void ServerInterface::concurrent_process_before_send_3() {
    generate_aa_seed();
//    generate_aa();
//    compute_hh_comb();
    generate_hh_comb_and_bb_and_aabb();
}

std::vector<unsigned char> ServerInterface::send_3_bytes(int i) {
    int size_encrypted_shares = num_clients * predicate.num_blinds_per_group_element * ENCRISTSCALBYTES;
    int size_check_strings =
            num_clients * predicate.num_blinds_per_group_element * (max_malicious_clients + 1) * RISTBYTES;
    int size_aa_seed = crypto_core_ed25519_HASHBYTES;
    int size_hh_comb = (predicate.num_samples + 1) * predicate.num_blinds_per_group_element * RISTBYTES;

    std::vector<unsigned char> ret(size_encrypted_shares + size_check_strings +
                                   size_aa_seed +
                                   size_hh_comb);

    auto it_ret = ret.begin();

    for (int j = 1; j <= num_clients; j++) {
        for (int k = 0; k < predicate.num_blinds_per_group_element; k++)
            encrypted_shamir_shares_collection[j][i][k].export_to_bytestream(it_ret);
    }

    for (int j = 1; j <= num_clients; j++) {
        shamir_check_string_collection[j].export_to_bytestream(it_ret);
    }

    predicate.aa_seed.export_to_bytestream(it_ret);

    predicate.hh_comb.export_to_bytestream(it_ret);

    assert(it_ret == ret.end());
    return ret;
}


std::string ServerInterface::send_3(int i) {
    return base64_encode(send_3_bytes(i));
}

void ServerInterface::receive_3_bytes(const std::vector<unsigned char> &bytes, int i) {
    if (bytes.size() != Proof::size(predicate.check_param.check_type,
                                    predicate.num_blinds_per_group_element,
                                    predicate.num_samples,
                                    predicate.inner_prod_bound_bits,
                                    predicate.max_bound_sq_bits)
                        + num_clients) {
        server_flags[i] = 1;
        return;
    }

    auto it_bytes = bytes.cbegin();

    proof_collection[i].import_from_bytestream(it_bytes);

    import_flags_start_from_1_from_bytestream(flags_collection[i], it_bytes);
    assert(it_bytes == bytes.end());

    check_proof(i);
}


void ServerInterface::receive_3(const std::string &bytes_str, int i) {
    receive_3_bytes(base64_decode(bytes_str), i);

}


void ServerInterface::process_before_send_4() {
    generate_dispute_table();
}

std::vector<unsigned char> ServerInterface::send_4_bytes(int i) {
    std::vector<unsigned char> ret(num_clients);
    auto it_ret = ret.begin();
    export_flags_start_from_1_to_bytestream(dispute_table[i], it_ret);
    assert(it_ret == ret.end());
    return ret;
}


std::string ServerInterface::send_4(int i) {
    return base64_encode(send_4_bytes(i));
}

void ServerInterface::receive_4_bytes(const std::vector<unsigned char> &bytes, int i) {
    if (bytes.size() != num_clients * predicate.num_blinds_per_group_element * RISTSCALBYTES) {
        server_flags[i] = 1;
        return;
    }
    auto it_bytes = bytes.cbegin();
    for (int j = 1; j <= num_clients; j++) {
        import_from_bytestream(dispute_shares_collection[i][j], it_bytes);
    }
    assert(it_bytes == bytes.end());
}


void ServerInterface::receive_4(const std::string &bytes_str, int i) {
    receive_4_bytes(base64_decode(bytes_str), i);
}

void ServerInterface::process_before_send_5() {
    check_disputes();
}


std::vector<unsigned char> ServerInterface::send_5_bytes(int i) {
    std::vector<unsigned char> ret(num_clients *predicate.num_blinds_per_group_element * RISTSCALBYTES + num_clients);
    auto it_ret = ret.begin();
    for (int j = 1; j <= num_clients; j++) {
        export_to_bytestream(dispute_shares_collection[j][i], it_ret);
    }
    export_flags_start_from_1_to_bytestream(server_flags, it_ret);
    assert(it_ret == ret.end());
    return ret;
}


std::string ServerInterface::send_5(int i) {
    return base64_encode(send_5_bytes(i));
}

void ServerInterface::receive_5_bytes(const std::vector<unsigned char> &bytes, int i) {
    if (bytes.size() != predicate.num_blinds_per_group_element * RISTSCALBYTES) {
        server_flags[i] = 1;
        return;
    }
    auto it_bytes = bytes.cbegin();
    import_from_bytestream(aggregates_collection[i], it_bytes);
    assert(it_bytes == bytes.end());
}


void ServerInterface::receive_5(const std::string &bytes_str, int i) {
    receive_5_bytes(base64_decode(bytes_str), i);
}

void ServerInterface::finish_iteration(bool parallel_on_clients) {
    compute_final_update(parallel_on_clients);
}

std::string ServerInterface::string_api_test(const std::string &a) {
    std::cout << "The input string is " << a << std::endl;
    std::string res = "Hey python";
    return res;
}