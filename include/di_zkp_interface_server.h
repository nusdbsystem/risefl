//
// Created by yizheng on 15/3/23.
//

#ifndef DI_ZKP_CRYPTO_DI_ZKP_INTERFACE_SERVER_H
#define DI_ZKP_CRYPTO_DI_ZKP_INTERFACE_SERVER_H

#include <vector>
#include <unordered_map>
#include "di_zkp_interface_common.h"
#include "utils.h"
#include "base64.h"

struct RistElemHasher {
    std::size_t operator()(const RistElem &h) const;
};

constexpr int MAX_MULT_BASE_TABLE_BIT_SIZE = 24;

class ServerInterface : public CommonInterface {
public:
    // (num_clients + 1, num_clients + 1)
    std::vector<std::vector<int>> flags_collection; // flags_collection[i][j] == 1 means that i has flagged j as malicious

    // (num_clients + 1, dim)
    std::vector<RistElemP3Vec> committed_updates_collection;

    // (num_clients + 1, num_clients + 1, num_blinds_per_group_element, SCALBYTES)
    std::vector<std::vector<std::vector<CipherWithNonce>>> encrypted_shamir_shares_collection; // row [i] is received from client i

    // (num_clients + 1)
    std::vector<Proof> proof_collection;

    // (num_norm_bound_samples + 1)
    RistScalVec bb;

    //
    RistScal bb_bias;

    // (dim)
    RistScalVec aa_bb;

    // (num_clients + 1, num_clients + 1)
    std::vector<std::vector<int>> dispute_table; // row [i] is received from client i

    // (num_clients + 1, num_clients + 1, num_blinds_per_group_element)
    std::vector<RistScalMat> dispute_shares_collection; // row [i] is received from client i

    // (num_clients + 1, num_blinds_per_group_element)
    RistScalMat aggregates_collection;

    int small_mult_base_table_bit_size; // small multiplications of g table has default size 2^(min(MAX_MULT_BASE_TABLE_BIT_SIZE, weight_bits))
    std::unordered_map<RistElem, long, RistElemHasher> small_mult_base_table;

    // (dim)
    std::vector<long> final_update_int;

    // (dim)
    std::vector<float> final_update_float;
    std::vector<float> final_update_float_avg;

    std::vector<std::vector<long>> updates_int_collection;
    std::vector<std::vector<float>> updates_float_collection;
    double float_l2_sq_multiplier;

    CheckParamFloat check_param_float;

    ServerInterface(int num_clients, int max_malicious_clients,
                    int dim, int num_blinds_per_group_element,
                    int weight_bits, int random_normal_bit_shifter,
                    int num_norm_bound_samples, int inner_prod_bound_bits, int max_bound_sq_bits,
                    CHECK_TYPE check_type,
                    bool b_precomp = true,
                    PROTOCOL_TYPE protocol_type = PROTOCOL_TYPE::PRIV) :
            CommonInterface(num_clients, max_malicious_clients,
                            dim, num_blinds_per_group_element,
                            weight_bits, random_normal_bit_shifter,
                            num_norm_bound_samples, inner_prod_bound_bits, max_bound_sq_bits,
                            check_type,
                            b_precomp, protocol_type),
            flags_collection(num_clients + 1, std::vector<int>(num_clients + 1, 0)),
            committed_updates_collection(num_clients + 1, RistElemP3Vec(dim)),
            encrypted_shamir_shares_collection(num_clients + 1,
                                               std::vector<std::vector<CipherWithNonce>>(num_clients + 1,
                                                                                         std::vector<CipherWithNonce>(
                                                                                                 num_blinds_per_group_element,
                                                                                                 CipherWithNonce(
                                                                                                         RISTSCALBYTES)))),
            proof_collection(num_clients + 1,
                             Proof(check_type, num_blinds_per_group_element, num_norm_bound_samples, inner_prod_bound_bits,
                                   max_bound_sq_bits)),
            bb((check_type == CHECK_TYPE::L2NORM || check_type == CHECK_TYPE::SPHERE) ?
                num_norm_bound_samples + 1 :
                num_norm_bound_samples + 2),
            aa_bb(dim),
            dispute_table(num_clients + 1, std::vector<int>(num_clients + 1, 0)),
            dispute_shares_collection(num_clients + 1,
                                      RistScalMat(num_clients + 1,
                                                  RistScalVec(num_blinds_per_group_element))),
            aggregates_collection(num_clients + 1, RistScalVec(num_blinds_per_group_element)),
            final_update_int(dim),
            final_update_float(dim),
            final_update_float_avg(dim),
            small_mult_base_table_bit_size(
                    (weight_bits < MAX_MULT_BASE_TABLE_BIT_SIZE) ? weight_bits : MAX_MULT_BASE_TABLE_BIT_SIZE),
            check_param_float(check_type) {
        if (protocol_type == PROTOCOL_TYPE::NON_PRIV_INT) {
            updates_int_collection = std::vector<std::vector<long>>(num_clients + 1,
                                                                    std::vector<long>(dim, 0));
        }
        if (protocol_type == PROTOCOL_TYPE::NON_PRIV_FLOAT) {
            updates_float_collection = std::vector<std::vector<float>>(num_clients + 1,
                                                                       std::vector<float>(dim, 0));
        }
    }

    void import_small_mult_base_table(const std::unordered_map<RistElem, long, RistElemHasher> &small_table) {
        small_mult_base_table = small_table;
    }

    void generate_small_mult_base_table() {
        auto curr = scalar_mult_base(RistScal(-(1 << (small_mult_base_table_bit_size - 1))));
        auto g = scalar_mult_base(c_scal_one);
        for (int i = -(1 << (small_mult_base_table_bit_size - 1));
             i < (1 << (small_mult_base_table_bit_size - 1)); i++) {
            small_mult_base_table.emplace(curr, i);
            add(curr, curr, g);
        }
    }

    void reset_dispute_table() {
        dispute_table = std::vector<std::vector<int>>(num_clients + 1, std::vector<int>(num_clients + 1, 0));
    }

//    void load_encrypted_shamir_shares_and_check_string(const std::vector<unsigned char>& s, int i);

//    void transpose_encrypted_shamir_shares();

    void import_committed_updates(std::vector<unsigned char> pp, int i) {
        p3_vec_from_bytes(committed_updates_collection[i], pp);
    }

    void generate_aa_seed() {
        rand_init(predicate.aa_seed);
    }

//    void compute_hh_comb() {
//        predicate.cofmpute_hh_comb();
//    }

    void generate_hh_comb_and_bb_and_aabb();

    void generate_dispute_table();

    void check_linear_comb_batch_commitments_with_bb(int i);

    void check_sq_bound_proof(int i);

    void check_proof(int i);

    void check_proofs();

    void check_disputes();

//    void transpose_dispute_shares();

    // if parallel_on_clients, parallel on clients, sequential on dim
    // otherwise, parallel on dim, sequential on clients
    void compute_final_update(bool parallel_on_clients = false);

    /* abstractions */

    // used at 1st iteration
    void initialize_from_seed(const std::string &seed) {
        std::vector<unsigned char> seed_bytes = base64_decode(seed);
        predicate.initialize_from_seed(seed_bytes);
        generate_small_mult_base_table();
    }

//    // used at every iteration
//    void initialize_new_iteration(const RistScal &B) {
//        reset_server_flags();
//        reset_dispute_table();
//        set_l2_bound_sq(B);
//    }

    void initialize_new_iteration(const CheckParamFloat &check_param_float) {
        if (protocol_type == PROTOCOL_TYPE::NON_PRIV_FLOAT) {
            this->check_param_float = check_param_float;
            generate_aa_seed();
            final_update_float = std::vector<float>(dim, 0);
            normalizing_factor = 1.0;
            float_l2_sq_multiplier = gammas_from_num_samples[predicate.num_samples - 1];
//            float_l2_sq_multiplier = predicate.num_samples + standard_deviation_factor * sqrt(2 * predicate.num_samples);
            return;
        }

        auto check_param = set_normalizing_factor_and_compute_check_param_from_gamma(check_param_float);
        reset_server_flags();
        reset_dispute_table();
        set_param(check_param);
//        initialize_new_iteration(B);

        if (protocol_type == PROTOCOL_TYPE::NON_PRIV_INT) {
            generate_aa_seed();
            final_update_int = std::vector<long>(dim, 0);
        }

    }

    void receive_1_bytes(const std::vector<unsigned char> &bytes, int i);

    void receive_1(const std::string &bytes_str,
                   int i);

    std::vector<unsigned char> send_2_bytes();

    std::string send_2();

    void receive_2_bytes(const std::vector<unsigned char> &bytes, int i);

    void receive_2(const std::string &bytes_str,
                   int i);

    void concurrent_process_before_send_3();

    std::vector<unsigned char> send_3_bytes(int i);

    std::string send_3(int i);

    void receive_3_bytes(const std::vector<unsigned char> &bytes, int i);

    void receive_3(const std::string &bytes_str,
                   int i);

    void process_before_send_4();

    std::vector<unsigned char> send_4_bytes(int i);

    std::string send_4(int i);

    void receive_4_bytes(const std::vector<unsigned char> &bytes, int i);

    void receive_4(const std::string &bytes_str,
                   int i);

    void process_before_send_5();

    std::vector<unsigned char> send_5_bytes(int i);

    std::string send_5(int i);

    void receive_5_bytes(const std::vector<unsigned char> &bytes, int i);

    void receive_5(const std::string &bytes_str,
                   int i);

    // if parallel_on_clients, parallel on clients, sequential on dim
    // otherwise, parallel on dim, sequential on clients
    // result stored in final_update_float
    void finish_iteration(bool parallel_on_clients = false);

    std::string string_api_test(const std::string &a);
};

long discrete_log(const RistElem &y,
                  const std::unordered_map<RistElem, long, RistElemHasher> &small_table,
                  int per_side_step_count);

#endif //DI_ZKP_CRYPTO_DI_ZKP_INTERFACE_SERVER_H
