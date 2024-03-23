//
// Created by yizheng on 15/3/23.
//

// risefl_interface_server.h
// 
// the RiseFL server interface header file
// both server and client interfaces are built on the common interface

#ifndef RISEFL_CRYPTO_RISEFL_INTERFACE_SERVER_H
#define RISEFL_CRYPTO_RISEFL_INTERFACE_SERVER_H

#include <vector>
#include <unordered_map>
#include "risefl_interface_common.h"
#include "utils.h"
#include "base64.h"

// compute the hash of a Ristretto group element
// used in computing discrete log at the aggregation stage
struct RistElemHasher {
    std::size_t operator()(const RistElem &h) const;
};

// to compute the discrete log, we use the baby-step giant-step algorithm
// for the fixed Ristretto group generator g, we store the table of (x g) for all -2^(b-1) <= x < 2^(b-1),
// in the computation of discretelog(g, y), we check for a match of y + 2^b a g in the table for every integer a,
// where b = min(weight_bits, MAX_MULT_BASE_TABLE_BIT_SIZE)
// with value MAX_MULT_BASE_TABLE_BIT_SIZE = 24, the table takes at most ~512MB of RAM
constexpr int MAX_MULT_BASE_TABLE_BIT_SIZE = 24;


class ServerInterface : public CommonInterface {
public:
    // the collection of client flags
    // flags_collection[i][j] == 1 means that i has flagged j as malicious
    // size: (num_clients + 1, num_clients + 1)
    // flags_collection[0] and flags_collection[i][0] are never used
    std::vector<std::vector<int>> flags_collection; 

    // the collection of all the clients' committed model updates
    // size: (num_clients + 1, dim)
    // committed_updates_collection[i] is client i's committed model update
    // committed_updates_collection[0] is never used
    std::vector<RistElemP3Vec> committed_updates_collection;

    // the collection of all the clients' encrypted shamir shares
    // size: (num_clients + 1, num_clients + 1, num_blinds_per_group_element, SCALBYTES)
    // encrypted_shamir_shares_collection[i] is the encrypted shamir shares sent to client i
    // encrypted_shamir_shares_collection[0] is never used
    std::vector<std::vector<std::vector<CipherWithNonce>>> encrypted_shamir_shares_collection; 

    // the collection of all the clients' proofs
    // size: (num_clients + 1)
    // proof_collection[i] is client i's proof
    // proof_collection[0] is never used
    std::vector<Proof> proof_collection;

    // the random numbers to check the correctness of clients' commitments of inner products <a_t, u>
    // size: (num_norm_bound_samples + 1)
    RistScalVec bb;

    // helper number computed from bb and random normal samples to check the correctness of clients' commitments if using sphere check
    RistScal bb_bias;

    // helper numbers computed from bb and random normal samples to check the correctness of clients' commitments if using sphere check
    // size: (dim)
    RistScalVec aa_bb;

    // the table of disputes
    // size: (num_clients + 1, num_clients + 1)
    // dispute_table[i] is the list of clients that marks i as malicious
    // dispute_table[0] is never used
    std::vector<std::vector<int>> dispute_table; 

    // the collection of shares sent at the dispute stage
    // size: (num_clients + 1, num_clients + 1, num_blinds_per_group_element)
    // dispute_shares_collection[i] is the clear-text shares sent by client i
    // dispute_shares_collection[0] is never used
    std::vector<RistScalMat> dispute_shares_collection;

    // the collection of aggregated shares
    // size: (num_clients + 1, num_blinds_per_group_element)
    // aggregates_collection[i] is the aggregated shares sent by client i
    // aggregates_collection[0] is never used
    RistScalMat aggregates_collection;

    // to compute the discrete log, we use the baby-step giant-step algorithm
    // for the fixed Ristretto group generator g, we store the table of (x g) for all -2^(b-1) <= x < 2^(b-1),
    // in the computation of discretelog(g, y), we check for a match of y + 2^b a g in the table for every integer a,
    // where b = small_mult_base_table_bit_size = min(weight_bits, MAX_MULT_BASE_TABLE_BIT_SIZE)
    int small_mult_base_table_bit_size; 

    // the table of (x g) for all 0 <= x < 2^small_mult_base_table_bit_size
    std::unordered_map<RistElem, long, RistElemHasher> small_mult_base_table;

    // the aggregated model update in integer format
    // size: (dim)
    std::vector<long> final_update_int;

    // the aggregated model update, converted back into float format
    // size: (dim)
    std::vector<float> final_update_float;

    // the average of the aggregated model update in float format
    // note: average is computed by dividing the aggregated model update by the number of ***honest clients*** (not num_clients)
    // size: (dim)
    std::vector<float> final_update_float_avg;

    // For test purposes only, used in non-private protocols (model updates transferred in clear text) 
    std::vector<std::vector<long>> updates_int_collection;
    std::vector<std::vector<float>> updates_float_collection;
    double float_l2_sq_multiplier;
    // For test purposes only end

    // the check parameters in float format
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

    // import the small multiplication table (a g) for all -2^(small_mult_base_table_bit_size - 1) <= a < 2^(small_mult_base_table_bit_size - 1)
    void import_small_mult_base_table(const std::unordered_map<RistElem, long, RistElemHasher> &small_table) {
        small_mult_base_table = small_table;
    }

    // import the small multiplication table (a g) for all -2^(small_mult_base_table_bit_size - 1) <= a < 2^(small_mult_base_table_bit_size - 1)
    void generate_small_mult_base_table() {
        auto curr = scalar_mult_base(RistScal(-(1 << (small_mult_base_table_bit_size - 1))));
        auto g = scalar_mult_base(c_scal_one);
        for (int i = -(1 << (small_mult_base_table_bit_size - 1));
             i < (1 << (small_mult_base_table_bit_size - 1)); i++) {
            small_mult_base_table.emplace(curr, i);
            add(curr, curr, g);
        }
    }

    // reset the dispute table to all 0 (no malicious flags)
    void reset_dispute_table() {
        dispute_table = std::vector<std::vector<int>>(num_clients + 1, std::vector<int>(num_clients + 1, 0));
    }

    // import committed updates pp which is sent from client i
    void import_committed_updates(std::vector<unsigned char> pp, int i) {
        p3_vec_from_bytes(committed_updates_collection[i], pp);
    }

    // generate the seed that is used to generate public group elements that are used in commitments and proofs
    void generate_aa_seed() {
        rand_init(predicate.aa_seed);
    }

    // ZKP preparation stage, compute the merged group elements and random numbers used to check correctness of commitmets of inner products
    void generate_hh_comb_and_bb_and_aabb();

    // generate the dispute table
    void generate_dispute_table();

    // check the correctness of client i's commitmet of inner products 
    void check_linear_comb_batch_commitments_with_bb(int i);

    // check the bound proof of client i's commitment of inner products
    void check_sq_bound_proof(int i);

    // check che proof of client i's commitments of model updates
    void check_proof(int i);

    // check all the clients' proofs
    void check_proofs();

    // check the clear-text Shamir shares at the dispute stage
    void check_disputes();


    // compute the aggregated model updates
    //  
    // if parallel_on_clients is true, the aggregation is computed parallel on clients, sequential on dim
    // otherwise, parallel on dim, sequential on clients
    // result stored in final_update_float
    //
    // remark: after benchmarking, the cost difference is small, whether parallel_on_clients is true or false
    void compute_final_update(bool parallel_on_clients = false);

    /* ******************************************************************************************************/
    /* abstractions *****************************************************************************************/
    /* see https://github.com/nusdbsystem/risefl/blob/main/doc/server.md#initialize-fl-training for details */
    /* ******************************************************************************************************/

    // used at 1st iteration
    void initialize_from_seed(const std::string &seed) {
        std::vector<unsigned char> seed_bytes = base64_decode(seed);
        predicate.initialize_from_seed(seed_bytes);
        generate_small_mult_base_table();
    }

    void initialize_new_iteration(const CheckParamFloat &check_param_float) {
        if (protocol_type == PROTOCOL_TYPE::NON_PRIV_FLOAT) {
            this->check_param_float = check_param_float;
            generate_aa_seed();
            final_update_float = std::vector<float>(dim, 0);
            normalizing_factor = 1.0;
            float_l2_sq_multiplier = gammas_from_num_samples[predicate.num_samples - 1];
            return;
        }

        auto check_param = set_normalizing_factor_and_compute_check_param_from_gamma(check_param_float);
        reset_server_flags();
        reset_dispute_table();
        set_param(check_param);

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

    // if parallel_on_clients is true, the aggregation is computed parallel on clients, sequential on dim
    // otherwise, parallel on dim, sequential on clients
    // result stored in final_update_float
    //
    // remark: after benchmarking, the cost difference is small, whether parallel_on_clients is true or false
    void finish_iteration(bool parallel_on_clients = false);

    std::string string_api_test(const std::string &a);
};


// compute the discrete log of y, i.e., the number c such that g^c = y
// using small_table as lookup at every step
// algorithm: check for match of (2^b a g) in small_table for every |a| <= per_side_step_count, where 2^b is the size of small_table
long discrete_log(const RistElem &y,
                  const std::unordered_map<RistElem, long, RistElemHasher> &small_table,
                  long per_side_step_count);

#endif //RISEFL_CRYPTO_RISEFL_INTERFACE_SERVER_H
