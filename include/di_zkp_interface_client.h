//
// Created by yizheng on 15/3/23.
//

#ifndef DI_ZKP_CRYPTO_DI_ZKP_INTERFACE_CLIENT_H
#define DI_ZKP_CRYPTO_DI_ZKP_INTERFACE_CLIENT_H

#include <vector>
#include "di_zkp_interface_common.h"
#include "zkp.h"
#include "exception.h"
#include "shamir.h"
#include "ristretto.h"
#include "ristretto_vector.h"
#include "utils.h"
#include "base64.h"
#include "bulletin.h"

class ClientInterface : public CommonInterface {
public:
    int client_id;

    // (num_clients + 1)
    std::vector<SignPubKey> bul_pub_keys;

    SignPrvKey bul_prv_key;

    // (dim)
    std::vector<long> weight_updates;
    // (num_blinds_per_group_element)
    RistScalVec blinds_to_share; // the random scalar used to encrypt weight updates, will be Shamir shared

    PublicKey dh_public_key;
    PrivateKey dh_private_key;

//    // (num_blinds_per_group_element, num_clients + 1)
//    RistScalMat shamir_shares; // shares of the client's secret
//
//    // (num_blinds_per_group_element, max_malicious_clients + 1)
//    std::vector<RistElemP3Vec> shamir_check_string;

    BatchShamirShareCheckString batch_shamir_share_with_check;

    // (num_clients + 1, num_blinds_per_group_element)
    std::vector<std::vector<CipherWithNonce>> encrypted_shamir_shares; // shares of the client's secret

    Proof proof;

    // (num_clients + 1, num_blinds_per_group_element)
    std::vector<std::vector<CipherWithNonce>> other_encrypted_shamir_shares; // encrypted shares of other client's secrets

    // (num_clients + 1, num_blinds_per_group_element)
    RistScalMat other_shamir_shares; // shares of other client's secrets

    // (num_clients + 1)
    std::vector<int> flags;

    // (num_clients + 1)
    std::vector<int> dispute_clients;

    // (num_clients + 1, num_blinds_per_group_element)
    RistScalMat dispute_shares;

    // (num_clients + 1, num_blinds_per_group_element)
    RistScalMat other_dispute_shares;

    // (num_blinds_per_group_element)
    RistScalVec aggregates;

    ClientInterface(int num_clients, int max_malicious_clients,
                    int dim, int num_blinds_per_group_element,
                    int weight_bits, int random_normal_bit_shifter,
                    int num_norm_bound_samples, int inner_prod_bound_bits, int max_bound_sq_bits,
                    CHECK_TYPE check_type,
                    int client_id,
                    const std::vector<SignPubKey> &bul_pub_keys = std::vector<SignPubKey>(),
                    const SignPrvKey &bul_prv_key = SignPrvKey(),
                    bool b_precomp = true,
                    PROTOCOL_TYPE protocol_type = PROTOCOL_TYPE::PRIV) :
            CommonInterface(num_clients, max_malicious_clients,
                            dim, num_blinds_per_group_element,
                            weight_bits, random_normal_bit_shifter,
                            num_norm_bound_samples, inner_prod_bound_bits, max_bound_sq_bits,
                            check_type,
                            b_precomp, protocol_type),
            client_id(client_id),
            bul_pub_keys(bul_pub_keys),
            bul_prv_key(bul_prv_key),
            weight_updates(dim),
            blinds_to_share(num_blinds_per_group_element),
            batch_shamir_share_with_check(num_blinds_per_group_element,
                                          num_clients,
                                          max_malicious_clients + 1),
            encrypted_shamir_shares(num_clients + 1, std::vector<CipherWithNonce>(num_blinds_per_group_element,
                                                                                  CipherWithNonce(RISTSCALBYTES))),
            proof(check_type, num_blinds_per_group_element, num_norm_bound_samples, inner_prod_bound_bits, max_bound_sq_bits),
            other_encrypted_shamir_shares(num_clients + 1, std::vector<CipherWithNonce>(num_blinds_per_group_element,
                                                                                        CipherWithNonce(
                                                                                                RISTSCALBYTES))),
            other_shamir_shares(num_clients + 1, RistScalVec(num_blinds_per_group_element)),
            flags(num_clients + 1),
            dispute_clients(num_clients + 1),
            dispute_shares(num_clients + 1, RistScalVec(num_blinds_per_group_element)),
            other_dispute_shares(num_clients + 1, RistScalVec(num_blinds_per_group_element)),
            aggregates(num_blinds_per_group_element) {}

    void import_weight_updates(const std::vector<long> &u) {
        weight_updates = u;
    }

    void reset_flags() {
        flags = std::vector<int>(num_clients + 1, 0);
    }

    void generate_dh_key_pair() {
        auto key_pair = ::generate_dh_key_pair();
        dh_public_key = key_pair.first;
        dh_private_key = key_pair.second;
    }

    void generate_batch_shares_and_check_string();

    void encrypt_shares();

//    std::vector<unsigned char> generate_encrypted_shares_and_check_string();

    void set_aa_seed(const RistHashbytes &s) {
        predicate.aa_seed = s;
    }

    void generate_proof();

    void decrypt_shamir_shares();

    void check_shamir_share_integrity();

    void generate_dispute_shares();

    void update_other_shamir_shares_with_dispute();

    void compute_aggegrated_share();


    /* abstractions */
    void initialize_from_seed(const std::string &seed) {
        std::vector<unsigned char> seed_bytes = base64_decode(seed);
        predicate.initialize_from_seed(seed_bytes);
    }

    // used at every iteration
    std::vector<unsigned char> send_1_internal(const CheckParam &check_param, const std::vector<long> &u);

    std::vector<unsigned char> send_1_bytes(const CheckParamFloat &check_param_float,
                                            const std::vector<float> &u_float);

    std::string send_1(const CheckParamFloat &check_param_float,
                       const std::vector<float> &u_float);

    std::vector<unsigned char> receive_and_send_2_bytes(const std::vector<unsigned char> &bytes);

    std::string receive_and_send_2(const std::string &bytes_str);

    std::vector<unsigned char> receive_and_send_3_bytes(const std::vector<unsigned char> &bytes);

    std::string receive_and_send_3(const std::string &bytes_str);

    std::vector<unsigned char> receive_and_send_4_bytes(const std::vector<unsigned char> &bytes);

    std::string receive_and_send_4(const std::string &bytes_str);

    std::vector<unsigned char> receive_and_send_5_bytes(const std::vector<unsigned char> &bytes);

    std::string receive_and_send_5(const std::string &bytes_str);

};

#endif //DI_ZKP_CRYPTO_DI_ZKP_INTERFACE_CLIENT_H
