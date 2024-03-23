//
// Created by yizheng on 15/3/23.
//

// risefl_interface_client.h
// 
// the RiseFL client interface header file
// both server and client interfaces are built on the common interface

#ifndef RISEFL_CRYPTO_RISEFL_INTERFACE_CLIENT_H
#define RISEFL_CRYPTO_RISEFL_INTERFACE_CLIENT_H

#include <vector>
#include "risefl_interface_common.h"
#include "zkp.h"
#include "exception.h"
#include "shamir.h"
#include "ristretto.h"
#include "ristretto_vector.h"
#include "utils.h"
#include "base64.h"
#include "bulletin.h"

// a ClientInterface instance is a CommomInterface instance with some additional members and functions
class ClientInterface : public CommonInterface {
public:

    // the id of the client 
    // NOTE: client_id starts from 1. There is no client 0.
    int client_id;

    // the public keys of the bulletin board. 
    // size: (num_clients + 1)
    // bul_pub_keys[i] is client i's public key, which is used to verify the authenticity of client i's messages
    // bul_pub_keys[0] is never used
    std::vector<SignPubKey> bul_pub_keys;

    // the private key of this client, which is used to sign this client's messages
    SignPrvKey bul_prv_key;

    // weight updates of this client, reshaped into a one-dimensional vector
    // size: (dim)
    std::vector<long> weight_updates;

    // the random scalars used to encrypt weight updates, which will be Shamir shared
    // size: (num_blinds_per_group_element)
    RistScalVec blinds_to_share; 

    // the Diffie–Hellman public key and private key used to establish secret connections with other clients
    PublicKey dh_public_key;
    PrivateKey dh_private_key;

    // the shamir check strings of blinds_to_share
    BatchShamirShareCheckString batch_shamir_share_with_check;

    // shares of the client's secrets blinds_to_share, encrypted with dh_private_key
    // size: (num_clients + 1, num_blinds_per_group_element)
    // encrypted_shamir_shares[i] is the dh_private_key-encrypted share of the blinds_to_share which will sent to client i
    // encrypted_shamir_shares[0] is never used
    std::vector<std::vector<CipherWithNonce>> encrypted_shamir_shares; 

    // the proof that the model update is within a bound
    Proof proof;

    // encrypted shares of other client's secrets
    // size: (num_clients + 1, num_blinds_per_group_element)
    // other_encrypted_shamir_shares[i] is the encrypted share of client i's blinds_to_share. It can be decrypted with dh_public_key_collection[i]
    // other_encrypted_shamir_shares[0] is never used
    std::vector<std::vector<CipherWithNonce>> other_encrypted_shamir_shares; 

    // decrypted shares of other client's secrets
    // size: (num_clients + 1, num_blinds_per_group_element)
    // other_shamir_shares[i] is the decrypted share of client i's blinds_to_share. 
    // other_shamir_shares[0] is never used
    RistScalMat other_shamir_shares; 

    // flags of malicious clients
    // size: (num_clients + 1)
    // flags[i] == 0 means i is not flagged as malicious
    // flags[i] == 1 means i is flagged as malicious
    // flags[0] is never used
    std::vector<int> flags;

    // the list of clients that flags the current client as malicious
    // (if the current client is honest, then everyone on the list must be malicious)
    // size: (num_clients + 1)
    // dispute_clients[i] == 0 means that i does not flags the current client as malicious
    // dispute_clients[i] == 1 means that i flags the current client as malicious
    // dispute_clients[0] is never used
    std::vector<int> dispute_clients;

    // the shares of the blinds sent to the list of clients flags the current client as malicious, in clear text
    // (this number must not exceed max_malicious_clients)
    // size: (num_clients + 1, num_blinds_per_group_element)
    // if dispute_clients[i] == 0, then dispute_shares[i] is the all-zero vector
    // if dispute_clients[i] == 1, then dispute_shares[i] is the clear text of the shares of the blinds sent to client i
    // dispute_clients[0] is the all-zero vector
    RistScalMat dispute_shares;

    // the shares of the blinds received from the list of clients that the current client flags as malicious, in clear text
    // size: (num_clients + 1, num_blinds_per_group_element)
    // if the current client flags i as malicious, other_dispute_shares[i] contains client i's shares in clear text, used for verification of its integrity against client i's Shamir check strings
    RistScalMat other_dispute_shares;

    // the aggregated shares of the blinds held by the current client 
    // the aggregation only sums client[i]'s shares for honest client i
    // size: (num_blinds_per_group_element)
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

    // import weight updates in long int format
    void import_weight_updates(const std::vector<long> &u) {
        weight_updates = u;
    }

    // reset all flags to 0 (i.e. no malicious clients)
    void reset_flags() {
        flags = std::vector<int>(num_clients + 1, 0);
    }

    // generate a Diffie-Hellman key pair, store them in dh_public_key and dh_private_key
    void generate_dh_key_pair() {
        auto key_pair = ::generate_dh_key_pair();
        dh_public_key = key_pair.first;
        dh_private_key = key_pair.second;
    }

    // generate Shamir shares of blinds_to_share and the corresponding check strings
    void generate_batch_shares_and_check_string();

    // encrypt the shares for client i using dh_public_key_collection[i] and dh_private_key, for all 1 <= i <= num_clients
    // store the results in encrypted_shamir_shares
    void encrypt_shares();

    // set the seed which is used to generate public group elements
    void set_aa_seed(const RistHashbytes &s) {
        predicate.aa_seed = s;
    }

    // generate a proof that the weight update is within the region, store the result in the class member proof
    void generate_proof();

    // decrypt all the Shamir shares received from all the other clients, store the result in other_shamir_shares
    void decrypt_shamir_shares();

    // check the integrity of the Shamir shares received from all the other clients, mark flags[i] = 1 if client i fails the check
    void check_shamir_share_integrity();

    // generate the shares in clear text for all i that marks myself as malicious, store the result in dispute_shares
    // if more than max_malicious_clients marks myself as malicious, quit the protocol
    void generate_dispute_shares();

    // for the clients that did not pass the Shamir integrity check, some sent the correct shares, and they are stored in other_shamir_shares
    void update_other_shamir_shares_with_dispute();

    // compute the aggregated share for all the clients marked as honest
    void compute_aggegrated_share();


    /* ******************************************************************************************************/
    /* abstractions *****************************************************************************************/
    /* see https://github.com/nusdbsystem/risefl/blob/main/doc/client.md#initialize-fl-training for details */
    /* ******************************************************************************************************/

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

#endif //RISEFL_CRYPTO_RISEFL_INTERFACE_CLIENT_H
