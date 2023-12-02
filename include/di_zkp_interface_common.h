//
// Created by yizheng on 15/3/23.
//

#ifndef DI_ZKP_CRYPTO_DI_ZKP_INTERFACE_COMMON_H
#define DI_ZKP_CRYPTO_DI_ZKP_INTERFACE_COMMON_H

#include <vector>
#include <array>
#include "ristretto.h"
#include "ristretto_vector.h"
#include "zkp.h"
#include "rist_fast_computation.h"
#include "exception.h"
#include "utils.h"
#include "bulletin.h"

constexpr int ENCRISTSCALBYTES = crypto_box_MACBYTES + RISTSCALBYTES + crypto_box_NONCEBYTES;
static const unsigned int NUM_THREADS = std::thread::hardware_concurrency();


// encryption wrappers
// https://doc.libsodium.org/public-key_cryptography/authenticated_encryption
struct PublicKey {
    std::array<unsigned char, crypto_box_PUBLICKEYBYTES> pub;

    inline bool operator==(const PublicKey &other) {
        return pub == other.pub;
    }

    void export_to_bytestream(std::vector<unsigned char>::iterator &it) const {
        copy_and_shift_dst(pub.begin(), pub.end(), it);
    }

    void import_from_bytestream(std::vector<unsigned char>::const_iterator &it_bytes) {
        copy_and_shift_src(it_bytes, pub.begin(), pub.end());
    }
};

struct PrivateKey {
    std::array<unsigned char, crypto_box_SECRETKEYBYTES> prv;

    inline bool operator==(const PrivateKey &other) {
        return prv == other.prv;
    }
};

struct Nonce {
    std::array<unsigned char, crypto_box_NONCEBYTES> non;

    inline bool operator==(const Nonce &other) {
        return non == other.non;
    }
};

struct CipherWithNonce {
    std::vector<unsigned char> cip;
    Nonce non;

    CipherWithNonce(int clear_text_length) : cip(clear_text_length + crypto_box_MACBYTES) {}

    inline bool operator==(const CipherWithNonce &other) {
        return cip == other.cip && non == other.non;
    }

    void export_to_bytestream(std::vector<unsigned char>::iterator &it) const {
        copy_and_shift_dst(cip.begin(), cip.end(), it);
        copy_and_shift_dst(non.non.begin(), non.non.end(), it);
    }

    void import_from_bytestream(std::vector<unsigned char>::const_iterator &it_bytes) {
        copy_and_shift_src(it_bytes, cip.begin(), cip.end());
        copy_and_shift_src(it_bytes, non.non.begin(), non.non.end());
    }
};

//void bytestream_from_cipher_with_nonce(std::vector<unsigned char>::iterator &it_bytes, const CipherWithNonce &r);
//
//void cipher_with_nonce_from_bytestream(CipherWithNonce &r, int cipher_length,
//                                       std::vector<unsigned char>::const_iterator &it_bytes);


// generates a Diffie-Hellman key pair
inline std::pair<PublicKey, PrivateKey> generate_dh_key_pair() {
    PublicKey pub;
    PrivateKey prv;
    crypto_box_keypair(pub.pub.data(), prv.prv.data());
    return {pub, prv};
}

inline int encrypt(std::vector<unsigned char> &cip,
                   const std::vector<unsigned char> &mes,
                   const Nonce &non,
                   const PublicKey &pub,
                   const PrivateKey &prv) {
    cip.resize(crypto_box_MACBYTES + mes.size());
    return crypto_box_easy(cip.data(), mes.data(), mes.size(), non.non.data(),
                           pub.pub.data(), prv.prv.data());
}

inline Nonce generate_random_nonce() {
    Nonce non;
    randombytes_buf(non.non.data(), crypto_box_NONCEBYTES);
    return non;
}

inline CipherWithNonce encrypt(const std::vector<unsigned char> &mes,
                               const PublicKey &pub,
                               const PrivateKey &prv) {
    CipherWithNonce ret(mes.size());
    ret.non = generate_random_nonce();
    encrypt(ret.cip, mes, ret.non, pub, prv);
    return ret;
}

// returns 0 if decrypt success, -1 if failure
inline int decrypt(std::vector<unsigned char> &mes,
                   const std::vector<unsigned char> &cip,
                   const Nonce &non,
                   const PublicKey &pub,
                   const PrivateKey &prv) {
    mes.resize(cip.size() - crypto_box_MACBYTES);
    return crypto_box_open_easy(mes.data(), cip.data(), cip.size(), non.non.data(),
                                pub.pub.data(), prv.prv.data());
}

inline std::vector<unsigned char> decrypt(const CipherWithNonce &cip_with_nonce,
                                          const PublicKey &pub,
                                          const PrivateKey &prv) {
    std::vector<unsigned char> mes;
    if (decrypt(mes, cip_with_nonce.cip, cip_with_nonce.non, pub, prv) == -1)
        throw DecryptionError();
    return mes;
}

enum class CHECK_TYPE {
    L2NORM, SPHERE, COSINE_SIM, ZENO
};

struct L2NormCheckParam {
    RistScal bound_sq;
};

struct SphereCheckParam {
    RistScal bound_sq;
    RistScalVec center;
};

struct CosineSimCheckParam {
    RistScal bound_sq;
    RistScalVec pivot;
    RistScal l2_sq_multiplier;
    RistScal inner_prod_sq_multiplier;
};

struct CheckParam {
    CHECK_TYPE check_type;
    L2NormCheckParam l2_param;
    SphereCheckParam sphere_param;
    CosineSimCheckParam cosine_param;

    CheckParam(CHECK_TYPE check_type) : check_type(check_type) {}
};

int extra_inner_prod(CHECK_TYPE c);
int extra_square(CHECK_TYPE c);

struct Predicate {
public:
    int dim;
    int num_weight_keys;
    int num_blinds_per_group_element;

    int num_samples;

    int random_normal_bit_shifter; // a_1, ..., a_m are obtained by multiplying coeff of random samples from Gaussian(0,1)
    // by 2^random_normal_bit_shifter and round to nearest integer

    // the bound of linear combinations obtained by sampling
    int inner_prod_bound_bits;

//    // the type of the check
//    CHECK_TYPE check_type;

    CheckParam check_param;
//    /* l2 norm & cosine similarity:  the bound of the square of the l2 norm */
//    /* sphere defense: the bound of the square of the radius of the sphere */
//    RistScal bound_sq;

    // the square of bound of l2 norm is bounded by max_bound_sq_bits bits
    int max_bound_sq_bits;

    // (num_samples * max_bound_sq_bits)
    RistP3VecAndBytes bound_elem_keys_1;

    // (num_samples * max_bound_sq_bits)
    RistP3VecAndBytes bound_elem_keys_2;

    RistHashbytes aa_seed;

    // (num_weight_keys)
    RistElemP3Vec hh; // the group elements used to encrypt all the clients' weight updates

    bool b_precomp;

    // (b_precomp ? num_weight_keys : 0)
    std::vector<RistElemPrecompTable> hh_precomp; // dimension num_weight_keys, the precomputation table of hh

    // (num_samples + 1 or +2, num_blinds_per_group_element)
    RistP3MatAndBytes hh_comb; // the linear combination of hh by (aa_0, aa_pos) in blocks of num_weight_keys

    RistP3AndBytes square_key;

    Predicate(int dim, int num_blinds_per_group_element,
              int num_samples,
              int random_normal_bit_shifter,
              CHECK_TYPE check_type,
//              CheckParam check_param,
              int inner_prod_bound_bits,
              int max_bound_sq_bits,
              bool b_precomp = true) :
            dim(dim),
            num_weight_keys((dim % num_blinds_per_group_element == 0) ? dim / num_blinds_per_group_element :
                            dim / num_blinds_per_group_element + 1),
            num_blinds_per_group_element(num_blinds_per_group_element),
            random_normal_bit_shifter(random_normal_bit_shifter),
            num_samples(num_samples),
            inner_prod_bound_bits(inner_prod_bound_bits),
            check_param(check_type),
            max_bound_sq_bits(max_bound_sq_bits),
            bound_elem_keys_1(num_samples * max_bound_sq_bits),
            bound_elem_keys_2(num_samples * max_bound_sq_bits),
            hh(num_weight_keys),
            b_precomp(b_precomp),
            hh_precomp(b_precomp ? num_weight_keys : 0),
            hh_comb(num_samples + 1 + extra_inner_prod(check_type),
                    num_blinds_per_group_element) {
//        switch (check_type) {
//            case CHECK_TYPE::L2NORM: this->check_param.l2_param = check_param.l2_param; break;
//            case CHECK_TYPE::SPHERE: this->check_param.sphere_param = check_param.sphere_param; break;
//            case CHECK_TYPE::COSINE_SIM: this->check_param.cosine_param = check_param.cosine_param; break;
//            case CHECK_TYPE::ZENO: this->check_param.zeno_param = check_param.zeno_param; break;
//        }
    }

//    void set_l2_bound_sq(const RistScal &B) {
//        check_param.l2_param.bound_sq = B;
//    }

//    // generate aa from a random seed
//    void generate_aa_from_seed(const RistHashbytes &seed);

//    // only server needs to compute hh_comb
//    void compute_hh_comb();
//
//    // only client needs to check hh_comb
//    bool check_correctness_hh_comb();

//    // both server and client need to compute hh_comb_nonzero_sum
//    void compute_hh_comb_nonzero_sum() {
////        sum_single_thread(hh_comb_nonzero_sum, RistElemP3Vec(hh_comb.begin() + 1, hh_comb.end()));
//        hh_comb_nonzero_sum = sum_single_thread(RistElemP3Vec(hh_comb.begin() + 1, hh_comb.end()));
//    }

    /* abstractions */
    void initialize_from_seed(const std::vector<unsigned char> &seed);
};


struct Proof {
    CHECK_TYPE check_type;

    // (num_norm_bound_samples + 1)
    RistP3VecAndBytes linear_comb_batch_commitments;

    // (num_norm_bound_samples + 1)
    RistP3VecAndBytes linear_comb_single_commitments;

    // (num_blinds_per_group_element, num_norm_bound_samples + 1)
    PedersenBatchWithZeroProof proof_well_formed;

    // (num_norm_bound_samples * inner_prod_bound_bits)
//    std::vector<RangeProofPowerTwoP3> proof_linear_comb_bound;
    RangeProofPowerTwoP3 proof_linear_comb_bound;

    // (num_norm_bound_samples)
    RistP3VecAndBytes square_commitments;

    // (num_norm_bound_samples), based on linear_comb_single_commitments
    BatchPedersenWithSquareProof proof_squares;

    RangeProofPowerTwoP3 proof_sum_range;

    // cosine similarity needs an extra proof for l2 bound
    RangeProofPowerTwoP3 proof_sum_range_for_cosine;

    Proof(CHECK_TYPE check_type, int num_blinds_per_group_element, int num_norm_bound_samples, int inner_prod_bound_bits, int sum_bound_bits) :
            check_type(check_type),
            linear_comb_batch_commitments(num_norm_bound_samples + 1 + extra_inner_prod(check_type)),
            linear_comb_single_commitments(num_norm_bound_samples + 1 + extra_inner_prod(check_type)),
            proof_well_formed(num_blinds_per_group_element, num_norm_bound_samples + 1 + extra_inner_prod(check_type)),
            proof_linear_comb_bound(num_norm_bound_samples * inner_prod_bound_bits),
            square_commitments(num_norm_bound_samples + extra_square(check_type)),
            proof_squares(num_norm_bound_samples + extra_square(check_type)),
            proof_sum_range(sum_bound_bits),
            proof_sum_range_for_cosine(check_type == CHECK_TYPE::COSINE_SIM ?
                sum_bound_bits :
                0){}


//    void set_sum_bound_bits(int sum_bound_bits) {
//        proof_sum_range.set_bound_bits(sum_bound_bits);
//    }

    static int
    size(CHECK_TYPE check_type,
         int num_blinds_per_group_element, int num_norm_bound_samples, int inner_prod_bound_bits,
         int sq_sum_bound_bits) {
        int num_norm_bound_samples_extra = num_norm_bound_samples + extra_inner_prod(check_type);
        int num_sq_commitments = num_norm_bound_samples + extra_square(check_type);
        return RISTBYTES * (num_norm_bound_samples_extra + 1) * 2 + // linear comb batch commitments + single commitments
                PedersenBatchWithZeroProof::size(num_blinds_per_group_element, num_norm_bound_samples_extra + 1) +
                // proof well formed
                RangeProofPowerTwoP3::size(num_norm_bound_samples * inner_prod_bound_bits) +
                // proof abs val of linear comb bound
                RISTBYTES * num_sq_commitments + // square commitments
                BatchPedersenWithSquareProof::size(num_sq_commitments) + // proof squares
                RangeProofPowerTwoP3::size(sq_sum_bound_bits) + // proof sum range
                RangeProofPowerTwoP3::size(sq_sum_bound_bits) * (check_type == CHECK_TYPE::COSINE_SIM); // proof sum range for cosine
    }

    void export_to_bytestream(std::vector<unsigned char>::iterator &it) const {
        linear_comb_batch_commitments.export_to_bytestream(it);
        linear_comb_single_commitments.export_to_bytestream(it);
        proof_well_formed.export_to_bytestream(it);
        proof_linear_comb_bound.export_to_bytestream(it);
        square_commitments.export_to_bytestream(it);
        proof_squares.export_to_bytestream(it);
        proof_sum_range.export_to_bytestream(it);
        if (check_type == CHECK_TYPE::COSINE_SIM) {
            proof_sum_range_for_cosine.export_to_bytestream(it);
        }
    }

    void import_from_bytestream(std::vector<unsigned char>::const_iterator &it) {
        linear_comb_batch_commitments.import_from_bytestream(it);
        linear_comb_single_commitments.import_from_bytestream(it);
        proof_well_formed.import_from_bytestream(it);
        proof_linear_comb_bound.import_from_bytestream(it);
        square_commitments.import_from_bytestream(it);
        proof_squares.import_from_bytestream(it);
        proof_sum_range.import_from_bytestream(it);
        if (check_type == CHECK_TYPE::COSINE_SIM) {
            proof_sum_range_for_cosine.import_from_bytestream(it);
        }
    }
};

enum class PROTOCOL_TYPE {
    PRIV, NON_PRIV_INT, NON_PRIV_FLOAT
};

struct L2NormCheckParamFloat {
    float bound;
};

struct SphereCheckParamFloat {
    float bound;
    std::vector<float> center;
};

struct CosineSimCheckParamFloat {
    float bound;
    std::vector<float> pivot;
    float cosine_bound;
};

struct CheckParamFloat {
    CHECK_TYPE check_type;
    L2NormCheckParamFloat l2_param;
    SphereCheckParamFloat sphere_param;
    CosineSimCheckParamFloat cosine_param;

    CheckParamFloat(CHECK_TYPE check_type) : check_type(check_type) {}
};

class CommonInterface {
public:
    PROTOCOL_TYPE protocol_type;

    float normalizing_factor;

    int num_clients;
    int max_malicious_clients;
    int dim;
    int weight_bits;

    // (num_clients + 1)
    std::vector<int> server_flags;

    // (num_clients + 1)
    std::vector<PublicKey> dh_public_key_collection;

    // (num_clients + 1)
    std::vector<std::vector<unsigned char>> bul_signed_pub_key_collection;

//    RistHashbytes aa_seed;

    Predicate predicate;
//
    // (num_clients + 1, num_blinds_per_group_element, max_malicious_clients + 1)
    std::vector<RistP3MatAndBytes> shamir_check_string_collection;
//
//    // (num_clients + 1, num_blinds_per_group_element)
//    RistScalMat shamir_batch_coeff_collection;

    CommonInterface(int num_clients, int max_malicious_clients,
                    int dim, int num_blinds_per_group_element,
                    int weight_bits, int random_normal_bit_shifter,
                    int num_samples, int inner_prod_bound_bits, int max_bound_sq_bits,
                    CHECK_TYPE check_type,
                    bool b_precomp = true,
                    PROTOCOL_TYPE protocol_type = PROTOCOL_TYPE::PRIV) :
            protocol_type(protocol_type),
            num_clients(num_clients),
            max_malicious_clients(max_malicious_clients),
            dim(dim),
            weight_bits(weight_bits),
            server_flags(num_clients + 1, 0),
            dh_public_key_collection(num_clients + 1),
            bul_signed_pub_key_collection(num_clients + 1,
                                          std::vector<unsigned char>(SIGNBYTES + crypto_box_PUBLICKEYBYTES)),
            predicate(dim, num_blinds_per_group_element, num_samples, random_normal_bit_shifter,
                      check_type,
                      inner_prod_bound_bits,
                      max_bound_sq_bits, b_precomp),
            shamir_check_string_collection(num_clients + 1,
                                           RistP3MatAndBytes(num_blinds_per_group_element,
                                                             max_malicious_clients + 1)) {}

//    void import_group_element_keys(const RistElemP3Vec &hh, const RistElemP3Vec &k1, const RistElemP3Vec &k2) {
//        predicate.hh = hh;
//        predicate.bound_elem_keys_1 = k1;
//        predicate.bound_elem_keys_2 = k2;
//    }

    void set_param(const CheckParam &c) {
        predicate.check_param = c;
    }

//    RistScal set_normalizing_factor_and_compute_bound_sq(float norm_bound, float standard_deviation_factor);
    CheckParam set_normalizing_factor_and_compute_check_param_from_gamma(CheckParamFloat check_param_float);

    void reset_server_flags() {
        server_flags = std::vector<int>(num_clients + 1, 0);
    }

    int valid_client_count() {
        int ret = 0;
        for (int i = 1; i <= num_clients; i++) {
            ret += (1 - server_flags[i]);
        }
        return ret;
    }

//    void generate_aa();

//    void compute_hh_comb_nonzero_sum() {
//        predicate.compute_hh_comb_nonzero_sum();
//    }
};

void import_flags_start_from_1_from_bytestream(std::vector<int> &flags, std::vector<unsigned char>::const_iterator &it);

void export_flags_start_from_1_to_bytestream(const std::vector<int> &flags, std::vector<unsigned char>::iterator &it);

#endif //DI_ZKP_CRYPTO_DI_ZKP_INTERFACE_COMMON_H
