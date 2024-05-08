//
// Created by yizheng on 15/3/23.
//

// risefl_interface_common.h
// 
// the RiseFL common interface header file
// both server and client interfaces are built on the common interface

#ifndef RISEFL_CRYPTO_RISEFL_INTERFACE_COMMON_H
#define RISEFL_CRYPTO_RISEFL_INTERFACE_COMMON_H

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


// Diffie-Hellman public key, C++ wrapper on 
// https://doc.libsodium.org/public-key_cryptography/authenticated_encryption
struct PublicKey {
    // the public key
    std::array<unsigned char, crypto_box_PUBLICKEYBYTES> pub;

    // overload ==
    inline bool operator==(const PublicKey &other) {
        return pub == other.pub;
    }

    // convert the public key to bytes, write the bytes to a vector starting from the place pointed by the iterator it, and move the iterator it to one place after the end
    void export_to_bytestream(std::vector<unsigned char>::iterator &it) const {
        copy_and_shift_dst(pub.begin(), pub.end(), it);
    }

    // read bytes from a vector starting from the place pointed by the iterator it, convert the bytes to the public key, and move the iterator it to one place after the end
    void import_from_bytestream(std::vector<unsigned char>::const_iterator &it_bytes) {
        copy_and_shift_src(it_bytes, pub.begin(), pub.end());
    }
};

// Diffie-Hellman private key, C++ wrapper on 
// https://doc.libsodium.org/public-key_cryptography/authenticated_encryption
struct PrivateKey {
    // the private key
    std::array<unsigned char, crypto_box_SECRETKEYBYTES> prv;

    // overload ==
    inline bool operator==(const PrivateKey &other) {
        return prv == other.prv;
    }
};

// Diffie-Hellman nonce, C++ wrapper on 
// https://doc.libsodium.org/public-key_cryptography/authenticated_encryption
struct Nonce {
    // the nonce
    std::array<unsigned char, crypto_box_NONCEBYTES> non;

    // overload ==
    inline bool operator==(const Nonce &other) {
        return non == other.non;
    }
};


// Diffie-Hellman cipher with nonce, i.e. the encrypted message
struct CipherWithNonce {
    // the cipher
    std::vector<unsigned char> cip;
    // the nonce
    Nonce non;

    // create an empty encrypted message from the length of the clear text
    CipherWithNonce(int clear_text_length) : cip(clear_text_length + crypto_box_MACBYTES) {}


    // overload ==
    inline bool operator==(const CipherWithNonce &other) {
        return cip == other.cip && non == other.non;
    }

    // convert the encrypted message to bytes, write the bytes to a vector starting from the place pointed by the iterator it, and move the iterator it to one place after the end
    void export_to_bytestream(std::vector<unsigned char>::iterator &it) const {
        copy_and_shift_dst(cip.begin(), cip.end(), it);
        copy_and_shift_dst(non.non.begin(), non.non.end(), it);
    }

    // read bytes from a vector starting from the place pointed by the iterator it, convert the bytes to the encrypted message, and move the iterator it to one place after the end
    void import_from_bytestream(std::vector<unsigned char>::const_iterator &it_bytes) {
        copy_and_shift_src(it_bytes, cip.begin(), cip.end());
        copy_and_shift_src(it_bytes, non.non.begin(), non.non.end());
    }
};


// generate a Diffie-Hellman key pair
inline std::pair<PublicKey, PrivateKey> generate_dh_key_pair() {
    PublicKey pub;
    PrivateKey prv;
    crypto_box_keypair(pub.pub.data(), prv.prv.data());
    return {pub, prv};
}

// helper function to encrypt a message 
// encrypt a message mes to a cipher text using nonce non, public key pub and private key prv, and write the cipher to cip
inline int encrypt(std::vector<unsigned char> &cip,
                   const std::vector<unsigned char> &mes,
                   const Nonce &non,
                   const PublicKey &pub,
                   const PrivateKey &prv) {
    cip.resize(crypto_box_MACBYTES + mes.size());
    return crypto_box_easy(cip.data(), mes.data(), mes.size(), non.non.data(),
                           pub.pub.data(), prv.prv.data());
}

// generate a random nonce
inline Nonce generate_random_nonce() {
    Nonce non;
    randombytes_buf(non.non.data(), crypto_box_NONCEBYTES);
    return non;
}

// encrypt a message mes using public key pub and private key prv, output the encrypted message
inline CipherWithNonce encrypt(const std::vector<unsigned char> &mes,
                               const PublicKey &pub,
                               const PrivateKey &prv) {
    CipherWithNonce ret(mes.size());
    ret.non = generate_random_nonce();
    encrypt(ret.cip, mes, ret.non, pub, prv);
    return ret;
}

// helper function to decrypt a cipher cip using nonce non, public key pub and private key prv
// if decryption successful, write the clear text to mes and return 0
// if decryption unsuccessful, return -1 
inline int decrypt(std::vector<unsigned char> &mes,
                   const std::vector<unsigned char> &cip,
                   const Nonce &non,
                   const PublicKey &pub,
                   const PrivateKey &prv) {
    mes.resize(cip.size() - crypto_box_MACBYTES);
    return crypto_box_open_easy(mes.data(), cip.data(), cip.size(), non.non.data(),
                                pub.pub.data(), prv.prv.data());
}

// decrypt a message cip_with_nonce using public key pub and private key prv, outputs the decrypted message
// if decryption fails, throw an exception DecryptionError
inline std::vector<unsigned char> decrypt(const CipherWithNonce &cip_with_nonce,
                                          const PublicKey &pub,
                                          const PrivateKey &prv) {
    std::vector<unsigned char> mes;
    if (decrypt(mes, cip_with_nonce.cip, cip_with_nonce.non, pub, prv) == -1)
        throw DecryptionError();
    return mes;
}

// the check type: L2 norm, sphere, cosine similarity
enum class CHECK_TYPE {
    L2NORM, SPHERE, COSINE_SIM
};

// the parameters of L2 norm check
struct L2NormCheckParam {
    // the square of the L2 norm bound
    RistScal bound_sq;
};

// the parameters of the sphere check
struct SphereCheckParam {
    // the square of the radius of the sphere
    RistScal bound_sq;
    // the center of the sphere
    RistScalVec center;
};

// the parameters of the cosine similarity check
// checking: ||u|| <= r AND cosine(u, v) >= alpha
struct CosineSimCheckParam {
    // r^2
    RistScal bound_sq;
    // v
    RistScalVec pivot;

    // for the next two parameters, we actually check l2_sq_multiplier * ||u||^2 ||v||^2 <= inner_prod_sq_multiplier * <u,v>^2
    RistScal l2_sq_multiplier;
    RistScal inner_prod_sq_multiplier;
};

// the check parameters
struct CheckParam {
    CHECK_TYPE check_type;
    L2NormCheckParam l2_param;
    SphereCheckParam sphere_param;
    CosineSimCheckParam cosine_param;

    CheckParam(CHECK_TYPE check_type) : check_type(check_type) {}
};

// if CHECK_TYPE is COSINE_SIM, return 1 because ZKP needs to compute the commitment of an extra inner product <u,v>
// otherwise, return 0
int extra_inner_prod(CHECK_TYPE c);

// if CHECK_TYPE is COSINE_SIM, ZKP needs to compute the commitment of an extra square ||u||^2
// otherwise, return 0
int extra_square(CHECK_TYPE c);


// the check predicate
struct Predicate {
public:
    // number of model parameters
    int dim;
    // how many public group elements are used to commit model updates (the paper uses dim)
    int num_weight_keys;
    // how many blinds are used per public group element that is used to commit model updates (the paper uses 1)
    int num_blinds_per_group_element;

    // the number of normal samples k used in the probablistic check
    int num_samples;

    // the random normal samples a_1, ..., a_k are obtained by multiplying random samples from Gaussian(0,1)
    // by 2^random_normal_bit_shifter and rounded to the nearest integer
    int random_normal_bit_shifter; 

    // the bit-bound of linear combinations <a_t, u> of random normal samples and model updates obtained by sampling
    //  i.e. every <a_t, u> is in the interval [-2^(inner_prod_bound_bits-1), 2^(inner_prod_bound_bits-1)) 
    int inner_prod_bound_bits;

    // the check parameters
    CheckParam check_param;

    // the square of sum of ||<a_t,u>||^2 is guraranteed to be bounded by max_bound_sq_bits bits
    //  i.e. sum of ||<a_t,u>||^2 is in the interval [-2^(max_bound_sq_bits-1), 2^(max_bound_sq_bits-1)) 
    int max_bound_sq_bits;

    // the 1st set of group elements used in the proof of bound of the inner product <a_t,u>, 
    // size: (num_samples * max_bound_sq_bits)
    RistP3VecAndBytes bound_elem_keys_1;

    // the 2nd set of group elements used in the proof of bound of the inner product <a_t,u>, 
    // size: (num_samples * max_bound_sq_bits)
    RistP3VecAndBytes bound_elem_keys_2;

    // the seed used to generate the public group elements
    RistHashbytes aa_seed;

    // the group elements used to encrypt all the clients' weight updates
    // size: (num_weight_keys)
    RistElemP3Vec hh; 

    // whether precomputation is used in commitment (the paper uses false)
    bool b_precomp;

    // if precomputation is used in commitment, the precomputation table of hh
    // size: (b_precomp ? num_weight_keys : 0)
    std::vector<RistElemPrecompTable> hh_precomp;

    // the linear combination of hh by random normal samples in blocks of num_weight_keys
    // size: if check_param.check_type == CHECK_TYPE::COSINE_SIM, then (num_samples + 2, num_blinds_per_group_element) 
    //  otherwise, (num_samples + 1, num_blinds_per_group_element) 
    RistP3MatAndBytes hh_comb; 

    // the group element used in proof of square     
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
    // initialize all the group elements using seed
    //      the group elements include: hh, bound_elem_keys_1, bound_elem_keys_2, square_key
    void initialize_from_seed(const std::vector<unsigned char> &seed);
};


// the ZKP proof that checks the predicate
struct Proof {
    // the type of the check
    CHECK_TYPE check_type;

    // the commitments of inner products 
    // size: (num_norm_bound_samples + 1) if not cosine similarity check, (num_norm_bound_samples + 2) if cosine similarity check, 
    RistP3VecAndBytes linear_comb_batch_commitments;

    // another set of commitments of inner products using square_key, used to prove square
    // size: (num_norm_bound_samples + 1) if not cosine similarity check, (num_norm_bound_samples + 2) if cosine similarity check, 
    RistP3VecAndBytes linear_comb_single_commitments;

    // proof that the two sets of commitments of inner products are from the same values
    // size: (num_blinds_per_group_element, num_norm_bound_samples + 1), if not cosine similarity check, (num_blinds_per_group_element, num_norm_bound_samples + 2) if cosine similarity check, 
    PedersenBatchWithZeroProof proof_well_formed;

    // proof that every inner product is bounded (this is necessary to prevent overflow in finite field arithmetic)
    // size: (num_norm_bound_samples * inner_prod_bound_bits)
    RangeProofPowerTwoP3 proof_linear_comb_bound;

    // the commitments of squares of inner products
    // size: (num_norm_bound_samples) if not cosine similarity check, (num_norm_bound_samples + 1) if cosine similarity check, 
    RistP3VecAndBytes square_commitments;

    // the proof of squares of inner products, based on linear_comb_single_commitments
    BatchPedersenWithSquareProof proof_squares;

    // the proof that the sum of squares of inner products is bounded
    RangeProofPowerTwoP3 proof_sum_range;

    // if cosine similarity, an extra bound proof
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

    // the size of the proof in bytes
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

    // convert the proof to bytes, write the bytes to a vector starting from the place pointed by the iterator it, and move the iterator it to one place after the end
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

    // read bytes from a vector starting from the place pointed by the iterator it, convert the bytes to the proof, and move the iterator it to one place after the end
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

// For test purposes only, using non-private protocols (model updates transferred in clear text) 
enum class PROTOCOL_TYPE {
    PRIV, NON_PRIV_INT, NON_PRIV_FLOAT
};

// the parameters of L2 norm check in float format, prior to rounding to int
struct L2NormCheckParamFloat {
    // the L2 norm bound
    float bound;
};

// the parameters of sphere check in float format, prior to rounding to int
struct SphereCheckParamFloat {
    // the L2 norm bound
    float bound;
    // the center of the sphere
    std::vector<float> center;
};

// the parameters of cosine similarity check in float format, prior to rounding to int
// checking: ||u|| <= r AND cosine(u, v) >= alpha
struct CosineSimCheckParamFloat {
    // r
    float bound;
    // v
    std::vector<float> pivot;
    // alpha
    float cosine_bound;
};

// the check parameters  in float format, prior to rounding to int
struct CheckParamFloat {
    CHECK_TYPE check_type;
    L2NormCheckParamFloat l2_param;
    SphereCheckParamFloat sphere_param;
    CosineSimCheckParamFloat cosine_param;

    CheckParamFloat(CHECK_TYPE check_type) : check_type(check_type) {}
};



// the common interface class used in both server and client
// used as the base class for ServerInterface and ClientInterface
class CommonInterface {
public:
    // For test purposes only, possibly using non-private protocols (model updates transferred in clear text)
    PROTOCOL_TYPE protocol_type;

    // every coordinate of the model update is divided by normalizing_factor, resulting in a number in [-1, 1]
    float normalizing_factor;

    // the total number of clients
    int num_clients;

    // the maximum number of malicious clients
    int max_malicious_clients;

    // the number of model parameters
    int dim;

    // the bit-length of weight updates
    int weight_bits;

    // the server's flags on malicious clients
    // size: (num_clients + 1)
    // for i > 0, server_flag[i] == 1 means client i is flagged as malicious by the server, 0 otherwise
    std::vector<int> server_flags;

    // the collection of Diffe-Hellman public keys
    // size: (num_clients + 1)
    // for i > 0, dh_public_key_collection[i] is client i's Diffe-Hellman public key
    std::vector<PublicKey> dh_public_key_collection;

    // the collection of public keys on the bulletin board that are used to verify authenticity of messages
    // size: (num_clients + 1)
    // for i > 0, bul_signed_pub_key_collection[i] is client i's public key on the bulletin board
    std::vector<std::vector<unsigned char>> bul_signed_pub_key_collection;

    // the check predicate
    Predicate predicate;

    // the collection of Shamir check strings
    // size: (num_clients + 1, num_blinds_per_group_element, max_malicious_clients + 1)
    // for i > 0, shamir_check_string_collection[i] is client i's Shamir check strings
    std::vector<RistP3MatAndBytes> shamir_check_string_collection;

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

    // set the check parameter
    void set_param(const CheckParam &c) {
        predicate.check_param = c;
    }
    
    // from the check parameters (in float format), set the normalizing factor and compute the internal check parameters (in integer format)
    CheckParam set_normalizing_factor_and_compute_check_param_from_gamma(CheckParamFloat check_param_float);

    // reset the server flags to all 0 (i.e. no malicious clients)
    void reset_server_flags() {
        server_flags = std::vector<int>(num_clients + 1, 0);
    }

    // get the number of valid clients
    int valid_client_count() {
        int ret = 0;
        for (int i = 1; i <= num_clients; i++) {
            ret += (1 - server_flags[i]);
        }
        return ret;
    }
};

// read bytes from a vector starting from the place pointed by the iterator it, convert the bytes to flags[1]..flags[num_clients], and move the iterator it to one place after the end
// flags[0] never changed
void import_flags_start_from_1_from_bytestream(std::vector<int> &flags, std::vector<unsigned char>::const_iterator &it);

// convert flags[1]..flags[num_clients] to bytes, write the bytes to a vector starting from the place pointed by the iterator it, and move the iterator it to one place after the end
// flags[0] never exported
void export_flags_start_from_1_to_bytestream(const std::vector<int> &flags, std::vector<unsigned char>::iterator &it);

/**
* This method converts the SignPubKey object to base64 encoded string that
* can be passed via SwigPyObject and used in python
*/
std::string convert_sign_pub_key_to_string(const SignPubKey &pub_key);

/**
* This method converts the SignPrvKey object to base64 encoded string that
* can be passed via SwigPyObject and used in python
*/
std::string convert_sign_prv_key_to_string(const SignPrvKey &prv_key);

/**
* This method converts the base64 string passed from SwigPyObject to SignPubKey
*/
SignPubKey convert_string_to_sign_pub_key(const std::string &key_str);

/**
* This method converts the base64 string passed from SwigPyObject to SignPrvKey
*/
SignPrvKey convert_string_to_sign_prv_key(const std::string &key_str);

#endif //RISEFL_CRYPTO_RISEFL_INTERFACE_COMMON_H
