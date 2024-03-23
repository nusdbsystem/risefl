//
// Created by yizheng on 8/3/23.
//

// shamir.h
//
// verifiable Shamir sharing tools

#ifndef RISEFL_CRYPTO_SHAMIR_H
#define RISEFL_CRYPTO_SHAMIR_H

#include <chrono>
#include <cassert>
#include <vector>

#include "ristretto.h"
#include "ristretto_vector.h"
#include "zkp.h"
#include "rist_fast_computation.h"

// compute the shamir shares of s with n clients, any t of which can recover s
// also compute the check string
// returns: {rr, zz} where :
//      rr is a vector of Ristretto scalars of size (n+1), rr[0] is never used, rr[i] is client i's share of s
//      zz is a vector of Ristretto group elements of size t
std::pair<RistScalVec, RistElemP3Vec> compute_shamir_shares_with_check_string(const RistScal &s, int n, int t);

// The shamir shares of a batch of scalars and the corresponding check strings
class BatchShamirShareCheckString {
public:
    // size: (batch_size, n + 1), where n is the number of clients
    // shares[j][0] is never used
    // shares[j][i] is client i's share of the j-th element
    std::vector<RistScalVec> shares;

    // size:(batch_size, t), where t is the number required to recover the shares
    // check_strings[j] is the check string of the j-th element
    RistP3MatAndBytes check_strings;

//    std::vector<RistElemP3Vec> check_strings;

//     (batch_size, t * RISTBYTES)
//    std::vector<std::vector<unsigned char>> check_strings_bytes;

//    // (batch_size)
//    RistScalVec check_string_coeffs;
//
//    // (t)
//    RistP3VecAndBytes batch_check_string;

    // (t * RISTBYTES)
//    std::vector<unsigned char> batch_check_string_bytes;

    BatchShamirShareCheckString() = default;

    BatchShamirShareCheckString(int batch_size, int n, int t) :
            shares(batch_size, RistScalVec(n + 1)),
            check_strings(batch_size, t)
//        check_strings_bytes(batch_size, std::vector<unsigned char>(t * RISTBYTES)),
//        check_string_coeffs(batch_size),
//        batch_check_string(t)
//        batch_check_string_bytes(t * RISTBYTES)
    {}
};

// compute the shamir shares of a vector ss with n clients, any t of which can recover ss
// also compute the check string
// returns: an BatchShamirShareCheckString instance containing the shares and the check strings
BatchShamirShareCheckString compute_batch_shamir_shares_with_check_string(const RistScalVec &ss, int n, int t);

//inline std::pair<RistScalVec, RistElemP3Vec> compute_shamir_shares_with_check_string(const RistScal &s, int n, int t) {
//    return compute_shamir_shares_with_check_string(ZZ_p_from_rist(s), n, t);
//}

// verifies the integrity of the shamir share with the check, given the number of client t that is required to recover the secret and the current client id i
// returns: true if valid, false if invalid
bool b_shamir_verify(const RistScal &share, const RistElemP3Vec &checks, int t, int i);

// verifies the integrity of the vector of shamir shares with the check, given the number of client t that is required to recover the secret and the current client id i
// returns: true if valid, false if invalid
bool b_shamir_verify(const RistScalVec &shares, const RistElemP3Mat &checks, int t, int i);

//RistScalVec rist_from_vec_zz_p(const vec_ZZ_p &vv);
//
//vec_ZZ_p vec_zz_p_from_rist(const RistScalVec &ss);

//ZZ_p shamir_recover(const vec_ZZ_p &shares, const std::vector<int> &indices, int t);

// recovers the shamir share from t shares and t indices
// indices must have length at least t
// shares must have length at least t
// shares[i] is client indices[i]'s share
RistScal shamir_recover(const RistScalVec &shares, const std::vector<int> &indices, int t);

// recovers the shamir shares from t shares and t indices
// indices must have length at least t
// shares[0] must have length at least t
// shares[j][i] is client indices[i]'s share of the j-th shared element
RistScalVec shamir_recover(const RistScalMat &shares, const std::vector<int> &indices, int t);

#endif //RISEFL_CRYPTO_SHAMIR_H
