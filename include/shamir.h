//
// Created by yizheng on 8/3/23.
//

#ifndef RISEFL_CRYPTO_SHAMIR_H
#define RISEFL_CRYPTO_SHAMIR_H

//#include <NTL/ZZ.h>
//#include <NTL/ZZ_p.h>
//#include <NTL/ZZ_pX.h>
#include <chrono>
#include <cassert>
#include <vector>

#include "ristretto.h"
#include "ristretto_vector.h"
//#include "zz_p_conversion.h"
#include "zkp.h"
#include "rist_fast_computation.h"
//
//using NTL::ZZ;
//using NTL::ZZ_p;
//using NTL::ZZ_pX;
//using NTL::vec_ZZ_p;

std::pair<RistScalVec, RistElemP3Vec> compute_shamir_shares_with_check_string(const RistScal &s, int n, int t);

class BatchShamirShareCheckString {
public:
    // (batch_size, n + 1)
    std::vector<RistScalVec> shares;

    // (batch_size, t)
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

BatchShamirShareCheckString compute_batch_shamir_shares_with_check_string(const RistScalVec &ss, int n, int t);

//inline std::pair<RistScalVec, RistElemP3Vec> compute_shamir_shares_with_check_string(const RistScal &s, int n, int t) {
//    return compute_shamir_shares_with_check_string(ZZ_p_from_rist(s), n, t);
//}

bool b_shamir_verify(const RistScal &share, const RistElemP3Vec &checks, int t, int i);

bool b_shamir_verify(const RistScalVec &shares, const RistElemP3Mat &checks, int t, int i);

//RistScalVec rist_from_vec_zz_p(const vec_ZZ_p &vv);
//
//vec_ZZ_p vec_zz_p_from_rist(const RistScalVec &ss);

//ZZ_p shamir_recover(const vec_ZZ_p &shares, const std::vector<int> &indices, int t);

RistScal shamir_recover(const RistScalVec &shares, const std::vector<int> &indices, int t);

RistScalVec shamir_recover(const RistScalMat &shares, const std::vector<int> &indices, int t);

#endif //RISEFL_CRYPTO_SHAMIR_H
