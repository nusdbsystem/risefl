//
// Created by yizheng on 8/3/23.
//

//#include <NTL/ZZ.h>
//#include <NTL/ZZ_p.h>
//#include <NTL/ZZ_pX.h>
//#include <NTL/vec_ZZ_p.h>
#include <vector>
#include <algorithm>
#include <execution>
//#include "zz_p_conversion.h"
#include "include/zkp.h"
#include "include/shamir.h"
#include "include/utils.h"

//using NTL::vec_ZZ_p;

/*
 * produces shamir'scalar_seed secret share
 * input:
 *  scalar_seed, the secret,
 *  num_clients, the number of people to share,
 *  t, the minimum number of people needed to recover the share
 * output:
 *  the ordered pair of
 *      num_clients+1 shares s0, ..., sn (s0=scalar_seed, s1r~sn are shares)
 *      the check strings Psi_0, ..., Psi_{t-1}
 *  note: si = f(i) for i = 0, ..., num_clients
 *        f has degree <= t-1
 */

RistScal evaluate(const RistScalVec &coeffs, const RistScal &p) {
    int deg = coeffs.size() - 1;
    RistScal out = coeffs[deg];
    for (int i = deg - 1; i >= 0; i--) {
        out = p * out + coeffs[i];
    }
    return out;
}

std::pair<RistScalVec, RistElemP3Vec> compute_shamir_shares_with_check_string(const RistScal &s, int n, int t) {
    assert(t > 0);

    // generate a random polynomial of degree at most t-1 with
    //  constant coeff = s
    RistScalVec coeffs(t);
    rand_init(coeffs);
    coeffs[0] = s;
    RistScalVec shares(n + 1);
    for (int i = 1; i <= n; i++) {
        shares[i] = evaluate(coeffs, RistScal(i));   // shares[0] is never used, never initialized
    }

    RistElemP3Vec checks(t);

    std::for_each(
            std::execution::par_unseq,
            checks.begin(),
            checks.end(),
            [&checks, &coeffs](auto &&check) {
                auto i = &check - checks.data();
                pedersen_zero_commit(check, coeffs[i]);
            });
    return {shares, checks};
}

BatchShamirShareCheckString compute_batch_shamir_shares_with_check_string(const RistScalVec &ss, int n, int t) {
    BatchShamirShareCheckString ret(ss.size(), n, t);
    std::for_each(std::execution::par_unseq, ss.begin(), ss.end(),
                  [&ss, &n, &t, &ret](auto &&s) {
                      auto k = &s - ss.data();
                      auto shares_with_check_string = compute_shamir_shares_with_check_string(s, n, t);
                      ret.shares[k] = shares_with_check_string.first;
                      ret.check_strings.elems[k] = shares_with_check_string.second;
                  });
    ret.check_strings.fill_bytes();
    //    std::for_each(std::execution::par_unseq, ss.begin(), ss.end(),
//                  [&ss, &n, &t, &ret](auto && s) {
//                      auto k = &s - ss.data();
//                      ret.check_string_coeffs[k] = rist_scalar_from_hash_from_bytes(ret.check_strings.bytes, bytes_from_int(k));
//                  });
//
//    for (int i = 0; i < ret.check_strings.size_rows(); i++) {
//        linear_comb(ret.batch_check_string.elems[i], ret.check_string_coeffs, ret.check_strings.slice_col_0().elems);
//    }
//    ret.batch_check_string.fill_bytes();

    return ret;
}

bool b_shamir_verify(const RistScal &share, const RistElemP3Vec &checks, int t, int i) {
    if (checks.size() != t)
        return false;
    RistScalVec powers = scalar_geometric_series(t, RistScal(i));
    RistElemP3 check_res = pedersen_zero_commit_p3(share);
    return check_res == linear_comb(powers, checks);
}

bool b_shamir_verify(const RistScalVec &shares, const RistElemP3Mat &checks, int t, int i) {
    int batch_size = shares.size();
    assert(batch_size == checks.size());
    assert(checks[0].size() == t);

    RistScalVec coeffs(batch_size);
    rand_init(coeffs);


    RistElemP3 left, right;
    pedersen_zero_commit(left, inner_prod(shares, coeffs));

    RistScalVec coeffs_dup(batch_size * t);
    RistScalVec powers = scalar_geometric_series(t, RistScal(i));

    for (int k = 0; k < batch_size; k++)
        for (int j = 0; j < t; j++)
            coeffs_dup[k * t + j] = coeffs[k] * powers[j];

    linear_comb(right, coeffs_dup, concatenate_double_vec(checks));

    return left == right;
}

RistScal shamir_recover(const RistScalVec &shares, const std::vector<int> &indices, int t) {
    RistScal indices_prod(1);
    for (int i = 0; i < indices.size(); i++) {
        indices_prod = indices_prod * RistScal(indices[i]);
    }
    RistScalVec from_single_share(t, indices_prod);
    RistScal to_divide;
    for (int i = 0; i < t; i++) {
        to_divide = RistScal(indices[i]);
        for (int j = 0; j < t; j++) {
            if (j != i) {
                to_divide = to_divide * RistScal(indices[j] - indices[i]);
            }
        }
        from_single_share[i] = from_single_share[i] / to_divide;
    }
    return inner_prod(shares, from_single_share);
}

RistScalVec shamir_recover(const RistScalMat &shares, const std::vector<int> &indices, int t) {
    RistScalVec ret;
    ret.reserve(shares.size());
    for (auto &&ss: shares) {
        ret.emplace_back(shamir_recover(ss, indices, t));
    }
    return ret;
}
