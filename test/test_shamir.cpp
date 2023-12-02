//
// Created by yizheng on 8/3/23.
//
#include <iostream>
#include <chrono>
#include <cassert>
#include <vector>
#include <string>
#include <algorithm>
#include <execution>
#include <thread>
#include <sodium.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include "../include/ristretto.h"
#include "../include/ristretto_vector.h"
#include "../include/shamir.h"

using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pX;

void test_shamir() {
    auto secret = RistScal(1);
    int n = 100, t = 30;

    auto share_res = compute_shamir_shares_with_check_string(secret, n, t);
    auto shares = share_res.first;
    auto checks = share_res.second;

    for (int i = 1; i <= n; i++) {
        assert(b_shamir_verify(shares[i], checks, t, i));
    }

    assert(!b_shamir_verify(shares[1], checks, t, 2));
    assert(!b_shamir_verify(shares[2], checks, t, 1));
    checks[1] = checks[1] + checks[2];
    assert(!b_shamir_verify(shares[1], checks, t, 1));

    RistScalVec shares_for_recovery;
    shares_for_recovery.resize(t);
    std::vector<int> indices(t);
    for (int i = 0; i < t; i++) {
        indices[i] = 50 + i;
        shares_for_recovery[i] = shares[50 + i];
    }
    auto recovered_secret = shamir_recover(shares_for_recovery, indices, t);
    assert(secret == recovered_secret);

    shares_for_recovery[0] = shares_for_recovery[0] + c_scal_one;
    assert(secret != shamir_recover(shares_for_recovery, indices, t));

    // check linearity of Shamir sharing
    RistScal secret2(88);
//    ZZ_p secret2 = random_ZZ_p();
    auto share_res2 = compute_shamir_shares_with_check_string(secret2, n, t);
    auto shares2 = share_res2.first;
    auto checks2 = share_res2.second;

    for (int i = 0; i < t; i++) {
        indices[i] = 60 + i;
        shares_for_recovery[i] = shares[60 + i] + shares2[60 + i];
    }

    auto recovered_secret_sum = shamir_recover(shares_for_recovery, indices, t);

    assert(secret + secret2 == recovered_secret_sum);


    RistScalVec rrr(100);
//    auto sss = vec_zz_p_from_rist(rrr);
//    assert (rrr == rist_from_vec_zz_p(sss));
//    assert (sss == vec_zz_p_from_rist(rist_from_vec_zz_p(sss)));

    std::cout << "shamir test passed!\n";
}

