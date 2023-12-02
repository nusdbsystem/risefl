//
// Created by yizheng on 5/4/23.
//
#include "../include/rist_fast_computation.h"
#include "../include/mpz_conversion.h"
#include <chrono>
#include <cassert>

#include <cmath>

void test_mpz() {
    int n1 = -100000;
    MyMpz m1(n1);
    RistScal r1;
    rist_from_mpz(r1, m1);
    assert(r1 == RistScal(n1));

    int n2 = 3456;
    MyMpz m2(n2);
    RistScal r2;
    rist_from_mpz(r2, m1 + m2);
    assert(r2 == RistScal(n1 + n2));
    rist_from_mpz(r2, m1 - m2);
    assert(r2 == RistScal(n1 - n2));
    rist_from_mpz(r2, m1 * m2);
    assert(r2 == RistScal(n1 * n2));
//    rist_from_mpz(r2, m1.remainder(m2));
//    assert(r2 == RistScal(std::remainder(n1, n2)));

    MyMpz m3 = c_mpz_rist_order;
    m3.reduce_to_rist();
    assert(mpz_sgn(m3.pointer) == 0);

    int length = 10;
    RistScalVec rr(length), ss(length);
    rand_init(rr);
    rand_init(ss);
    MpzVec rr_mpz(length);
    MpzVec ss_mpz(length);
    mpz_vec_from_rist(rr_mpz, rr);
    mpz_vec_from_rist(ss_mpz, ss);
    MyMpz inner_prod_mpz;
    for (int i = 0; i < length; i++) {
        inner_prod_mpz += rr_mpz[i] * ss_mpz[i];
    }
    inner_prod_mpz.reduce_to_rist();
    RistScal inner_prod_rist = inner_prod(rr, ss);
    RistScal inner_prod_conv;
    rist_from_mpz(inner_prod_conv, inner_prod_mpz);
    assert(inner_prod_rist == inner_prod_conv);



    std::cout << "mpz test complete!\n";
}