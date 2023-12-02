//
// Created by yizheng on 6/3/23.
//
#include <sodium.h>
#include <iostream>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <chrono>
#include <cassert>
#include <vector>
#include <string>
#include <algorithm>
#include <execution>
#include <thread>
#include "../include/ristretto.h"
#include "../include/ristretto_vector.h"


using namespace std::chrono_literals;
typedef std::chrono::high_resolution_clock Clock;

using namespace NTL;

void test_ristretto_element() {
    RistElem p, p_prime;
    rand_init(p);
    rand_init(p_prime);
    assert(is_valid(p) == 1);
    assert(p != p_prime);

    RistElem q;
    RistHashbytes r;
    hash_init(q, r);
    assert(is_valid(q) == 1);

    RistScal n, n_prime;
    rand_init(n);
    rand_init(n_prime);
    assert (n != n_prime);
    RistElem p2 = n * q;
    assert(is_valid(p2) == 1);

    RistScal m;
    rand_init(m);
    RistElem p3 = scalar_mult_base(m);
    assert(is_valid(p3) == 1);

    assert(p + q == q + p);
    assert((p - p2) - p3 == p - (p2 + p3));
    assert((p - p2) + p3 == p - (p2 - p3));
    assert(p + p3 - p3 == p);
    assert(p - p3 + p3 == p);

    assert((m + n) * p2 == m * p2 + n * p2);
    assert((m - n) * p2 == m * p2 - n * p2);

    RistScal k;
    rand_init(k);
    assert(n + m == m + n);
    assert((n + m) + k == n + (m + k));
    assert((n - m) + k == n - (m - k));
    assert(n - k + k == n);
    assert(n + m - m == n);
    assert(n - m == -(m - n));
    assert(-(m + n) == -m - n);

    RistScal scalar_zero = (n + (-n));
    RistScal scalar_one = (n + one_minus(n));
    assert(n / n == scalar_one);
    assert(n * scalar_zero == scalar_zero);
    assert(scalar_zero + scalar_one == scalar_one);
    assert(scalar_zero / scalar_one == scalar_zero);
    assert(scalar_one * scalar_zero == scalar_zero);

    assert(n * m == m * n);
    assert(n * (m * k) == (n * m) * k);
    assert(n * (m + k) == n * m + n * k);
    assert(n * (m - k) == n * m - n * k);
    assert(n / m == scalar_invert(m / n));
    assert(m + n / k == (m * k + n) / k);
    RistScal j;
    rand_init(j);
    assert(m / j + n / k == (m * k + n * j) / (j * k));

    assert(scalar_one * p == p);
    assert(p + p == (scalar_one + scalar_one) * p);
    RistElem zero = scalar_zero * p;
    assert(p3 - p3 == zero);
    assert(n * zero == zero);
    assert(zero == c_elem_zero);

    RistScal n4, n5;
    n4 = scalar_reduce(r);
    scalar_reduce(n5, r);
    assert(n4 == n5);

    assert(scalar_one == RistScal(1));
    assert(scalar_one + RistScal(-1) == RistScal(0));
    assert(scalar_zero == RistScal(0));
    assert(scalar_one + scalar_one == RistScal(2));

    auto power = RistScal(2);

    for (int i = 1; i < 10; i++) {
        power = power + power;
    }
    assert(power == RistScal(1 << 10));
    for (int i = 10; i < 20; i++) {
        power = power + power;
    }
    assert(power == RistScal(1 << 20));

    std::cout << "Ristretto element test passed!" << std::endl;
}

void test_ristretto_vector() {
    int length = 10;
    RistElemVec rr(length), pp(length), qq(length);
    assert(rr.size() == length);
    rand_init(pp);
    rand_init(qq);
    assert (pp == rist_elem_vec_from_bytes(bytes_from_rist_elem_vec(pp)));

    for (auto p: pp) {
        assert(is_valid(p) == 1);
    }
    add(rr, pp, qq);
    for (auto p: rr) {
        assert(is_valid(p) == 1);
    }
    for (int i = 0; i < length; i++) {
        assert(rr[i] == pp[i] + qq[i]);
    }

    RistHashbytesVec hh(length);
    RistElemVec pp1(length), pp2(length);
    hash_init(pp1, hh);
    hash_init(pp2, hh);
    for (auto p: pp1) {
        assert(is_valid(p) == 1);
    }
    RistHashbytesVec hh3(length);
    for (int i = 0; i < length; i++) {
        hh3[i].hashbytes[0] = hh[i].hashbytes[0] + 1;
    }
    RistElemVec pp3(length);
    hash_init(pp3, hh3);
    for (int i = 0; i < length; i++) {
        assert(!(pp3[i] == pp2[i]));
    }

    RistScalVec nn(length), nn2(length), mm(length), kk(length);
    rand_init(mm);
    assert (mm == rist_scal_vec_from_bytes(bytes_from_rist_scal_vec(mm)));
    rand_init(kk);
    scalar_add(nn, mm, kk);
    for (int i = 0; i < length; i++) {
        assert(nn[i] == mm[i] + kk[i]);
    }
    nn2 = mm + kk;
    assert(nn == nn2);
    scalar_sub(nn, mm, kk);
    for (int i = 0; i < length; i++) {
        assert(nn[i] == mm[i] - kk[i]);
    }
    nn2 = mm - kk;
    assert(nn == nn2);
    scalar_mul(nn, mm, kk);
    for (int i = 0; i < length; i++) {
        assert(nn[i] == mm[i] * kk[i]);
    }
    nn2 = mm * kk;
    assert(nn == nn2);

    scalar_mult(pp3, nn, pp2);
    for (int i = 0; i < length; i++) {
        assert(pp3[i] == nn[i] * pp2[i]);
    }
    scalar_mult_base(pp3, kk);
    for (int i = 0; i < length; i++) {
        assert(pp3[i] == scalar_mult_base(kk[i]));
    }

    RistNonreducedScalarVec ss(length);
    scalar_reduce(nn, ss);
    scalar_reduce(mm, *reinterpret_cast<const RistHashbytesVec *>(&ss));
    assert(nn == mm);

    RistElemVec to_sum(length);
    rand_init(to_sum);
    auto sum_elem = sum(to_sum);
    RistElem sum_seq = to_sum[0];

    for (int i = 1; i < length; i++) {
        add(sum_seq, sum_seq, to_sum[i]);
    }
    assert(sum_seq == sum_elem);


    //benchmark

    RistScalVec to_sum_scalar(length);
    rand_init(to_sum_scalar);
    auto sum_elem_scalar = sum(to_sum_scalar);
    RistScal sum_seq_scalar = to_sum_scalar[0];

    for (int i = 1; i < length; i++) {
        scalar_add(sum_seq_scalar, sum_seq_scalar, to_sum_scalar[i]);
    }
    assert(sum_seq_scalar == sum_elem_scalar);

    auto t1 = Clock::now();
    auto comb = linear_comb(to_sum_scalar, to_sum);
    auto t2 = Clock::now();

    long double mytime;
    mytime = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / double(1000000000);

    std::cout << "linear comb cost:" << mytime << std::endl;

    auto comb_seq = to_sum_scalar[0] * to_sum[0];
    for (int i = 1; i < length; i++) {
        add(comb_seq, comb_seq, to_sum_scalar[i] * to_sum[i]);
    }
    assert(comb == comb_seq);


    RistScalVec to_sum_scalar_2(length);
    rand_init(to_sum_scalar_2);
    auto prod = inner_prod(to_sum_scalar, to_sum_scalar_2);
    auto prod_seq = to_sum_scalar[0] * to_sum_scalar_2[0];
    for (int i = 1; i < length; i++) {
        scalar_add(prod_seq, prod_seq, to_sum_scalar[i] * to_sum_scalar_2[i]);
    }
    assert(prod == prod_seq);


    t1 = Clock::now();
    for (int dddd = 0; dddd < 1000; dddd++)
        auto ggggg = to_sum_scalar[0] * to_sum[0];
    t2 = Clock::now();

    mytime = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / double(1000000000);

    std::cout << "1000 scalar mul takes:" << mytime << std::endl;


    t1 = Clock::now();
    for (int dddd = 0; dddd < 1000; dddd++)
        to_sum_scalar[0] = to_sum_scalar[0] * to_sum_scalar[1];
    t2 = Clock::now();

    mytime = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / double(1000000000);

    std::cout << "1000 scalar product takes:" << mytime << std::endl;


    t1 = Clock::now();
    for (int dddd = 0; dddd < 1000; dddd++)
        auto pptppfw = to_sum[0] + to_sum[1];
    t2 = Clock::now();

    mytime = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / double(1000000000);

    std::cout << "1000 adds takes:" << mytime << std::endl;


    std::cout << "Ristretto vector test passed!" << std::endl;
}

