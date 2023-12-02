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
#include "../include/ristretto.h"
#include "../include/ristretto_vector.h"
#include "../zz_p_conversion.h"
#include "../include/zkp.h"

void test_zkp() {
    int length = 1000;
    RistScalVec xx(length);
    RistScal r;
    RistElemVec hh(length);
    rand_init(xx);
    rand_init(r);
    rand_init(hh);

    typedef std::chrono::high_resolution_clock Clock;

    auto t1 = Clock::now();
    auto yy = pedersen_commit(xx, hh, r);
    auto t2 = Clock::now();

    long double mytime;
    mytime = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / double(1000000000);
    std::cout << "Pedersen commitments took " << mytime * 1000 << " microseconds on average\n";

    for (int i = 0; i < length; i++) {
        assert(yy[i] == pedersen_commit(xx[i], hh[i], r));
    }

    RistElemVec yy2(length);
    t1 = Clock::now();
    pedersen_commit(yy2, xx, hh, r);
    t2 = Clock::now();
    mytime = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / double(1000000000);

    std::cout << "Pedersen commitments on existing Ristretto element vector took " << mytime * 1000
              << " microseconds on average\n";
    assert(yy == yy2);

    pedersen_commit(yy[0], xx[0], hh[0], r);
    assert(yy == yy2);

    assert(yy[0] == pedersen_zero_commit(xx[0]) + r * hh[0]);
    RistElem y;
    pedersen_zero_commit(y, xx[0]);
    assert(yy[0] == y + r * hh[0]);


    auto y0 = pedersen_zero_commit(r);
    auto y1 = pedersen_commit(xx[0], hh[0], r);
    auto proof = generate_pedersen_with_zero_proof(hh[0], y0, y1, r, xx[0]);
    assert(b_verify_pedersen_with_zero(hh[0], y0, y1, proof));

    y = pedersen_zero_commit(r);
    yy = pedersen_commit(xx, hh, r);
    auto vec_proof = generate_pedersen_vec_with_zero_proof(hh, y, yy, r, xx);
    assert(b_verify_pedersen_vec_with_zero(hh, y, yy, vec_proof));

    RistScal r1, r2;
    rand_init(r1);
    rand_init(r2);
    y1 = pedersen_commit(xx[1], hh[0], r1);
    auto y2 = pedersen_commit(xx[1] * xx[1], hh[0], r2);
    auto proof2 = generate_pedersen_with_square_proof(hh[0], y1, y2, xx[1], r1, r2);
    assert(b_verify_pedersen_with_square(hh[0], y1, y2, proof2));


    // range proof
    {
        int n = 8;
        RistElemVec gg(n), hh(n);
        rand_init(gg);
        rand_init(hh);
        RistElem h;
        rand_init(h);
        RistScal gamma;
        rand_init(gamma);
        auto v = RistScal(200);
        auto V = pedersen_commit(v, h, gamma);
        auto proof = generate_range_proof_power_two(h, V, gg, hh, n, gamma, v);
        assert(b_verify_range_power_two(h, V, gg, hh, n, proof));

        for (int i = 0; i < 256; i++) {
            v = RistScal(i);
            V = pedersen_commit(v, h, gamma);
            assert(b_verify_range_power_two(h, V, gg, hh, n,
                                            generate_range_proof_power_two(h, V, gg, hh, n, gamma, v)));
        }
        for (int i = 256; i < 290; i++) {
            v = RistScal(i);
            V = pedersen_commit(v, h, gamma);
            assert(!b_verify_range_power_two(h, V, gg, hh, n,
                                             generate_range_proof_power_two(h, V, gg, hh, n, gamma, v)));
        }

        n = 14;
        RistElemVec gg1(n), hh1(n);
        rand_init(gg1);
        rand_init(hh1);
        v = RistScal(1 << 12);
        V = pedersen_commit(v, h, gamma);
        assert(b_verify_range_power_two(h, V, gg1, hh1, n,
                                        generate_range_proof_power_two(h, V, gg1, hh1, n, gamma, v)));
        v = RistScal(1 << 14);
        V = pedersen_commit(v, h, gamma);
        assert(!b_verify_range_power_two(h, V, gg1, hh1, n,
                                         generate_range_proof_power_two(h, V, gg1, hh1, n, gamma, v)));

        n = 24;
        RistElemVec gg2(n), hh2(n);
        rand_init(gg2);
        rand_init(hh2);
        v = RistScal(1 << 23);
        V = pedersen_commit(v, h, gamma);
        assert(b_verify_range_power_two(h, V, gg2, hh2, n,
                                        generate_range_proof_power_two(h, V, gg2, hh2, n, gamma, v)));
        v = RistScal(1 << 24);
        V = pedersen_commit(v, h, gamma);
        assert(!b_verify_range_power_two(h, V, gg2, hh2, n,
                                         generate_range_proof_power_two(h, V, gg1, hh2, n, gamma, v)));


        n = 1;
        RistElemVec gg3(n), hh3(n);
        rand_init(gg3);
        rand_init(hh3);
        for (int i = 0; i < 2; i++) {
            v = RistScal(i);
            V = pedersen_commit(v, h, gamma);
            assert(b_verify_range_power_two(h, V, gg3, hh3, n,
                                            generate_range_proof_power_two(h, V, gg3, hh3, n, gamma, v)));
        }
        for (int i = 2; i < 10; i++) {
            v = RistScal(i);
            V = pedersen_commit(v, h, gamma);
            assert(!b_verify_range_power_two(h, V, gg3, hh3, n,
                                             generate_range_proof_power_two(h, V, gg3, hh3, n, gamma, v)));
        }

        n = 4;
        RistElemVec gg4(n), hh4(n);
        rand_init(gg4);
        rand_init(hh4);
        RistScal bound = RistScal(12);
        for (int i = 0; i < 12; i++) {
            v = RistScal(i);
            V = pedersen_commit(v, h, gamma);
            assert(b_verify_range(h, V, gg4, hh4, bound, generate_range_proof(h, V, gg4, hh4, bound, gamma, v)));
        }
        for (int i = 12; i < 22; i++) {
            v = RistScal(i);
            V = pedersen_commit(v, h, gamma);
            assert(!b_verify_range(h, V, gg4, hh4, bound, generate_range_proof(h, V, gg4, hh4, bound, gamma, v)));
        }
        for (int i = 1; i < 12; i++) {
            v = RistScal(-i);
            V = pedersen_commit(v, h, gamma);
            assert(!b_verify_range(h, V, gg4, hh4, bound, generate_range_proof(h, V, gg4, hh4, bound, gamma, v)));
        }

        std::cout << "range proof success!" << std::endl;
    }
    std::cout << "zkp test success!\n";
}