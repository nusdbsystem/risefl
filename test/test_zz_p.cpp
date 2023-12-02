//
// Created by yizheng on 28/2/23.
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

using NTL::ZZ;
using NTL::ZZ_p;

void test_ZZ_p_conversion() {
    RistScal zero_scalar, w;
    rand_init(w);
    zero_scalar = w - w;
    ZZ_p zero_number = ZZ_p_from_rist(zero_scalar);
    RistScal zero_scalar_2 = rist_from_ZZ_p(zero_number);
    assert(zero_scalar == zero_scalar_2);

    RistScal r, r2;
    rand_init(r);
    ZZ_p q = ZZ_p_from_rist(r);
    r2 = rist_from_ZZ_p(q);
    assert(r == r2);

    ZZ_p q3 = NTL::random_ZZ_p();
    RistScal r3 = rist_from_ZZ_p(q3);
    assert(q3 == ZZ_p_from_rist(r3));
    ZZ_p q4;
    zz_p_from_rist(q4, r3);
    assert(q3 == q4);

    assert(ZZ_p_from_rist(r2 + r3) == ZZ_p_from_rist(r2) + ZZ_p_from_rist(r3));
    assert(ZZ_p_from_rist(r2 - r3) == ZZ_p_from_rist(r2) - ZZ_p_from_rist(r3));
    assert(ZZ_p_from_rist(r2 * r3) == ZZ_p_from_rist(r2) * ZZ_p_from_rist(r3));
    assert(ZZ_p_from_rist(r2 / r3) == ZZ_p_from_rist(r2) / ZZ_p_from_rist(r3));

    std::cout << "ZZ_p test passed!" << std::endl;


}

