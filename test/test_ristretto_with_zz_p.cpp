//
// Created by yizheng on 10/3/23.
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
#include "../ristretto_with_zz_p.h"

void test_ristretto_with_zz_p() {
    RistElem q;
    rand_init(q);
    RistScal n;
    rand_init(n);
    ZZ_p m = ZZ_p_from_rist(n);
    assert(n * q == m * q);
    RistElem res;
    scalar_mult(res, m, q);
    assert(res == n * q);

    scalar_mult_base(res, m);
    assert(res == scalar_mult_base(n));
    assert(res == scalar_mult_base(m));

    std::cout << "Ristretto with ZZ_p test passed!" << std::endl;
}