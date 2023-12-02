//
// Created by yizheng on 27/3/23.
//
#include "../include/rist_fast_computation.h"
#include <iostream>
#include <chrono>
#include <cassert>
#include "../bench/simulate.h"


void test_precomp() {
    RistElem r;
    rand_init(r);

    RistElemP3 r_p3;
    p3_from_bytes(r_p3, r.element);

    RistElemPrecompTable table;
    generate_precomp_table(table, r_p3);

    RistScal n;
    rand_init(n);

    RistElemP3 s1, s2;
    std::cout << "normal scalar mult: " << measure_time([&](){    scalar_mult(s1, n, r_p3); })
     << "\nprecomp scalar mult:" << measure_time([&]() { scalar_mult_precomp(s2, n, table); }) << "\n";

    assert(s1 == s2);

    std::cout << "precomp test complete!\n";
}