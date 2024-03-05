#include <iostream>
#include <chrono>
#include <cassert>
#include <vector>
#include <string>
#include <algorithm>
#include <execution>
#include <thread>
#include <sodium.h>
#include "../include/rist_fast_computation.h"
#include "../include/risefl_interface_client.h"
#include "../include/risefl_interface_common.h"
#include "../include/risefl_interface_server.h"
#include "simulate.h"


void bench_precomp() {
    int length = 1000;
    RistElemP3Vec rr_p3(length);
    RistElemVec rr(length);
    for (int i = 0; i < rr_p3.size(); i++) {
        p3_from_bytes(rr_p3[i], rr[i].element);
    }

    int per_elem_scal = 1000;
    RistScalMat vv1(length, RistScalVec(per_elem_scal));
    RistScalMat vv2(length, RistScalVec(per_elem_scal));
    for (auto && v : vv1)
        rand_init(v);
    for (auto && v : vv2)
        rand_init(v);
    std::vector<std::vector<int>> vv1_int(length, std::vector<int>(per_elem_scal));

    std::mt19937_64 generator;
    generator.seed(333);
    std::uniform_int_distribution<unsigned long> distrib(-30000, 30000);
    for (auto && row : vv1_int)
        for (auto && ele : row)
            ele = distrib(generator);

    RistElemP3Mat result(length, RistElemP3Vec(per_elem_scal));
    RistElemP3Mat result_precomp(length, RistElemP3Vec(per_elem_scal));

    std::cout << "mult without precomp: " << measure_time([&](){
        for (int i = 0; i < length; i++) {
            std::for_each(std::execution::par_unseq, result_precomp[i].begin(), result_precomp[i].end(),
                          [&](auto && res){
                              auto j = &res - result_precomp[i].data();
//                              result[i][j] = vv2[i][j] * rr_p3[i];
                              pedersen_commit(result[i][j], RistScal(vv1_int[i][j]), rr_p3[i], vv2[i][j]);
                          });
        }
    });

    std::vector<RistElemPrecompTable> precomp(length);
    for (int i = 0; i < length; i++)
        generate_precomp_table(precomp[i], rr_p3[i]);


    std::cout << "\nmult with precomp: " << measure_time([&](){
        for (int i = 0; i < length; i++) {
            std::for_each(std::execution::par_unseq, result_precomp[i].begin(), result_precomp[i].end(),
                          [&](auto && res){
                auto j = &res - result_precomp[i].data();
//                scalar_mult_precomp(result_precomp[i][j], vv2[i][j], precomp[i]);
                pedersen_commit_precomp(result_precomp[i][j], vv1_int[i][j], 16, precomp[i], vv2[i][j]);
            });
        }
    });

    assert(result == result_precomp);
}