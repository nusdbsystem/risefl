//
// Created by yizheng on 4/4/23.
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
#include "../include/di_zkp_interface_client.h"
#include "../include/di_zkp_interface_common.h"
#include "../include/di_zkp_interface_server.h"
#include "../include/ristretto.h"
#include "../include/ristretto_vector.h"
#include "simulate.h"
//#include "NTL/ZZ.h"
//#include "NTL/ZZ_p.h"
//#include <NTL/vec_ZZ.h>
//#include <NTL/vec_ZZ_p.h>
#include <random>

#include <gmpxx.h>
#include <gmp.h>

#include "../include/mpz_conversion.h"

//using NTL::ZZ;
//using NTL::ZZ_p;
//using NTL::vec_ZZ_p;



void bench_inner_prod_mod_p(const RistScalVec& aa_bb, const RistScalVec& bb, const std::vector<std::vector<int>>& aa_pos,
                            MpzVec& aa_bb_mpz, const MpzVec& bb_mpz, const MpzMat& aa_pos_mpz) {
    auto k = aa_pos.size();
    auto d = aa_pos[0].size();

//    MpzVec aa_bb_mpz(d);
//    MpzVec bb_mpz(k);
//    MpzMat aa_pos_mpz(k, MpzVec(d));
//    mpz_vec_from_rist(aa_bb_mpz, aa_bb);
//    mpz_vec_from_rist(bb_mpz, bb);
//    for (int i = 0; i < k; i++) {
//        mpz_vec_from_rist(aa_pos_mpz[i], aa_pos[i]);
//    }
//
    std::for_each(
            std::execution::par_unseq,
            aa_bb.begin(),
            aa_bb.end(),
            [&](auto&& v) {
                auto i = &v - aa_bb.data();
                for (int m = 0; m < k; m++) {
                    aa_bb_mpz[i] += MyMpz(aa_pos[m][i]) * bb_mpz[m];
                }
            }
    );
//
//
//
//
//    ZZ_p c_zz;
//    InnerProduct(c_zz, aa_zz, bb_zz);
//
//    RistScal c;
//    rist_from_zz_p(c, c_zz);
//
//    RistScal u;
//    rand_init(u);
//    auto v = u * aa;
//    auto v2 = u * bb;
}

void bench_inner_prod_mod_p() {
    int d = 1e5;
    int k = 25;
    RistScalVec aa_bb(d);
    RistScalVec bb(k);
    std::vector<std::vector<int>> aa_pos(k, std::vector<int>(d));
    std::cout << "1" << std::endl;

    rand_init(bb);

    for (auto &&v : aa_pos){
        std::mt19937_64 generator;
        generator.seed(111);
        std::uniform_int_distribution<int> unif_int;


//    for (auto &&a: aa_0) {
//        std::uniform_int_distribution<unsigned long> unif_unsigned_long;
//        RistHashbytes hash;
//        for (int j = 0; j < hash.hashbytes.size(); j += sizeof(unsigned long)) {
//            auto u = unif_unsigned_long(generator);


        for (auto && u : v) {
            u = unif_int(generator);
            u >>= 16;
        }
    }

    std::cout << "2" << std::endl;

    MpzVec aa_bb_mpz(d);
    std::cout << "2.1" << std::endl;

    MpzVec bb_mpz(k);
    std::cout << "2.2" << std::endl;
    MpzMat aa_pos_mpz(k, MpzVec(d));
//    mpz_vec_from_rist(aa_bb_mpz, aa_bb);

    std::cout << "3" << std::endl;

    mpz_vec_from_rist(bb_mpz, bb);
//    mpz_mat_from_int(aa_pos_mpz, aa_pos);

    std::cout << "4" << std::endl;

    std::cout << "\n" << measure_time([&](){ bench_inner_prod_mod_p(aa_bb, bb, aa_pos, aa_bb_mpz, bb_mpz, aa_pos_mpz); }) << "\n";
}