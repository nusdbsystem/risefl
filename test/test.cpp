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
#include "../include/ristretto.h"
#include "../include/ristretto_vector.h"
//#include "../zz_p_conversion.h"
//#include "test_ristretto.cpp"
//#include "test_zz_p.cpp"
//#include "test_shamir.cpp"
//#include "test_zkp.cpp"
//#include "test_zkp_hash.cpp"
//#include "test_ristretto_with_zz_p.cpp"
//#include "test_interface_common.cpp"
//#include "test_interface.cpp"
//#include "test_fast.cpp"
#include "test_interface_abstract.cpp"
//#include "test_mpz.cpp"
//#include "test_precomp.cpp"
#include "test_zkp_batch.cpp"
#include <limits>
#include <random>


void test() {
    // initialize the Ristretto group and the order p of ZZ_p
    assert(sodium_init() >= 0);
//    test_ristretto_element();
//    test_ristretto_vector();
//    test_ZZ_p_conversion();
//    test_ristretto_with_zz_p();
//    test_shamir();
//    test_zkp();
//    test_zkp_hash();
//
//    test_interface_common();
//    test_fast();
//    test_mpz();
//    test_precomp();
//    test_interface();
//    test_zkp_batch();
    test_interface_abstract();

    std::cout << "All tests passed!" << std::endl;

}