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
#include "../include/zkp_hash.h"


void test_zkp_hash() {

    std::vector<unsigned char> a{12, 34, 5}, b{66}, c{};
    auto h1 = hash_from_bytes(a);
    auto h2 = hash_from_bytes(a, c);
    assert(h1 == h2);

    assert(hash_from_bytes() == hash_from_bytes(c));
    assert(hash_from_bytes(a, a, b, c, c) == hash_from_bytes(a, a, b));

    std::vector<unsigned char> a1{12, 34}, b1{5, 66};
    assert(hash_from_bytes(a1, c, b1) == hash_from_bytes(c, a, b));

    std::vector<unsigned char> key{88, 0};
    assert(hash_from_bytes_with_key(key, a1, c, b1) == hash_from_bytes_with_key(key, c, a, b, c));

    assert(rist_scalar_from_hash_from_bytes(c, a, b) == rist_scalar_from_hash_from_bytes(a, c, c, b));
    assert(rist_element_from_hash_from_bytes(c, a, b) == rist_element_from_hash_from_bytes(a, c, c, b));

    std::cout << "ZKP hash test passed!" << std::endl;

}