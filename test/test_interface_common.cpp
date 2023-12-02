//
// Created by yizheng on 16/3/23.
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
#include "../include/di_zkp_interface_common.h"


using namespace std::chrono_literals;
typedef std::chrono::high_resolution_clock Clock;

void test_interface_common() {
    auto key_pair = generate_dh_key_pair();
    auto pub = key_pair.first;
    auto prv = key_pair.second;

    std::vector<unsigned char> mes = {1, 3, 9, 222};
    auto cip_with_nonce = encrypt(mes, pub, prv);
    assert (decrypt(cip_with_nonce, pub, prv) == mes);

    auto key_pair2 = generate_dh_key_pair();
    auto pub2 = key_pair2.first;
    auto prv2 = key_pair2.second;
    assert(pub.pub != pub2.pub);
    assert(prv.prv != prv2.prv);
    cip_with_nonce = encrypt(mes, pub, prv2);
    assert (decrypt(cip_with_nonce, pub2, prv) == mes);
    assert (decrypt(cip_with_nonce, pub2, prv2) != mes);
    auto cip_with_nonce2 = encrypt(std::vector<unsigned char>(mes.begin(), mes.end()), pub, prv2);
    assert (decrypt(cip_with_nonce, pub2, prv) == mes);

    std::cout << "Diffie-Hellman encryption test passed!\n";

    std::cout << "interface commom test passed!\n";
}