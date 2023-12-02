//
// Created by yizheng on 14/6/23.
//

#include "../include/bulletin.h"

void test_bulletin() {
    auto key_pair = gen_sign_key_pair();
    auto pub = key_pair.first;
    auto prv = key_pair.second;
    std::vector<unsigned char> mes{10, 30, 40, 9};
    auto signed_mes = sign(mes, prv);
    assert(sign_open(signed_mes, pub) == mes);

    signed_mes.emplace_back(2);
    bool exception_thrown = false;
    try {
        sign_open(signed_mes, pub);
    }
    catch (InvalidSign&) {
        exception_thrown = true;
    }
    assert(exception_thrown);

    std::cout << "bulletin test complete!\n";
}
