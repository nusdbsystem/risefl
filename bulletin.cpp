//
// Created by yizheng on 14/6/23.
//

#include "include/bulletin.h"

std::pair<SignPubKey, SignPrvKey> gen_sign_key_pair() {
    SignPubKey pub;
    SignPrvKey prv;
    crypto_sign_keypair(pub.key.data(), prv.key.data());
    return {pub, prv};
}

std::vector<unsigned char> sign(const std::vector<unsigned char> &mes, const SignPrvKey &prv) {
    std::vector<unsigned char> ret(SIGNBYTES + mes.size());
    crypto_sign(ret.data(), NULL, mes.data(), mes.size(), prv.key.data());
    return ret;
}

std::vector<unsigned char> sign_open(const std::vector<unsigned char> &signed_mes, const SignPubKey &pub) {
    if (signed_mes.size() - SIGNBYTES < 0) {
        throw InvalidSign();
    }
    std::vector<unsigned char> ret(signed_mes.size() - SIGNBYTES);
    if (crypto_sign_open(ret.data(), NULL, signed_mes.data(), signed_mes.size(), pub.key.data()) < 0) {
        throw InvalidSign();
    }
    return ret;
}
