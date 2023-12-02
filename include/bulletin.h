//
// Created by yizheng on 14/6/23.
//

#ifndef DI_ZKP_CRYPTO_BULLETIN_H
#define DI_ZKP_CRYPTO_BULLETIN_H

#include <array>
#include <vector>
#include "sodium.h"
#include "exception.h"

constexpr int SIGNPUBKEYBYTES = crypto_sign_PUBLICKEYBYTES;
constexpr int SIGNPRVKEYBYTES = crypto_sign_SECRETKEYBYTES;
constexpr int SIGNBYTES = crypto_sign_BYTES;

struct SignPubKey {
    std::array<unsigned char, SIGNPUBKEYBYTES> key;
};

struct SignPrvKey {
    std::array<unsigned char, SIGNPRVKEYBYTES> key;
};

std::pair<SignPubKey, SignPrvKey> gen_sign_key_pair();

std::vector<unsigned char> sign(const std::vector<unsigned char> &mes,
                                const SignPrvKey &prv);

// throws InvalidSign() exception if sign is invalid
std::vector<unsigned char> sign_open(const std::vector<unsigned char> &signed_mes,
                                     const SignPubKey &pub);

#endif //DI_ZKP_CRYPTO_BULLETIN_H
