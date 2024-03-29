//
// Created by yizheng on 10/3/23.
//

// zkp_hash.h
//
// the hash functions used in ZKP proof generation and verification

#ifndef RISEFL_CRYPTO_ZKP_HASH_H
#define RISEFL_CRYPTO_ZKP_HASH_H

#include <sodium.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include "ristretto.h"
#include "ristretto_vector.h"

/*
 * hash function wrapper
 * https://doc.libsodium.org/hashing/generic_hashing
*/

// concatenate an unknown number of vectors into one
template<typename... Ts>
std::vector<unsigned char> concatenate(const Ts &... ts) {
    std::vector<unsigned char> bytes;
    (bytes.insert(bytes.end(), ts.begin(), ts.end()), ...);
    return bytes;
}

// returns hash from a key and an unknown number of vectors
template<typename... Ts>
RistHashbytes hash_from_bytes_with_key(const std::vector<unsigned char> &key, const Ts &... ts) {
    std::vector<unsigned char> bytes = concatenate(std::forward<const Ts>(ts)...);
//    (bytes.insert(bytes.end(), ts.begin(), ts.end()), ...);
    RistHashbytes hash;
    crypto_generichash(hash.hashbytes.data(), crypto_core_ristretto255_HASHBYTES,
                       bytes.data(), bytes.size(), key.data(), key.size());
    return hash;
}

// returns hash from empty key and an unknown number of vectors
template<typename... Ts>
RistHashbytes hash_from_bytes(const Ts &... ts) {
    return hash_from_bytes_with_key(std::vector<unsigned char>{}, std::forward<const Ts>(ts)...);
}

// hash an unknown number of vectors into a Ristretto scalar
template<typename... Ts>
RistScal rist_scalar_from_hash_from_bytes(const Ts &... ts) {
    return scalar_reduce(hash_from_bytes(std::forward<const Ts>(ts)...));
}

// hash an unknown number of vectors into a Ristretto group member
template<typename... Ts>
RistElem rist_element_from_hash_from_bytes(const Ts &... ts) {
    RistElem h;
    hash_init(h, hash_from_bytes(std::forward<const Ts>(ts)...));
    return h;
}

#endif //RISEFL_CRYPTO_ZKP_HASH_H
