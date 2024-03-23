//
// Created by yizheng on 14/6/23.
//

// bulletin.h 
//
// On a bulletin board, every client generates a pair consisting of a private key and a public key. Every message sent by the client can be signed with the private key, which generates a signed message. The signed message can be verified with the public key for its authenticity. Authenticity means that the message is indeed sent from the client that holds the private key. If the authenticity is verified, the signed message can be converted back into the original message.

// We use libsodium to implement this mechanism. It is essentially a C++ wrapper of:
// https://doc.libsodium.org/public-key_cryptography/public-key_signatures


#ifndef RISEFL_CRYPTO_BULLETIN_H
#define RISEFL_CRYPTO_BULLETIN_H

#include <array>
#include <vector>
#include "sodium.h"
#include "exception.h"

constexpr int SIGNPUBKEYBYTES = crypto_sign_PUBLICKEYBYTES;
constexpr int SIGNPRVKEYBYTES = crypto_sign_SECRETKEYBYTES;
constexpr int SIGNBYTES = crypto_sign_BYTES;

// public key which is used to verify the authenticity of messages
struct SignPubKey {
    std::array<unsigned char, SIGNPUBKEYBYTES> key;
};

// private key which is used to sign messages
struct SignPrvKey {
    std::array<unsigned char, SIGNPRVKEYBYTES> key;
};

// generates a pair consisting of a public key and a private key
std::pair<SignPubKey, SignPrvKey> gen_sign_key_pair();

// sign a message using a private key
std::vector<unsigned char> sign(const std::vector<unsigned char> &mes,
                                const SignPrvKey &prv);

// verify the authenticity of the signed message 
//      if yes: return the original message prior to signing
//      if no: throws InvalidSign() exception 
std::vector<unsigned char> sign_open(const std::vector<unsigned char> &signed_mes,
                                     const SignPubKey &pub);

#endif //RISEFL_CRYPTO_BULLETIN_H
