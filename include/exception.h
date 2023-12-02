//
// Created by yizheng on 29/3/23.
//

#ifndef DI_ZKP_CRYPTO_EXCEPTION_H
#define DI_ZKP_CRYPTO_EXCEPTION_H

#include <exception>

class DeserializationError : std::exception {
public:
    DeserializationError() = default;
};

class Abort : std::exception {
public:
    Abort() = default;
};

class DecryptionError : std::exception {
public:
    DecryptionError() = default;
};

class InvalidSign : std::exception {
public:
    InvalidSign() = default;
};

#endif //DI_ZKP_CRYPTO_EXCEPTION_H
