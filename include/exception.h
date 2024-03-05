//
// Created by yizheng on 29/3/23.
//

#ifndef RISEFL_CRYPTO_EXCEPTION_H
#define RISEFL_CRYPTO_EXCEPTION_H

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

#endif //RISEFL_CRYPTO_EXCEPTION_H
