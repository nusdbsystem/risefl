//
// Created by yizheng on 4/4/23.
//

#ifndef RISEFL_CRYPTO_MPZ_CONVERSION_H
#define RISEFL_CRYPTO_MPZ_CONVERSION_H

#include <iostream>
#include <chrono>
#include <cassert>
#include <vector>
#include <string>
#include <algorithm>
#include <execution>
#include <thread>
#include <sodium.h>

//#include <gmpxx.h>
#include <gmp.h>
#include "ristretto.h"
#include "ristretto_vector.h"

class MyMpz {
public:
    mpz_t pointer;

    MyMpz() { mpz_init(pointer); }

    ~MyMpz() { mpz_clear(pointer); }

    MyMpz(const MyMpz &x) { mpz_init_set(pointer, x.pointer); }

    friend void swap(MyMpz &x, MyMpz &y) { std::swap(x.pointer, y.pointer); }

    MyMpz(MyMpz &&x) noexcept { swap(*this, x); }

    MyMpz &operator=(MyMpz x) {
        swap(*this, x);
        return *this;
    }


    MyMpz(int n) {
        mpz_init_set_si(pointer, n);
    }

    MyMpz(std::string s, int base) {
        mpz_init_set_str(pointer, s.data(), base);
    }

    MyMpz(const RistScal &s) : MyMpz() {
        mpz_import(pointer, 32, -1, 1, 1, 0, s.scalar.data());
    }

    MyMpz operator+(const MyMpz &x) const {
        MyMpz res;
        mpz_add(res.pointer, pointer, x.pointer);
        return res;
    }

    MyMpz operator-(const MyMpz &x) const {
        MyMpz res;
        mpz_sub(res.pointer, pointer, x.pointer);
        return res;
    }

    MyMpz operator*(const MyMpz &x) const {
        MyMpz res;
        mpz_mul(res.pointer, pointer, x.pointer);
        return res;
    }

    void operator+=(const MyMpz &x) {
        mpz_add(pointer, pointer, x.pointer);
    }

//    MyMpz remainder(const MyMpz &x) const{
//        MyMpz res;
//        mpz_mod(res.pointer, pointer, x.pointer);
//        return res;
//    }
    void reduce_to_rist();
};

static const MyMpz c_mpz_rist_order("7237005577332262213973186563042994240857116359379907606001950938285454250989", 10);

void mpz_from_rist(MyMpz &x, const RistScal &s);

void mpz_from_int(MyMpz &x, int n);

void rist_from_mpz(RistScal &s, const MyMpz &x);

using MpzVec = std::vector<MyMpz>;
using MpzMat = std::vector<MpzVec>;

void mpz_vec_from_int(MpzVec &xx, const std::vector<int> &nn);

//void mpz_mat_from_int(MpzMat& xx, const std::vector<std::vector<int>> &nn);
void mpz_vec_from_rist(MpzVec &xx, const RistScalVec &ss);

void mpz_mat_from_rist(MpzMat &xx, const RistScalMat &ss);

void rist_from_mpz_vec(RistScalVec &ss, const MpzVec &xx);

void rist_from_mpz_mat(RistScalMat &ss, const MpzMat &xx);


#endif //RISEFL_CRYPTO_MPZ_CONVERSION_H
