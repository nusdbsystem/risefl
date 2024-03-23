//
// Created by yizheng on 4/4/23.
//

// mpz_conversion.h
// 
// C++ wrapper based on the C-library GMP
// https://gmplib.org/
// Implements C++ classes that manipulates arbitrary length integer arithmetic

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

// An instance of the class MyMpz contains an arbitrary length integer
class MyMpz {
public:

    // the arbitrary length integer
    mpz_t pointer;

    // default constructor: initialize the arbitrary length integer to 0
    MyMpz() { mpz_init(pointer); }

    // destructor
    ~MyMpz() { mpz_clear(pointer); }

    // copy constructor
    MyMpz(const MyMpz &x) { mpz_init_set(pointer, x.pointer); }

    friend void swap(MyMpz &x, MyMpz &y) { std::swap(x.pointer, y.pointer); }

    MyMpz(MyMpz &&x) noexcept { swap(*this, x); }

    // assign operator
    MyMpz &operator=(MyMpz x) {
        swap(*this, x);
        return *this;
    }


    // constructor using a C++ built-in int n as parameter: initialize the arbitrary length integer to n
    MyMpz(int n) {
        mpz_init_set_si(pointer, n);
    }

    // constructor using a string s and an int base: initialize the arbitrary length integer to the number intepreted by s in base-n
    // example: MyMpz("334", 8) initializes the arbitrary length integer to 3 * 64 + 3 * 8 + 4
    MyMpz(std::string s, int base) {
        mpz_init_set_str(pointer, s.data(), base);
    }

    // constructor using a RistScal s: initialize the arbitrary length integer to the 256-bit integer encoded by s
    MyMpz(const RistScal &s) : MyMpz() {
        mpz_import(pointer, 32, -1, 1, 1, 0, s.scalar.data());
    }

    // If a and b are MyMpz instances, a + b returns an MyMpz instance that contains the sum of the ones in a and b 
    MyMpz operator+(const MyMpz &x) const {
        MyMpz res;
        mpz_add(res.pointer, pointer, x.pointer);
        return res;
    }

    // If a and b are MyMpz instances, a + b returns an MyMpz instance that contains the difference of the ones in a and b 
    MyMpz operator-(const MyMpz &x) const {
        MyMpz res;
        mpz_sub(res.pointer, pointer, x.pointer);
        return res;
    }

    // If a and b are MyMpz instances, a + b returns an MyMpz instance that contains the product of the ones in a and b 
    MyMpz operator*(const MyMpz &x) const {
        MyMpz res;
        mpz_mul(res.pointer, pointer, x.pointer);
        return res;
    }

    // If a and b are MyMpz instances, a += b adds the arbitrary length integer in b to the one in a
    void operator+=(const MyMpz &x) {
        mpz_add(pointer, pointer, x.pointer);
    }

//    MyMpz remainder(const MyMpz &x) const{
//        MyMpz res;
//        mpz_mod(res.pointer, pointer, x.pointer);
//        return res;
//    }

    // changes the arbitrary length integer (say n) to (n mod r), where r is the order of the Ristretto group 7237005577332262213973186563042994240857116359379907606001950938285454250989
    // as a result, the arbitrary length integer must be between 0 and r-1 (both inclusive)
    void reduce_to_rist();
};

// the order of the Ristretto group
static const MyMpz c_mpz_rist_order("7237005577332262213973186563042994240857116359379907606001950938285454250989", 10);

// convert a RistScal to an MyMpz
void mpz_from_rist(MyMpz &x, const RistScal &s);

// convert a C++ built-in int to an MyMpz
void mpz_from_int(MyMpz &x, int n);

// convert an MyMpz x to a RistScal (x mod r), where r is the order of the Ristretto group 7237005577332262213973186563042994240857116359379907606001950938285454250989
void rist_from_mpz(RistScal &s, const MyMpz &x);

using MpzVec = std::vector<MyMpz>;
using MpzMat = std::vector<MpzVec>;

// converts a vector of C++ built-in ints to a vector of MyMpz's
void mpz_vec_from_int(MpzVec &xx, const std::vector<int> &nn);

// converts a vector of RistScal's to a vector of MyMpz's
void mpz_vec_from_rist(MpzVec &xx, const RistScalVec &ss);

// converts a matrix (double vector) of RistScal's to a matrix (double vector) of MyMpz's
void mpz_mat_from_rist(MpzMat &xx, const RistScalMat &ss);

// converts a vector of MyMpz's to a vector of RistScal's
void rist_from_mpz_vec(RistScalVec &ss, const MpzVec &xx);

// converts a matrix (double vector) of MyMpz's to a matrix (double vector) of RistScal's
void rist_from_mpz_mat(RistScalMat &ss, const MpzMat &xx);


#endif //RISEFL_CRYPTO_MPZ_CONVERSION_H
