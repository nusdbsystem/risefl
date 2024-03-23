//
// Created by yizheng on 6/3/23.
//

// ristretto.h
// 
// Ristretto group operations.
// C++ wrapper on https://doc.libsodium.org/advanced/point-arithmetic/ristretto

#ifndef RISEFL_CRYPTO_RISTRETTO_H
#define RISEFL_CRYPTO_RISTRETTO_H

#include <array>
#include <vector>
#include <sodium.h>
#include <cstring>

// the length of an element of the Ristretto group in bytes (which is 32)
constexpr int RISTBYTES = crypto_core_ristretto255_BYTES;
// the length of Ristretto scalar (i.e. a nonnegative number < the order of the Ristretto group) in bytes (which is 32)
constexpr int RISTSCALBYTES = crypto_core_ristretto255_SCALARBYTES;
// the length of the Ristretto hash in bytes (which is 64)
constexpr int RISTHASHBYTES = crypto_core_ristretto255_HASHBYTES;

// Ristretto group element
struct RistElem {
    std::array<unsigned char, RISTBYTES> element;
};

// Ristretto scalar (i.e. a nonnegative number < the order of the Ristretto group)
struct RistScal {
    std::array<unsigned char, RISTSCALBYTES> scalar;

    RistScal() = default;

    // initialize the number to a C++ built-in long int
    explicit RistScal(long i);

    // convert the scalar to bytes, write the bytes to a vector starting from the place pointed by the iterator it, and move the iterator it to one place after the end
    void export_to_bytestream(std::vector<unsigned char>::iterator &it) const;

    // read bytes from a vector starting from the place pointed by the iterator it, convert the bytes to the scalar, and move the iterator it to one place after the end
    void import_from_bytestream(std::vector<unsigned char>::const_iterator &it);
};

// Ristretto hash
struct RistHashbytes {
    std::array<unsigned char, crypto_core_ristretto255_HASHBYTES> hashbytes;

    // convert the hash to bytes, write the bytes to a vector starting from the place pointed by the iterator it, and move the iterator it to one place after the end
    void export_to_bytestream(std::vector<unsigned char>::iterator &it) const;

    // read bytes from a vector starting from the place pointed by the iterator it, convert the bytes to the hash, and move the iterator it to one place after the end
    void import_from_bytestream(std::vector<unsigned char>::const_iterator &it);
};

// nonreduced Ristretto scalar (in 64 bytes), never used
struct RistNonreducedScalar {
    std::array<unsigned char, crypto_core_ristretto255_NONREDUCEDSCALARBYTES> nonreduced_scalar;
};

inline bool operator==(const RistElem &p, const RistElem &q) {
    return p.element == q.element;
}

inline bool operator==(const RistScal &p, const RistScal &q) {
    return p.scalar == q.scalar;
}

inline bool operator==(const RistHashbytes &p, const RistHashbytes &q) {
    return p.hashbytes == q.hashbytes;
}

inline bool operator==(const RistNonreducedScalar &p, const RistNonreducedScalar &q) {
    return p.nonreduced_scalar == q.nonreduced_scalar;
}

inline bool operator!=(const RistElem &p, const RistElem &q) {
    return p.element != q.element;
}

inline bool operator!=(const RistScal &p, const RistScal &q) {
    return p.scalar != q.scalar;
}

inline bool operator!=(const RistHashbytes &p, const RistHashbytes &q) {
    return p.hashbytes != q.hashbytes;
}

inline bool operator!=(const RistNonreducedScalar &p, const RistNonreducedScalar &q) {
    return p.nonreduced_scalar != q.nonreduced_scalar;
}

/* ***************** */
/* group operations */
/* ****************** */


// check if p is a valid point in the Ristretto group. returns 1 if true, 0 if false
int is_valid(const RistElem &p);

// initialize p with a random point in the Ristretto group
void rand_init(RistElem &p);

// initialize p with the hash value r
int hash_init(RistElem &p, const RistHashbytes &r);

// fills q with n * p. returns 0 on success, or -1 if q is identity
int scalar_mult(RistElem &q, const RistScal &n, const RistElem &p);

// returns n * p
RistElem operator*(const RistScal &n, const RistElem &p);

// fills q with n * (the fixed generator). returns 0 on success, or -1 if n is 0
int scalar_mult_base(RistElem &q, const RistScal &n);

// returns n * (the fixed generator)
RistElem scalar_mult_base(const RistScal &n);

// fill r with p + q. returns 0 on success, -1 if p and/or q are invalid
int add(RistElem &r, const RistElem &p, const RistElem &q);

// returns p + q
RistElem operator+(const RistElem &p, const RistElem &q);

// fill r with p - q. returns 0 on success, -1 if p and/or q are invalid
int sub(RistElem &r, const RistElem &p, const RistElem &q);

// return p - q
RistElem operator-(const RistElem &p, const RistElem &q);



/* ***************** */
/* scalar operations */
/* ****************** */

// helper function. the top bit of a RistScal should always be 0. used in carrying out arithmetic
void clear_top_bit(RistScal &r);

// initializes r with a random scalar
void rand_init(RistScal &r);

// reduces the 64-byte scalar to 32-byte scalar in [0, ..., L), L order of the Ristretto group
void scalar_reduce(RistScal &r, const RistNonreducedScalar &s);

// reduces a 64-byte hash to a RistScal
void scalar_reduce(RistScal &r, const RistHashbytes &s);

// reduces a 64-byte number to a RistScal (never used)
RistScal scalar_reduce(const RistNonreducedScalar &s);

// reduces a 64-byte hash s to a Ristscal
RistScal scalar_reduce(const RistHashbytes &s);

// fills recip with the multiplicative inverse of s
int scalar_invert(RistScal &recip, const RistScal &s);

// computes the multiplicative inverse of s
RistScal scalar_invert(const RistScal &s);

// fills neg with the additive inverse of s
void scalar_neg(RistScal &neg, const RistScal &s);

// computes the additive inverse of s
RistScal operator-(const RistScal &s);

// fills comp with s such that s + comp = 1
void scalar_complement(RistScal &comp, const RistScal &s);

// computes 1 - s
RistScal one_minus(const RistScal &s);

//fills z with x + y
void scalar_add(RistScal &z, const RistScal &x, const RistScal &y);

RistScal operator+(const RistScal &x, const RistScal &y);

RistScal &operator+=(RistScal &x, const RistScal &y);

RistScal &operator*=(RistScal &x, const RistScal &y);

// fills z with x - y
void scalar_sub(RistScal &z, const RistScal &x, const RistScal &y);

RistScal operator-(const RistScal &x, const RistScal &y);

// fills z with x * y
void scalar_mul(RistScal &z, const RistScal &x, const RistScal &y);

RistScal operator*(const RistScal &x, const RistScal &y);

RistScal operator*(const RistScal &x, int y);

RistScal operator/(const RistScal &x, const RistScal &y);

// the scalar 0
const RistScal c_scal_zero = RistScal(0);
// the scalar 1
const RistScal c_scal_one = RistScal(1);

// randomize a 64-byte hash r
void rand_init(RistHashbytes &r);

// the identity Ristretto group element
const RistElem c_elem_zero = scalar_mult_base(c_scal_zero);

#endif //RISEFL_CRYPTO_RISTRETTO_H
