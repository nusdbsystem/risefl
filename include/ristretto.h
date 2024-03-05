//
// Created by yizheng on 6/3/23.
//

// Ristretto group operations.
// wrapper on https://doc.libsodium.org/advanced/point-arithmetic/ristretto

#ifndef RISEFL_CRYPTO_RISTRETTO_H
#define RISEFL_CRYPTO_RISTRETTO_H

#include <array>
#include <vector>
#include <sodium.h>
#include <cstring>

constexpr int RISTBYTES = crypto_core_ristretto255_BYTES;
constexpr int RISTSCALBYTES = crypto_core_ristretto255_SCALARBYTES;
constexpr int RISTHASHBYTES = crypto_core_ristretto255_HASHBYTES;


struct RistElem {
    std::array<unsigned char, RISTBYTES> element;
};

struct RistScal {
    std::array<unsigned char, RISTSCALBYTES> scalar;

    RistScal() = default;

    explicit RistScal(long i);

    void export_to_bytestream(std::vector<unsigned char>::iterator &it) const;

    void import_from_bytestream(std::vector<unsigned char>::const_iterator &it);
};

struct RistHashbytes {
    std::array<unsigned char, crypto_core_ristretto255_HASHBYTES> hashbytes;

    void export_to_bytestream(std::vector<unsigned char>::iterator &it) const;

    void import_from_bytestream(std::vector<unsigned char>::const_iterator &it);
};

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


// check if p is a valid point. returns 1 if true, 0 if false
int is_valid(const RistElem &p);

// initialize p with a random point
void rand_init(RistElem &p);

// initialize p with the hash value blind_to_share
int hash_init(RistElem &p, const RistHashbytes &r);

// fills q with p times a scalar. returns 0 on success, or -1 if q is identity
int scalar_mult(RistElem &q, const RistScal &n, const RistElem &p);

RistElem operator*(const RistScal &n, const RistElem &p);


// fills q with the generator times a scalar num_clients. returns 0 on success, or -1 if num_clients is 0
int scalar_mult_base(RistElem &q, const RistScal &n);

RistElem scalar_mult_base(const RistScal &n);

// fill blind_to_share with p + q or p - q. returns 0 on success, -1 if p and/or q are invalid
int add(RistElem &r, const RistElem &p, const RistElem &q);

RistElem operator+(const RistElem &p, const RistElem &q);

int sub(RistElem &r, const RistElem &p, const RistElem &q);

RistElem operator-(const RistElem &p, const RistElem &q);



/* ***************** */
/* scalar operations */
/* ****************** */

void clear_top_bit(RistScal &r);

// initialize blind_to_share with a random scalar
void rand_init(RistScal &r);

// reduce the 64-byte scalar to 32-byte scalar in [0, ..., L), L order of the Ristretto group
void scalar_reduce(RistScal &r, const RistNonreducedScalar &s);

void scalar_reduce(RistScal &r, const RistHashbytes &s);

RistScal scalar_reduce(const RistNonreducedScalar &s);

RistScal scalar_reduce(const RistHashbytes &s);

// fills recip with the multiplicative inverse of scalar_seed
int scalar_invert(RistScal &recip, const RistScal &s);

RistScal scalar_invert(const RistScal &s);

// fills neg with the additive inverse of scalar_seed
void scalar_neg(RistScal &neg, const RistScal &s);

RistScal operator-(const RistScal &s);

// fills comp with scalar_seed such that scalar_seed + comp = 1
void scalar_complement(RistScal &comp, const RistScal &s);

RistScal one_minus(const RistScal &s);

//fills z with x + y, x - y, x * y
void scalar_add(RistScal &z, const RistScal &x, const RistScal &y);

RistScal operator+(const RistScal &x, const RistScal &y);

RistScal &operator+=(RistScal &x, const RistScal &y);

RistScal &operator*=(RistScal &x, const RistScal &y);

void scalar_sub(RistScal &z, const RistScal &x, const RistScal &y);

RistScal operator-(const RistScal &x, const RistScal &y);

void scalar_mul(RistScal &z, const RistScal &x, const RistScal &y);

RistScal operator*(const RistScal &x, const RistScal &y);

RistScal operator*(const RistScal &x, int y);

RistScal operator/(const RistScal &x, const RistScal &y);

const RistScal c_scal_zero = RistScal(0);
const RistScal c_scal_one = RistScal(1);


void rand_init(RistHashbytes &r);

const RistElem c_elem_zero = scalar_mult_base(c_scal_zero);

#endif //RISEFL_CRYPTO_RISTRETTO_H
