//
// Created by yizheng on 14/3/23.
//
#include "include/ristretto.h"
#include "include/utils.h"
#include <cstring>
#include <memory>

RistScal::RistScal(long c) {
    if (c >= 0) {
        memcpy(scalar.data(), &c, sizeof c);
        memset(scalar.data() + sizeof c, 0, scalar.size() - sizeof c);
    } else {
        *this = -RistScal(-c);
    }
}

void RistScal::export_to_bytestream(std::vector<unsigned char>::iterator &it) const {
    copy_and_shift_dst(scalar.begin(), scalar.end(), it);
}

void RistScal::import_from_bytestream(std::vector<unsigned char>::const_iterator &it) {
    copy_and_shift_src(it, scalar.begin(), scalar.end());
}

void RistHashbytes::export_to_bytestream(std::vector<unsigned char>::iterator &it) const {
    copy_and_shift_dst(hashbytes.begin(), hashbytes.end(), it);
}

void RistHashbytes::import_from_bytestream(std::vector<unsigned char>::const_iterator &it) {
    copy_and_shift_src(it, hashbytes.begin(), hashbytes.end());
}

int is_valid(const RistElem &p) {
    return crypto_core_ristretto255_is_valid_point(p.element.data());
}

void rand_init(RistElem &p) {
    crypto_core_ristretto255_random(p.element.data());
}

int hash_init(RistElem &p, const RistHashbytes &r) {
    return crypto_core_ristretto255_from_hash(p.element.data(), r.hashbytes.data());
}

int scalar_mult(RistElem &q, const RistScal &n, const RistElem &p) {
    return crypto_scalarmult_ristretto255(q.element.data(), n.scalar.data(), p.element.data());
}

RistElem operator*(const RistScal &n, const RistElem &p) {
    RistElem q;
    scalar_mult(q, n, p);
    return q;
}

int scalar_mult_base(RistElem &q, const RistScal &n) {
    return crypto_scalarmult_ristretto255_base(q.element.data(), n.scalar.data());
}

RistElem scalar_mult_base(const RistScal &n) {
    RistElem q;
    scalar_mult_base(q, n);
    return q;
}

int add(RistElem &r, const RistElem &p, const RistElem &q) {
    return crypto_core_ristretto255_add(r.element.data(), p.element.data(), q.element.data());
}

RistElem operator+(const RistElem &p, const RistElem &q) {
    RistElem r;
    add(r, p, q);
    return r;
}

int sub(RistElem &r, const RistElem &p, const RistElem &q) {
    return crypto_core_ristretto255_sub(r.element.data(), p.element.data(), q.element.data());
}

RistElem operator-(const RistElem &p, const RistElem &q) {
    RistElem r;
    sub(r, p, q);
    return r;
}

void clear_top_bit(RistScal &r) {
    r.scalar[31] &= 127;
}

void rand_init(RistScal &r) {
    crypto_core_ristretto255_scalar_random(r.scalar.data());
    clear_top_bit(r);
}

void scalar_reduce(RistScal &r, const RistNonreducedScalar &s) {
    crypto_core_ristretto255_scalar_reduce(r.scalar.data(), s.nonreduced_scalar.data());
    clear_top_bit(r);
}

void scalar_reduce(RistScal &r, const RistHashbytes &s) {
    crypto_core_ristretto255_scalar_reduce(r.scalar.data(),
                                           reinterpret_cast<const RistNonreducedScalar *>(&s)->nonreduced_scalar.data());
    clear_top_bit(r);
}

RistScal scalar_reduce(const RistNonreducedScalar &s) {
    RistScal r;
    scalar_reduce(r, s);
    return r;
}

RistScal scalar_reduce(const RistHashbytes &s) {
    return scalar_reduce(*reinterpret_cast<const RistNonreducedScalar *>(&s));
}

int scalar_invert(RistScal &recip, const RistScal &s) {
    int i = crypto_core_ristretto255_scalar_invert(recip.scalar.data(), s.scalar.data());
    clear_top_bit(recip);
    return i;
}

RistScal scalar_invert(const RistScal &s) {
    RistScal recip;
    scalar_invert(recip, s);
    return recip;
}

void scalar_neg(RistScal &neg, const RistScal &s) {
    crypto_core_ristretto255_scalar_negate(neg.scalar.data(), s.scalar.data());
    clear_top_bit(neg);
}

RistScal operator-(const RistScal &s) {
    RistScal neg;
    scalar_neg(neg, s);
    return neg;
}

void scalar_complement(RistScal &comp, const RistScal &s) {
    crypto_core_ristretto255_scalar_complement(comp.scalar.data(), s.scalar.data());
    clear_top_bit(comp);
}

RistScal one_minus(const RistScal &s) {
    RistScal comp;
    scalar_complement(comp, s);
    return comp;
}

void scalar_add(RistScal &z, const RistScal &x, const RistScal &y) {
    crypto_core_ristretto255_scalar_add(z.scalar.data(), x.scalar.data(), y.scalar.data());
    clear_top_bit(z);
}

RistScal operator+(const RistScal &x, const RistScal &y) {
    RistScal z;
    scalar_add(z, x, y);
    return z;
}

RistScal &operator+=(RistScal &x, const RistScal &y) {
    scalar_add(x, x, y);
    return x;
}

RistScal &operator*=(RistScal &x, const RistScal &y) {
    scalar_mul(x, x, y);
    return x;
}

void scalar_sub(RistScal &z, const RistScal &x, const RistScal &y) {
    crypto_core_ristretto255_scalar_sub(z.scalar.data(), x.scalar.data(), y.scalar.data());
    clear_top_bit(z);
}

RistScal operator-(const RistScal &x, const RistScal &y) {
    RistScal z;
    scalar_sub(z, x, y);
    return z;
}

void scalar_mul(RistScal &z, const RistScal &x, const RistScal &y) {
    crypto_core_ristretto255_scalar_mul(z.scalar.data(), x.scalar.data(), y.scalar.data());
    clear_top_bit(z);
}

RistScal operator*(const RistScal &x, const RistScal &y) {
    RistScal z;
    scalar_mul(z, x, y);
    return z;
}

RistScal operator*(const RistScal &x, int y) {
    return x * RistScal(y);
}

RistScal operator/(const RistScal &x, const RistScal &y) {
    return x * scalar_invert(y);
}

void rand_init(RistHashbytes &r) {
    randombytes_buf(r.hashbytes.data(), r.hashbytes.size());
}
