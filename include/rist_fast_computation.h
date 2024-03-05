//
// Created by yizheng on 26/3/23.
//

#ifndef RISEFL_CRYPTO_RIST_FAST_COMPUTATION_H
#define RISEFL_CRYPTO_RIST_FAST_COMPUTATION_H


//#include "sodium/ed25519_ref10.h"
#include "sodium/ed25519_ref10.h"

#include <array>
#include <cassert>
#include <vector>
#include "ristretto.h"
#include "ristretto_vector.h"
#include "utils.h"

struct RistElemP1 {
    ge25519_p1p1 element;
};
struct RistElemP2 {
    ge25519_p2 element;
};

struct RistElemP3 {
    ge25519_p3 element;

    RistElemP3() = default;

    RistElemP3 &operator=(RistElemP3 other);

    RistElemP3(const RistElemP3 &other);

    RistElemP3(RistElemP3 &&other) noexcept;

    bool operator==(const RistElemP3 &other) const;

    bool operator!=(const RistElemP3 &other) const {
        return !(*this == other);
    }
};

void swap(RistElemP3 &p, RistElemP3 &q);

struct RistElemCached {
    ge25519_cached element;
};

struct RistElemPrecomp {
    ge25519_precomp element;
};

using RistElemPrecompTable = std::array<ge25519_precomp[8], 32>;

using RistElemP3Vec = std::vector<RistElemP3>;
using RistElemP3Mat = std::vector<RistElemP3Vec>;

int p3_from_bytes(RistElemP3 &p, const std::array<unsigned char, RISTBYTES> &s);

int p3_from_bytes(RistElemP3 &p, const std::vector<unsigned char> &s);

int p3_from_bytestream(RistElemP3 &p, std::vector<unsigned char>::const_iterator &it_bytes);

void p3_vec_from_bytestream(RistElemP3Vec &rr, std::vector<unsigned char>::const_iterator &it_bytes);

void p3_vec_from_bytes(RistElemP3Vec &rr, const std::vector<unsigned char> &pp);

void bytes_from_p3(std::array<unsigned char, RISTBYTES> &s, const RistElemP3 &p);

void bytes_from_p3(std::vector<unsigned char> &s, const RistElemP3 &p);

void bytestream_from_p3(std::vector<unsigned char>::iterator &it_bytes, const RistElemP3 &p);

std::array<unsigned char, RISTBYTES> bytes_from_p3(const RistElemP3 &p);

void bytestream_from_p3_vec(std::vector<unsigned char>::iterator &it_bytes, const RistElemP3Vec &pp);

std::vector<unsigned char> bytes_from_p3_vec(const RistElemP3Vec &pp);

void rist_elem_vec_from_p3_vec(RistElemVec &rr, RistElemP3Vec &pp);

void p3_copy(RistElemP3 &r, const RistElemP3 &p);
//void p3_copy(RistElemP3Vec& r, const RistElemP3Vec& p);

inline void p1_to_p2(RistElemP2 &r, const RistElemP1 &p) {
    ge25519_p1p1_to_p2(&r.element, &p.element);
}

inline void p1_to_p3(RistElemP3 &r, const RistElemP1 &p) {
    ge25519_p1p1_to_p3(&r.element, &p.element);
}

inline void p2_to_p3(RistElemP3 &r, const RistElemP2 &p) {
    ge25519_p2_to_p3(&r.element, &p.element);
}

inline void p3_to_p2(RistElemP2 &r, const RistElemP3 &p) {
    ge25519_p3_to_p2(&r.element, &p.element);
}

inline void p3_to_cached(RistElemCached &r, const RistElemP3 &p) {
    ge25519_p3_to_cached(&r.element, &p.element);
}

//inline void add(RistElemP1& r, const RistElemP3& p, const RistElemCached& q) {
//    ge25519_add(&r, &p, &q);
//}
//inline void sub(RistElemP1& r, const RistElemP3& p, const RistElemCached& q) {
//    ge25519_sub(&r, &p, &q);
//}

inline void p3_add(RistElemP3 &r, const RistElemP3 &p, const RistElemP3 &q) {
    ge25519_p3_add(&r.element, &p.element, &q.element);
}

inline void p3_add(RistElemP3 &r, const RistElemP3 &p, const RistElemCached &q) {
    RistElemP1 p1;
    ge25519_add(&p1.element, &p.element, &q.element);
    ge25519_p1p1_to_p3(&r.element, &p1.element);
}

inline void p3_sub(RistElemP3 &r, const RistElemP3 &p, const RistElemP3 &q) {
    ge25519_p3_sub(&r.element, &p.element, &q.element);
}

inline void p3_sub(RistElemP3 &r, const RistElemP3 &p, const RistElemCached &q) {
    RistElemP1 p1;
    ge25519_sub(&p1.element, &p.element, &q.element);
    ge25519_p1p1_to_p3(&r.element, &p1.element);
}

inline void p3_0(RistElemP3 &r) {
    ge25519_p3_0(&r.element);
}

inline bool is_0(const RistElemP3 &r) {
    RistElemP3 s;
    p3_0(s);
    return r == s;
}

inline RistElemP3 operator+(const RistElemP3 &p, const RistElemP3 &q) {
    RistElemP3 r;
    p3_add(r, p, q);
    return r;
}

inline RistElemP3 &operator+=(RistElemP3 &p, const RistElemP3 &q) {
    p3_add(p, p, q);
    return p;
}

inline RistElemP3 operator-(const RistElemP3 &p, const RistElemP3 &q) {
    RistElemP3 r;
    p3_sub(r, p, q);
    return r;
}

inline RistElemP3 p3_0() {
    RistElemP3 r;
    p3_0(r);
    return r;
}

/* r = r*(2^n)+q */
inline void power_2_add(RistElemP3 &r, int n, const RistElemP3 &q) {
    return ge25519_p3_dbladd(&r.element, n, &q.element);
}

void pedersen_commit(RistElemP3 &y, const RistScal &x, const RistElemP3 &h, const RistScal &r);

void pedersen_commit(RistElemP3 &y, long x, int x_bits, const RistElemP3 &h, const RistScal &r);

RistElemP3 pedersen_commit(const RistScal &x, const RistElemP3 &h, const RistScal &r);

void pedersen_commit_precomp(RistElemP3 &y, long x, int x_bits, const RistElemPrecompTable &table, const RistScal &r);

void pedersen_commit_precomp(RistElemP3 &y, const RistScal &x, const RistElemPrecompTable &table, const RistScal &r);

void pedersen_zero_commit(RistElemP3 &y, const RistScal &x);

void pedersen_zero_commit(RistElemP3 &y, long x);

RistElemP3 pedersen_zero_commit_p3(const RistScal &x);

RistElemP3 pedersen_zero_commit_p3(long x);

void pedersen_commit(RistElemP3Vec &yy, const RistScalVec &xx, const RistElemP3Vec &hh, const RistScal &r);

void pedersen_commit(RistElemP3Vec &yy, const RistScalVec &xx, const RistElemP3 &h, const RistScalVec &rr);

void pedersen_zero_commit(RistElemP3Vec &yy, const RistScalVec &xx);

RistElemP3Vec pedersen_commit(const RistScalVec &xx, const RistElemP3Vec &hh, const RistScal &r);

void pedersen_commit_p3_to_bytes(std::vector<unsigned char>::iterator &it_bytes, const RistScalVec &xx,
                                 const RistElemP3Vec &hh, const RistScal &r);

void
pedersen_commit_p3_to_bytes(std::vector<unsigned char>::iterator &it_bytes, const std::vector<long> &xx, int x_bits,
                            const RistElemP3Vec &hh, const RistScal &r);

void
pedersen_commit_p3_to_bytes_from_precomp(std::vector<unsigned char>::iterator &it_bytes, const std::vector<long> &xx,
                                         int xx_bits,
                                         const std::vector<RistElemPrecompTable> &tables,
                                         const RistScalVec &rr);

//std::vector<unsigned char> pedersen_commit_p3_to_bytes_from_precomp(const RistScalVec& xx, const RistElemP3Vec& hh, const RistScal& r);
void scalar_mult(RistElemP3 &r, const RistScal &n, const RistElemP3 &p);

RistElemP3 operator*(const RistScal &n, const RistElemP3 &p);

void scalar_mult(RistElemP3Vec &rr, const RistScal &n, const RistElemP3Vec &pp);

RistElemP3Vec operator*(const RistScal &n, const RistElemP3Vec &pp);

//void p3_add(RistElemP3Vec &rr, const RistElemP3Vec &pp, const RistElemP3Vec &qq);

void p3_add(RistElemP3Vec &rr, const RistElemP3Vec &pp, const RistElemP3Vec &qq);

RistElemP3Vec operator+(const RistElemP3Vec &pp, const RistElemP3Vec &qq);

void sum_single_thread(RistElemP3 &s, const RistElemP3Vec &vv);

void sum_single_thread(RistElemP3 &s, const RistElemP3Mat &vv);

//RistElemP3 sum_parallel(const RistElemP3Vec& vv);

void linear_comb(RistElemP3 &s, const RistScalVec &rr, const RistElemP3Vec &pp);

void linear_comb(RistElemP3 &s, std::vector<int> &rr, int bit_bound, const RistElemP3Vec &pp);

void linear_comb(RistElemP3 &s, const RistScalMat &rr, const RistElemP3Mat &pp);

RistElemP3 linear_comb(const RistScalVec &rr, const RistElemP3Vec &pp);

RistElemP3 linear_comb(const RistScalMat &rr, const RistElemP3Mat &pp);

class LinearCombCalculator {
public:
    RistScalVec coeffs;
    RistElemP3Vec group_elements;

    LinearCombCalculator() = default;

    void add(const RistScal &a, const RistElemP3 &r) {
        coeffs.emplace_back(a);
        group_elements.emplace_back(r);
    }

    void add(const RistElemP3 &r) {
        add(c_scal_one, r);
    }

    void sub(const RistScal &a, const RistElemP3 &r) {
        add(-a, r);
    }

    void sub(const RistElemP3 &r) {
        sub(c_scal_one, r);
    }

    void add_base(const RistScal &a) {
        add(a, pedersen_zero_commit_p3(1));
    }

    void sub_base(const RistScal &a) {
        sub(a, pedersen_zero_commit_p3(1));
    }

    void add(const RistScalVec &aa, const RistElemP3Vec &rr) {
        coeffs.insert(coeffs.end(), aa.begin(), aa.end());
        group_elements.insert(group_elements.end(), rr.begin(), rr.end());
    }

    void sub(const RistScalVec &bb, const RistElemP3Vec &rr) {
        add(-bb, rr);
    }

    void add(const RistScalMat &aa, const RistElemP3Mat &rr) {
        assert(aa.size() == rr.size());
        for (int i = 0; i < aa.size(); i++)
            coeffs.insert(coeffs.end(), aa[i].begin(), aa[i].end());
        for (int i = 0; i < rr.size(); i++)
            group_elements.insert(group_elements.end(), rr[i].begin(), rr[i].end());
    }

    void sub(const RistScalMat &aa, const RistElemP3Mat &rr) {
        add(-aa, rr);
    }

    RistElemP3 result() const {
        RistElemP3 res;
        linear_comb(res, coeffs, group_elements);
        return res;
    }

    bool is_0() const {
        return ::is_0(result());
    };
};

// ss[0] = rr[0] * pp[0] + rr[n] * pp[1] + rr[2 * n] * pp[2] + ...
// ss[1] = rr[1] * pp[0] + rr[n + 1] * pp[1] + rr[2 * n + 1] * pp[2] + ...
// ...
// ss[n - 1] = rr[n - 1] * pp[0] + rr[2 * n - 1] * pp[1] + rr[3 * n - 1] * pp[2] + ...
void linear_comb_block(RistElemP3Vec &ss, const RistScalVec &rr, const RistElemP3Vec &pp);

void linear_comb_block(RistElemP3Vec &ss, const std::vector<int> &rr, int bit_bound, const RistElemP3Vec &pp);

void p3_from_hash(RistElemP3 &r, const RistHashbytes &h);

RistElemP3 p3_from_hash(const RistHashbytes &h);

template<typename... Ts>
RistElemP3 p3_from_hash_from_bytes(const Ts &... ts) {
    return p3_from_hash(hash_from_bytes(std::forward<const Ts>(ts)...));
}

//void p3_to_precomp(RistElemPrecomp& s, const RistElemP3& r) {
//    ge25519_p3_to_precomp(&s.element, &r.element);
//}

void generate_precomp_table(RistElemPrecompTable &table, const RistElemP3 &r);

void scalar_mult_precomp(RistElemP3 &r, const RistScal &n, const RistElemPrecompTable &table);

void
scalar_double_mult_precomp(RistElemP3 &r, long x, int x_bits, const RistScal &n, const RistElemPrecompTable &table);

struct RistP3AndBytes {
    RistElemP3 elem;
    std::vector<unsigned char> bytes;

    RistP3AndBytes() : bytes(RISTBYTES) {}

    void fill_bytes() {
        bytes_from_p3(bytes, elem);
    }

    void fill_elem() {
        p3_from_bytes(elem, bytes);
    }

    void export_to_bytestream(std::vector<unsigned char>::iterator &it) const {
        copy_and_shift_dst(bytes.begin(), bytes.end(), it);
    }

    void import_from_bytestream(std::vector<unsigned char>::const_iterator &it) {
        copy_and_shift_src(it, bytes.begin(), bytes.end());
        fill_elem();
    }

    RistP3AndBytes operator+(const RistElemP3 &r) const {
        RistP3AndBytes ret;
        ret.elem = elem + r;
        ret.fill_bytes();
        return ret;
    }
};

struct RistP3VecAndBytes {
    RistElemP3Vec elems;
    std::vector<unsigned char> bytes;

    RistP3VecAndBytes(int n = 0) : elems(n), bytes(n * RISTBYTES) {}

    void fill_bytes() {
        auto it_bytes = bytes.begin();
        bytestream_from_p3_vec(it_bytes, elems);
        assert(it_bytes == bytes.end());
    }

    void fill_bytes(int i) {
        auto it_bytes = bytes.begin() + i * RISTBYTES;
        bytestream_from_p3(it_bytes, elems[i]);
    }

    void fill_elems() {
        auto it_bytes = bytes.cbegin();
        p3_vec_from_bytestream(elems, it_bytes);
        assert(it_bytes == bytes.end());
    }

    void export_to_bytestream(std::vector<unsigned char>::iterator &it) const {
        copy_and_shift_dst(bytes.begin(), bytes.end(), it);
    }

    void import_from_bytestream(std::vector<unsigned char>::const_iterator &it) {
        copy_and_shift_src(it, bytes.begin(), bytes.end());
        fill_elems();
    }

    auto size() const -> decltype(elems.size()) {
        return elems.size();
    }

    void resize(int n) {
        elems.resize(n);
        bytes.resize(n * RISTBYTES);
    }

    RistP3VecAndBytes segment(int n) const {
        RistP3VecAndBytes s(n);
        std::copy_n(elems.begin(), n, s.elems.begin());
        std::copy_n(bytes.begin(), n * RISTBYTES, s.bytes.begin());
        return s;
    }

    RistP3AndBytes operator[](int i) const {
        RistP3AndBytes ret;
        ret.elem = elems[i];
        std::copy_n(bytes.begin() + i * RISTBYTES, RISTBYTES, ret.bytes.begin());
        return ret;
    }

    void emplace_back(const RistP3AndBytes &r) {
        elems.emplace_back(r.elem);
        bytes.insert(bytes.end(), r.bytes.begin(), r.bytes.end());
    }

    RistP3VecAndBytes no_head() {
        RistP3VecAndBytes ret;
        ret.elems = RistElemP3Vec(elems.begin() + 1, elems.end());
        ret.bytes = std::vector<unsigned char>(bytes.begin() + RISTBYTES, bytes.end());
        return ret;
    }

    RistP3VecAndBytes no_head_no_tail() {
        RistP3VecAndBytes ret;
        ret.elems = RistElemP3Vec(elems.begin() + 1, elems.end() - 1);
        ret.bytes = std::vector<unsigned char>(bytes.begin() + RISTBYTES, bytes.end() - RISTBYTES);
        return ret;
    }
};

class RistP3MatAndBytes {
public:
    std::vector<RistElemP3Vec> elems;
    std::vector<unsigned char> bytes;

    RistP3MatAndBytes() = default;

    RistP3MatAndBytes(int m, int n) :
            elems(m, RistElemP3Vec(n)),
            bytes(m * n * RISTBYTES) {}

    void fill_bytes() {
        auto it_bytes = bytes.begin();
        for (int i = 0; i < elems.size(); i++) {
            bytestream_from_p3_vec(it_bytes, elems[i]);
        }
        assert(it_bytes == bytes.end());
    }

    void fill_elems() {
        auto it_bytes = bytes.cbegin();
        for (int i = 0; i < elems.size(); i++) {
            p3_vec_from_bytestream(elems[i], it_bytes);
        }
        assert(it_bytes == bytes.end());
    }

    void export_to_bytestream(std::vector<unsigned char>::iterator &it) const {
        copy_and_shift_dst(bytes.begin(), bytes.end(), it);
    }

    void import_from_bytestream(std::vector<unsigned char>::const_iterator &it) {
        copy_and_shift_src(it, bytes.begin(), bytes.end());
        fill_elems();
    }

    auto size_rows() const -> decltype(elems.size()) {
        return elems.size();
    }

    auto size_cols() const -> decltype(elems[0].size()) {
        return elems[0].size();
    }

    RistP3VecAndBytes slice_col_0() {
        RistP3VecAndBytes ret(size_rows());
        for (int i = 0; i < ret.size(); i++) {
            ret.elems[i] = elems[i][0];
            std::copy_n(bytes.begin() + i * size_cols() * RISTBYTES,
                        RISTBYTES,
                        ret.bytes.begin() + i * RISTBYTES
            );
        }
        return ret;
    }

    RistP3VecAndBytes operator[](int i) {
        RistP3VecAndBytes ret(size_cols());
        ret.elems = elems[i];
        std::copy_n(bytes.begin() + i * size_cols() * RISTBYTES, RISTBYTES * size_cols(), ret.bytes.begin());
        return ret;
    }
};


RistP3AndBytes zero_rist_p3_bytes();

RistP3AndBytes one_rist_p3_bytes();

const RistP3AndBytes c_p3_bytes_zero = zero_rist_p3_bytes();
const RistP3AndBytes c_p3_bytes_one = one_rist_p3_bytes();

RistElemP3 linear_comb_trim(const RistScalVec &rr, const RistElemP3Vec &ss);

RistElemP3Vec operator*(const RistScalVec &nn, const RistElemP3Vec &pp);

RistElemP3 sum_single_thread(const RistElemP3Vec &vv);

void rand_init(RistP3AndBytes &h);

#endif //RISEFL_CRYPTO_RIST_FAST_COMPUTATION_H
