//
// Created by yizheng on 8/3/23.
//

// zkp.h
//
// ZKP toolkit 

#ifndef RISEFL_CRYPTO_ZKP_H
#define RISEFL_CRYPTO_ZKP_H

//#include <NTL/ZZ.h>
//#include <NTL/ZZ_p.h>
//#include <NTL/ZZ_pX.h>
#include <chrono>
#include <cassert>
#include <vector>

#include "ristretto.h"
#include "ristretto_vector.h"
#include "zkp_hash.h"
#include "rist_fast_computation.h"
#include "utils.h"

// compute the pedersen commitment x * g + r * h, where g is the fixed group generator
inline RistElem pedersen_commit(const RistScal &x, const RistElem &h, const RistScal &r) {
    return scalar_mult_base(x) + r * h;
}

// compute the pedersen commitment x * g + r * h, and write the result into y
inline void pedersen_commit(RistElem &y, const RistScal &x, const RistElem &h, const RistScal &r) {
    add(y, scalar_mult_base(x), r * h);
}

// compute the pedersen commitment x * g, where g is the fixed group generator
inline RistElem pedersen_zero_commit(const RistScal &x) {
    return scalar_mult_base(x);
}

// compute the pedersen commitment x * g, and write the result into y
inline void pedersen_zero_commit(RistElem &y, const RistScal &x) {
    scalar_mult_base(y, x);
}

// for each i, compute the pedersen commitment xx[i] * g + r * hh[i], where g is the fixed group generator, and write the result into yy[i]
inline void pedersen_commit(RistElemVec &yy, const RistScalVec &xx, const RistElemVec &hh, const RistScal &r) {
    std::for_each(
            std::execution::par_unseq,
            yy.begin(),
            yy.end(),
            [&yy, &xx, &hh, &r](auto &&y) {
                auto i = &y - yy.data();
                pedersen_commit(y, xx[i], hh[i], r);
            }
    );
}

// for each i, compute the pedersen commitment xx[i] * g, where g is the fixed group generator, and write the result into yy[i]
inline void pedersen_zero_commit(RistElemVec &yy, const RistScalVec &xx) {
    std::for_each(
            std::execution::par_unseq,
            yy.begin(),
            yy.end(),
            [&yy, &xx](auto &&y) {
                auto i = &y - yy.data();
                pedersen_zero_commit(y, xx[i]);
            }
    );
}

// for each i, compute the pedersen commitment xx[i] * g + x * hh[i], where g is the fixed group generator
inline RistElemVec pedersen_commit(const RistScalVec &xx, const RistElemVec &hh, const RistScal &r) {
    RistElemVec yy(xx.size());
    pedersen_commit(yy, xx, hh, r);
    return yy;
}

// a PedersenWithZeroProof instance is a proof that the prover has knowledge of y0 = r * g and y1 = x * g + r * h, for some x and r.
class PedersenWithZeroProof {
public:
    RistScal c;
    RistScal s0;
    RistScal s1;

    PedersenWithZeroProof() = default;

    void export_to_bytestream(std::vector<unsigned char>::iterator &it) const {
        c.export_to_bytestream(it);
        s0.export_to_bytestream(it);
        s1.export_to_bytestream(it);
    }

    void import_from_bytestream(std::vector<unsigned char>::const_iterator &it) {
        c.import_from_bytestream(it);
        s0.import_from_bytestream(it);
        s1.import_from_bytestream(it);
    }

    static int size() {
        return RISTSCALBYTES * 3;
    }
};

// produce a proof that the prover has knowledge of y0 = r * g and y1 = x * g + r * h, for some x and r.
PedersenWithZeroProof
generate_pedersen_with_zero_proof(const RistP3AndBytes &h, const RistP3AndBytes &y0, const RistP3AndBytes &y1,
                                  const RistScal &r, const RistScal &x);

// verifies the proof that the prover has knowledge of y0 = r * g and y1 = x * g + r * h, for some x and r.
bool b_verify_pedersen_with_zero(const RistP3AndBytes &h, const RistP3AndBytes &y0, const RistP3AndBytes &y1,
                                 const PedersenWithZeroProof &proof);


// a PedersenBatchWithZeroProof instance is a proof that the prover has knowledge of zz[j] = rr[j] * g,  yy[i] = xx[i] * g + sum_j (rr[j] * hh[i][j]), and yy_prime[i] = xx[i] * g +  qq[i] * f, for some rr, xx, qq
class PedersenBatchWithZeroProof {
public:
    RistP3VecAndBytes uu;
    RistP3VecAndBytes tt;
    RistP3VecAndBytes tt_star;

    RistScalVec ss;
    RistScalVec ss_star;
    RistScalVec ss_prime;

    PedersenBatchWithZeroProof() = default;

    PedersenBatchWithZeroProof(int num_blinds_per_group_element,
                               int num_samples_plus_one) :
            uu(num_blinds_per_group_element),
            tt(num_samples_plus_one),
            tt_star(num_samples_plus_one),
            ss(num_blinds_per_group_element),
            ss_star(num_samples_plus_one),
            ss_prime(num_samples_plus_one) {}

    void export_to_bytestream(std::vector<unsigned char>::iterator &it) const {
        uu.export_to_bytestream(it);
        tt.export_to_bytestream(it);
        tt_star.export_to_bytestream(it);
        ::export_to_bytestream(ss, it);
        ::export_to_bytestream(ss_star, it);
        ::export_to_bytestream(ss_prime, it);
    }

    void import_from_bytestream(std::vector<unsigned char>::const_iterator &it) {
        uu.import_from_bytestream(it);
        tt.import_from_bytestream(it);
        tt_star.import_from_bytestream(it);
        ::import_from_bytestream(ss, it);
        ::import_from_bytestream(ss_star, it);
        ::import_from_bytestream(ss_prime, it);
    }

    static int size(int num_blinds_per_group_element,
                    int num_samples_plus_one) {
        return (RISTBYTES + RISTSCALBYTES) *
               (num_blinds_per_group_element + num_samples_plus_one + num_samples_plus_one);
    }
};


// produce a proof that the prove has knowledge of zz[j] = rr[j] * g,  yy[i] = xx[i] * g + sum_j (rr[j] * hh[i][j]), and yy_prime[i] = xx[i] * g +  qq[i] * f, for some rr, xx, qq
PedersenBatchWithZeroProof
generate_pedersen_batch_with_zero_proof(const RistP3MatAndBytes &hh, const RistP3VecAndBytes &zz,
                                        const RistP3VecAndBytes &yy,
                                        const RistP3AndBytes &f,
                                        const RistP3VecAndBytes &yy_prime,
                                        const RistScalVec &rr, const RistScalVec &xx, const RistScalVec &qq);

// verify the proof that the prover has knowledge of zz[j] = rr[j] * g,  yy[i] = xx[i] * g + sum_j (rr[j] * hh[i][j]), and yy_prime[i] = xx[i] * g +  qq[i] * f, for some rr, xx, qq
bool
b_verify_pedersen_batch_with_zero(const RistP3MatAndBytes &hh, const RistP3VecAndBytes &zz, const RistP3VecAndBytes &yy,
                                  const RistP3AndBytes &f,
                                  const RistP3VecAndBytes &yy_prime,
                                  const PedersenBatchWithZeroProof &proof);


// a PedersenWithSquareProof instance is a proof that the prover has knowledge of y1 = x * g + r1 * h, y2 = x * x * g + r2 * h, for some x, r1, r2
struct PedersenWithSquareProof {
    RistScal c;
    RistScal s1;
    RistScal s2;
    RistScal s3;

    void export_to_bytestream(std::vector<unsigned char>::iterator &it) const {
        c.export_to_bytestream(it);
        s1.export_to_bytestream(it);
        s2.export_to_bytestream(it);
        s3.export_to_bytestream(it);
    }

    void import_from_bytestream(std::vector<unsigned char>::const_iterator &it) {
        c.import_from_bytestream(it);
        s1.import_from_bytestream(it);
        s2.import_from_bytestream(it);
        s3.import_from_bytestream(it);
    }

    static int size() {
        return RISTSCALBYTES * 4;
    }
};


// produce a proof that the prove has knowledge of y1 = x * g + r1 * h, y2 = x * x * g + r2 * h, for some x, r1, r2
PedersenWithSquareProof
generate_pedersen_with_square_proof(const RistP3AndBytes &h, const RistP3AndBytes &y1, const RistP3AndBytes &y2,
                                    const RistScal &x, const RistScal &r1, const RistScal &r2);

// verify the proof that the prove has knowledge of y1 = x * g + r1 * h, y2 = x * x * g + r2 * h, for some x, r1, r2
bool b_verify_pedersen_with_square(const RistP3AndBytes &h, const RistP3AndBytes &y1, const RistP3AndBytes &y2,
                                   const PedersenWithSquareProof &proof);


// a PedersenWithSquareProof instance is a proof that the prover has knowledge of yy1[i] = xx[i] * g + rr1[i] * h, yy2[i] = xx[i] * xx[i] * g + rr2[i] * h, for some xx, rr1, rr2
struct BatchPedersenWithSquareProof {
    RistP3VecAndBytes tt1;
    RistP3VecAndBytes tt2;
    RistScalVec ss1;
    RistScalVec ss2;
    RistScalVec ss3;

    void export_to_bytestream(std::vector<unsigned char>::iterator &it) const {
        tt1.export_to_bytestream(it);
        tt2.export_to_bytestream(it);
        ::export_to_bytestream(ss1, it);
        ::export_to_bytestream(ss2, it);
        ::export_to_bytestream(ss3, it);
    }

    void import_from_bytestream(std::vector<unsigned char>::const_iterator &it) {
        tt1.import_from_bytestream(it);
        tt2.import_from_bytestream(it);
        ::import_from_bytestream(ss1, it);
        ::import_from_bytestream(ss2, it);
        ::import_from_bytestream(ss3, it);
    }

    BatchPedersenWithSquareProof() = default;

    BatchPedersenWithSquareProof(int n) :
            tt1(n), tt2(n), ss1(n), ss2(n), ss3(n) {}

    static int size(int n) {
        return n * (RISTBYTES * 2 + RISTSCALBYTES * 3);
    }
};

// produce a proof that the prove has knowledge of yy1[i] = xx[i] * g + rr1[i] * h, yy2[i] = xx[i] * xx[i] * g + rr2[i] * h, for some xx, rr1, rr2
BatchPedersenWithSquareProof
generate_batch_pedersen_with_square_proof(const RistP3AndBytes &h, const RistP3VecAndBytes &yy1,
                                          const RistP3VecAndBytes &yy2,
                                          const RistScalVec &xx, const RistScalVec &rr1, const RistScalVec &rr2);

// verify the proof that the prove has knowledge of yy1[i] = xx[i] * g + rr1[i] * h, yy2[i] = xx[i] * xx[i] * g + rr2[i] * h, for some xx, rr1, rr2
bool
b_verify_batch_pedersen_with_square(const RistP3AndBytes &h, const RistP3VecAndBytes &yy1, const RistP3VecAndBytes &yy2,
                                    const BatchPedersenWithSquareProof &proof);



/* ************************/
/*    range proof tools   */
/* ************************/

// a bulletproofs proof with length bound_bits inner product protocol (not length log(bound_bits))
//  as length is not the primary concern
//  not compressing bound_bits to log(bound_bits) saves both prover time and verifier measure_time
// see "Bulletproofs: Short Proofs for Confidential Transactions and More"
// https://eprint.iacr.org/2017/1066.pdf


class RangeProofPowerTwoP3 {
public:
    RistP3AndBytes A, S;
    RistP3AndBytes T1, T2;
    RistScal tau_x, mu;
    RistScal t_hat;
    RistScalVec ll, rr;

    RangeProofPowerTwoP3() = default;

    RangeProofPowerTwoP3(int bound_bits) : ll(bound_bits), rr(bound_bits) {}

    void set_bound_bits(int bound_bits) {
        ll.resize(bound_bits);
        rr.resize(bound_bits);
    }

    void export_to_bytestream(std::vector<unsigned char>::iterator &it) const {
        A.export_to_bytestream(it);
        S.export_to_bytestream(it);
        T1.export_to_bytestream(it);
        T2.export_to_bytestream(it);
        tau_x.export_to_bytestream(it);
        mu.export_to_bytestream(it);
        t_hat.export_to_bytestream(it);
        ::export_to_bytestream(ll, it);
        ::export_to_bytestream(rr, it);
    }

    void import_from_bytestream(std::vector<unsigned char>::const_iterator &it) {
        A.import_from_bytestream(it);
        S.import_from_bytestream(it);
        T1.import_from_bytestream(it);
        T2.import_from_bytestream(it);
        tau_x.import_from_bytestream(it);
        mu.import_from_bytestream(it);
        t_hat.import_from_bytestream(it);
        ::import_from_bytestream(ll, it);
        ::import_from_bytestream(rr, it);

    }

    static int size(int bound_bits) {
        return RISTBYTES * 4 + RISTSCALBYTES * 3 + RISTSCALBYTES * bound_bits * 2;
    }
};

RangeProofPowerTwoP3 generate_range_proof_power_two(const RistP3AndBytes &h, const RistP3AndBytes &V,
                                                    const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                                    int n,
                                                    const RistScal &gamma, const RistScal &v);

bool b_verify_range_power_two(const RistP3AndBytes &h, const RistP3AndBytes &V,
                              const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                              int n,
                              const RangeProofPowerTwoP3 &proof);

RangeProofPowerTwoP3 generate_abs_val_range_proof_power_two(const RistP3AndBytes &h, const RistP3AndBytes &V,
                                                            const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                                            int n,
                                                            const RistScal &gamma, const RistScal &v);

bool b_verify_abs_val_range_power_two(const RistP3AndBytes &h, const RistP3AndBytes &V,
                                      const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                      int n,
                                      const RangeProofPowerTwoP3 &proof);


struct RangeProofP3 {
    RangeProofPowerTwoP3 lower;
    RangeProofPowerTwoP3 upper;

    RangeProofP3(int bound_bits) : lower(bound_bits), upper(bound_bits) {}

    void set_bound_bits(int bound_bits) {
        lower.set_bound_bits(bound_bits);
        upper.set_bound_bits(bound_bits);
    }

    void export_to_bytestream(std::vector<unsigned char>::iterator &it) const {
        lower.export_to_bytestream(it);
        upper.export_to_bytestream(it);
    }

    void import_from_bytestream(std::vector<unsigned char>::const_iterator &it) {
        lower.import_from_bytestream(it);
        upper.import_from_bytestream(it);
    }

    static int size(int bound_bits) {
        return 2 * RangeProofPowerTwoP3::size(bound_bits);
    }
};

// produces a proof that a commitment of v, V = v * g + gamma * h, satisfies that v is in the interval [0, bound)
RangeProofP3 generate_range_proof(const RistP3AndBytes &h, const RistP3AndBytes &V,
                                  const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                  const RistScal &bound,
                                  const RistScal &gamma, const RistScal &v);

bool b_verify_range(const RistP3AndBytes &h, const RistP3AndBytes &V,
                    const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                    const RistScal &bound,
                    const RangeProofP3 &proof);

RangeProofPowerTwoP3 generate_range_proof_power_two_agg(const RistP3AndBytes &h, const RistP3VecAndBytes &VV,
                                                        const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                                        int n, int m, // m values of range [0, 2^n)
                                                        const RistScalVec &gamma, const RistScalVec &vv);

bool b_verify_range_power_two_agg(const RistP3AndBytes &h, const RistP3VecAndBytes &VV,
                                  const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                  int n, int m, // m values of range [0, 2^n)
                                  const RangeProofPowerTwoP3 &proof);

RangeProofPowerTwoP3 generate_abs_val_range_proof_power_two_agg(const RistP3AndBytes &h, const RistP3VecAndBytes &VV,
                                                                const RistP3VecAndBytes &gg,
                                                                const RistP3VecAndBytes &hh,
                                                                int n, int m, // m values of range [-2^(n-1), 2^(n-1))
                                                                const RistScalVec &gamma, const RistScalVec &vv);

bool b_verify_abs_val_range_power_two_agg(const RistP3AndBytes &h, const RistP3VecAndBytes &VV,
                                          const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                          int n, int m, // m values of range [-2^(n-1), 2^(n-1))
                                          const RangeProofPowerTwoP3 &proof);

class InnerProdProof {
public:
    RistP3VecAndBytes LL, RR;
    RistScal a, b;

    InnerProdProof() = default;

    InnerProdProof(int n) :
            LL(get_bit_length(n - 1)),
            RR(get_bit_length(n - 1)) {}

    void export_to_bytestream(std::vector<unsigned char>::iterator &it) const {
        LL.export_to_bytestream(it);
        RR.export_to_bytestream(it);
        a.export_to_bytestream(it);
        b.export_to_bytestream(it);
    }

    void import_from_bytestream(std::vector<unsigned char>::const_iterator &it) {
        LL.import_from_bytestream(it);
        RR.import_from_bytestream(it);
        a.import_from_bytestream(it);
        b.import_from_bytestream(it);
    }

    static int size(int n) {
        return RISTBYTES * 2 * get_bit_length(n - 1) + RISTSCALBYTES * 2;
    }
};

InnerProdProof generate_inner_prod_proof(const RistP3VecAndBytes &gg,
                                         const RistP3VecAndBytes &hh,
                                         const RistP3AndBytes &u,
                                         const RistP3AndBytes &P,
                                         const RistScalVec &aa, const RistScalVec &bb);

bool b_verify_inner_prod_proof(const RistP3VecAndBytes &gg,
                               const RistP3VecAndBytes &hh,
                               const RistP3AndBytes &u,
                               const RistP3AndBytes &P,
                               const InnerProdProof &proof);

class RangeProofPowerTwoP3Log {
public:
    RistP3AndBytes A, S;
    RistP3AndBytes T1, T2;
    RistScal tau_x, mu;
    RistScal t_hat;
    InnerProdProof inner_prod_proof;

    RangeProofPowerTwoP3Log() = default;

    RangeProofPowerTwoP3Log(int n) :
            inner_prod_proof(n) {}

    void export_to_bytestream(std::vector<unsigned char>::iterator &it) const {
        A.export_to_bytestream(it);
        S.export_to_bytestream(it);
        T1.export_to_bytestream(it);
        T2.export_to_bytestream(it);
        tau_x.export_to_bytestream(it);
        mu.export_to_bytestream(it);
        t_hat.export_to_bytestream(it);
        inner_prod_proof.export_to_bytestream(it);
    }

    void import_from_bytestream(std::vector<unsigned char>::const_iterator &it) {
        A.import_from_bytestream(it);
        S.import_from_bytestream(it);
        T1.import_from_bytestream(it);
        T2.import_from_bytestream(it);
        tau_x.import_from_bytestream(it);
        mu.import_from_bytestream(it);
        t_hat.import_from_bytestream(it);
        inner_prod_proof.import_from_bytestream(it);
    }

    static int size(int bound_bits) {
        return RISTBYTES * 4 + RISTSCALBYTES * 3 + InnerProdProof::size(bound_bits);
    }
};

RangeProofPowerTwoP3Log generate_range_proof_power_two_log(const RistP3AndBytes &h, const RistP3AndBytes &V,
                                                           const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                                           int n,
                                                           const RistScal &gamma, const RistScal &v);

bool b_verify_range_power_two_log(const RistP3AndBytes &h, const RistP3AndBytes &V,
                                  const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                  int n,
                                  const RangeProofPowerTwoP3Log &proof);


RangeProofPowerTwoP3Log generate_range_proof_power_two_agg_log(const RistP3AndBytes &h, const RistP3VecAndBytes &VV,
                                                               const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                                               int n, int m, // m values of range [0, 2^n)
                                                               const RistScalVec &gamma, const RistScalVec &vv);


bool b_verify_range_power_two_agg_log(const RistP3AndBytes &h, const RistP3VecAndBytes &VV,
                                      const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                      int n, int m, // m values of range [0, 2^n)
                                      const RangeProofPowerTwoP3Log &proof);


RangeProofPowerTwoP3Log
generate_abs_val_range_proof_power_two_agg_log(const RistP3AndBytes &h, const RistP3VecAndBytes &VV,
                                               const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                               int n, int m, // m values of range [-2^(n-1), 2^(n-1))
                                               const RistScalVec &gamma, const RistScalVec &vv);


bool b_verify_abs_val_range_power_two_agg_log(const RistP3AndBytes &h, const RistP3VecAndBytes &VV,
                                              const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                              int n, int m, // m values of range [-2^(n-1), 2^(n-1))
                                              const RangeProofPowerTwoP3Log &proof);

#endif //RISEFL_CRYPTO_ZKP_H
