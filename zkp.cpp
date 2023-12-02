//
// Created by yizheng on 8/3/23.
//

#include "include/zkp.h"
#include "include/zkp_hash.h"
#include "include/ristretto.h"
#include "include/ristretto_vector.h"
#include "include/rist_fast_computation.h"

PedersenWithZeroProof
generate_pedersen_with_zero_proof(const RistP3AndBytes &h,
                                  const RistP3AndBytes &y0, const RistP3AndBytes &y1,
                                  const RistScal &r, const RistScal &x) {
    RistScal v0, v1;
    rand_init(v0);
    rand_init(v1);

    RistP3AndBytes t0;
    t0.elem = pedersen_zero_commit_p3(v0);
    t0.fill_bytes();

    RistP3AndBytes t1;
    pedersen_commit(t1.elem, v1, h.elem, v0);
    t1.fill_bytes();

    PedersenWithZeroProof proof;
    proof.c = rist_scalar_from_hash_from_bytes(h.bytes, y0.bytes, y1.bytes, t0.bytes, t1.bytes);
    proof.s0 = v0 - proof.c * r;
    proof.s1 = v1 - proof.c * x;

    return proof;
}

bool b_verify_pedersen_with_zero(const RistP3AndBytes &h,
                                 const RistP3AndBytes &y0, const RistP3AndBytes &y1,
                                 const PedersenWithZeroProof &proof) {
    RistP3AndBytes t0prime;
    t0prime.elem = pedersen_commit(proof.s0, y0.elem, proof.c);
    t0prime.fill_bytes();

    RistP3AndBytes t1prime;
    pedersen_commit(t1prime.elem, proof.s1, h.elem, proof.s0);
    t1prime = t1prime + proof.c * y1.elem;
    t1prime.fill_bytes();

    return proof.c ==
           rist_scalar_from_hash_from_bytes(h.bytes, y0.bytes, y1.bytes, t0prime.bytes, t1prime.bytes);
}


PedersenBatchWithZeroProof
generate_pedersen_batch_with_zero_proof(const RistP3MatAndBytes &hh, const RistP3VecAndBytes &zz,
                                        const RistP3VecAndBytes &yy,
                                        const RistP3AndBytes &f,
                                        const RistP3VecAndBytes &yy_prime,
                                        const RistScalVec &rr, const RistScalVec &xx, const RistScalVec &qq) {
    int num_blinds_per_group_element = zz.size();
    assert(num_blinds_per_group_element == hh.size_cols());
    assert(num_blinds_per_group_element == rr.size());
    int num_samples_plus_one = yy.size();
    assert(num_samples_plus_one == hh.size_rows());
    assert(num_samples_plus_one == xx.size());
    assert(num_samples_plus_one == yy_prime.size());
    assert(num_samples_plus_one == qq.size());

    RistScalVec ww(num_blinds_per_group_element);
    RistScalVec vv(num_samples_plus_one);
    RistScalVec vv_prime(num_samples_plus_one);
    rand_init(ww);
    rand_init(vv);
    rand_init(vv_prime);

    PedersenBatchWithZeroProof proof(num_blinds_per_group_element, num_samples_plus_one);

//    RistP3VecAndBytes uu(num_blinds_per_group_element);
    pedersen_zero_commit(proof.uu.elems, ww);
    proof.uu.fill_bytes();

//    RistP3VecAndBytes tt(num_samples_plus_one);
    for (int i = 0; i < proof.tt.size(); i++) {
        pedersen_zero_commit(proof.tt.elems[i], vv[i]);
        proof.tt.elems[i] += linear_comb(ww, hh.elems[i]);
    }
    proof.tt.fill_bytes();

//    RistP3VecAndBytes tt_star(num_samples_plus_one);
    pedersen_commit(proof.tt_star.elems, vv, f.elem, vv_prime);
    proof.tt_star.fill_bytes();

    auto c = rist_scalar_from_hash_from_bytes(hh.bytes, f.bytes, zz.bytes, yy.bytes, yy_prime.bytes,
                                              proof.uu.bytes, proof.tt.bytes, proof.tt_star.bytes);
    proof.ss = ww - c * rr;
    proof.ss_star = vv - c * xx;
    proof.ss_prime = vv_prime - c * qq;

    return proof;
}

bool
b_verify_pedersen_batch_with_zero(const RistP3MatAndBytes &hh, const RistP3VecAndBytes &zz, const RistP3VecAndBytes &yy,
                                  const RistP3AndBytes &f,
                                  const RistP3VecAndBytes &yy_prime,
                                  const PedersenBatchWithZeroProof &proof) {
    int num_blinds_per_group_element = proof.uu.size();
    int num_samples_plus_one = proof.tt.size();
    assert(num_samples_plus_one == proof.tt_star.size());

    RistScalVec alpha(num_blinds_per_group_element), beta(num_samples_plus_one), gamma(num_samples_plus_one);
    rand_init(alpha);
    rand_init(beta);
    rand_init(gamma);

    auto c = rist_scalar_from_hash_from_bytes(hh.bytes, f.bytes, zz.bytes, yy.bytes, yy_prime.bytes,
                                              proof.uu.bytes, proof.tt.bytes, proof.tt_star.bytes);

    LinearCombCalculator lcc;
    lcc.add(alpha, proof.uu.elems);
    lcc.add(beta, proof.tt.elems);
    lcc.add(gamma, proof.tt_star.elems);

    lcc.sub_base(inner_prod(alpha, proof.ss) + inner_prod(beta + gamma, proof.ss_star));
    lcc.sub(inner_prod(gamma, proof.ss_prime), f.elem);
    lcc.sub(c * alpha, zz.elems);
    lcc.sub(c * beta, yy.elems);
    lcc.sub(c * gamma, yy_prime.elems);
    RistScalMat hh_coeffs(hh.size_rows());
    for (int i = 0; i < hh_coeffs.size(); i++) {
        hh_coeffs[i] = beta[i] * proof.ss;
    }
    lcc.sub(hh_coeffs, hh.elems);

    return lcc.is_0();

//    auto left = linear_comb(alpha, proof.uu.elems) + linear_comb(beta, proof.tt.elems) +
//                            linear_comb(gamma, proof.tt_star.elems);
//    auto right = pedersen_commit(inner_prod(alpha, proof.ss) + inner_prod(beta + gamma, proof.ss_star),
//                                 f.elem,
//                                 inner_prod(gamma, proof.ss_prime))
//                 +
//                 linear_comb(c * alpha, zz.elems)
//                 +
//                linear_comb(c * beta, yy.elems)
//                +
//            linear_comb(c * gamma, yy_prime.elems);
//
//    right += linear_comb(hh_coeffs, hh.elems);
//
//    return left == right;
//
//    RistP3VecAndBytes ttprime(yy.size());
//    for (int i = 0; i < ttprime.size(); i++) {
//        pedersen_commit(ttprime.elems[i], proof.ss_star[i], yy.elems[i], proof.c);
//        ttprime.elems[i] += linear_comb(proof.ss, hh.elems[i]);
//    }
//    ttprime.fill_bytes();
//
//    RistP3VecAndBytes uuprime(zz.size());
//    pedersen_commit(uuprime.elems, proof.ss, zz.elems, proof.c);
//    uuprime.fill_bytes();
//
//    RistP3VecAndBytes tt_star_prime(yy.size());
//    for (int i = 0; i < tt_star_prime.size(); i++) {
//        pedersen_commit(tt_star_prime.elems[i], proof.ss_star[i], yy_prime.elems[i], proof.c);
////        assert(tt_star_prime.elems[i] == pedersen_zero_commit_p3(proof.ss_star[i]) + proof.c * yy_prime.elems[i]);
//        tt_star_prime.elems[i] += proof.ss_prime[i] * f.elem;
//    }
//    tt_star_prime.fill_bytes();
//
//    return proof.c ==
//           rist_scalar_from_hash_from_bytes(hh.bytes, f.bytes, zz.bytes, yy.bytes, yy_prime.bytes,
//                                            uuprime.bytes, ttprime.bytes, tt_star_prime.bytes);
}

PedersenWithSquareProof
generate_pedersen_with_square_proof(const RistP3AndBytes &h, const RistP3AndBytes &y1, const RistP3AndBytes &y2,
                                    const RistScal &x, const RistScal &r1, const RistScal &r2) {
    RistScal v1, v2, v3;
    rand_init(v1);
    rand_init(v2);
    rand_init(v3);

    RistP3AndBytes t1;
    t1.elem = pedersen_commit(v1, h.elem, v2);
    t1.fill_bytes();

    RistP3AndBytes t2;
    t2.elem = v1 * y1.elem + v3 * h.elem;
    t2.fill_bytes();

    PedersenWithSquareProof proof;
    proof.c = rist_scalar_from_hash_from_bytes(h.bytes, y1.bytes, y2.bytes,
                                               t1.bytes, t2.bytes);
    proof.s1 = v1 - proof.c * x;
    proof.s2 = v2 - proof.c * r1;
    proof.s3 = v3 - proof.c * (r2 - r1 * x);
    return proof;
}

bool b_verify_pedersen_with_square(const RistP3AndBytes &h, const RistP3AndBytes &y1, const RistP3AndBytes &y2,
                                   const PedersenWithSquareProof &proof) {
    RistP3AndBytes t1prime;
    t1prime.elem = pedersen_commit(proof.s1, h.elem, proof.s2) + proof.c * y1.elem;
    t1prime.fill_bytes();

    RistP3AndBytes t2prime;
    t2prime.elem = proof.s1 * y1.elem + proof.s3 * h.elem + proof.c * y2.elem;
    t2prime.fill_bytes();

    return proof.c ==
           rist_scalar_from_hash_from_bytes(h.bytes, y1.bytes, y2.bytes,
                                            t1prime.bytes, t2prime.bytes);
}


BatchPedersenWithSquareProof
generate_batch_pedersen_with_square_proof(const RistP3AndBytes &h, const RistP3VecAndBytes &yy1,
                                          const RistP3VecAndBytes &yy2,
                                          const RistScalVec &xx, const RistScalVec &rr1, const RistScalVec &rr2) {
    int n = yy1.size();
    assert(n == yy2.size());
    assert(n == xx.size());
    assert(n == rr1.size());
    assert(n == rr2.size());

    RistScalVec vv1(n), vv2(n), vv3(n);
    rand_init(vv1);
    rand_init(vv2);
    rand_init(vv3);

    BatchPedersenWithSquareProof proof(n);

    pedersen_commit(proof.tt1.elems, vv1, h.elem, vv2);
    proof.tt1.fill_bytes();

    proof.tt2.elems = vv1 * yy1.elems + vv3 * RistElemP3Vec(n, h.elem);
    proof.tt2.fill_bytes();

    auto c = rist_scalar_from_hash_from_bytes(h.bytes, yy1.bytes, yy2.bytes,
                                              proof.tt1.bytes, proof.tt2.bytes);
    proof.ss1 = vv1 - c * xx;
    proof.ss2 = vv2 - c * rr1;
    proof.ss3 = vv3 - c * (rr2 - rr1 * xx);
    return proof;
}

bool
b_verify_batch_pedersen_with_square(const RistP3AndBytes &h, const RistP3VecAndBytes &yy1, const RistP3VecAndBytes &yy2,
                                    const BatchPedersenWithSquareProof &proof) {
    int n = yy1.size();
    assert(n == yy2.size());
    assert(n == proof.tt1.size());

    RistScalVec alpha(n), beta(n);
    rand_init(alpha);
    rand_init(beta);
    auto c = rist_scalar_from_hash_from_bytes(h.bytes, yy1.bytes, yy2.bytes,
                                              proof.tt1.bytes, proof.tt2.bytes);

    LinearCombCalculator lcc;
    lcc.add(alpha, proof.tt1.elems);
    lcc.add(beta, proof.tt2.elems);
    lcc.sub_base(inner_prod(alpha, proof.ss1));
    lcc.sub(inner_prod(alpha, proof.ss2) + inner_prod(beta, proof.ss3), h.elem);
    lcc.sub(c * alpha + beta * proof.ss1, yy1.elems);
    lcc.sub(c * beta, yy2.elems);

    return lcc.is_0();

//    auto left = linear_comb(alpha, proof.tt1.elems) + linear_comb(beta, proof.tt2.elems);
//
//    auto right = pedersen_commit(inner_prod(alpha, proof.ss1),
//                                 h.elem,
//                                 inner_prod(alpha, proof.ss2) + inner_prod(beta, proof.ss3))
//                 +
//                 linear_comb(c * alpha + beta * proof.ss1, yy1.elems)
//                 +
//                 linear_comb(c * beta, yy2.elems);
//    return left == right;
}

RistScalVec scalar_geometric_series_power_two(int n) {
    RistScalVec series(n);
    for (int i = 0; i < n; i++) {
        series[i] = power_of_two(i);
    }
    return series;
}

RangeProofPowerTwoP3 generate_range_proof_power_two(const RistP3AndBytes &h, const RistP3AndBytes &V,
                                                    const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                                    int n,
                                                    const RistScal &gamma, const RistScal &v) {
    RangeProofPowerTwoP3 proof(n);
    RistScalVec aa_L(n), aa_R(n);
    RistP3VecAndBytes gg_segment = gg.segment(n);
    RistP3VecAndBytes hh_segment = hh.segment(n);
    for (int b = 0; b < n; b++) {
        aa_L[b] = b_get_bit(v, b) ? c_scal_one : c_scal_zero;
        aa_R[b] = aa_L[b] - c_scal_one;
    }

    RistScal alpha;
    rand_init(alpha);
    LinearCombCalculator lcc_A;
    lcc_A.add(alpha, h.elem);
    lcc_A.add(aa_L, gg_segment.elems);
    lcc_A.add(aa_R, hh_segment.elems);
    proof.A.elem = lcc_A.result();
//    proof.A.elem = alpha * h.elem + linear_comb(aa_L, gg_segment.elems) + linear_comb(aa_R, hh_segment.elems);
    proof.A.fill_bytes();

    RistScalVec ss_L(n), ss_R(n);
    RistScal rho;
    rand_init(ss_L);
    rand_init(ss_R);
    rand_init(rho);
    LinearCombCalculator lcc_S;
    lcc_S.add(rho, h.elem);
    lcc_S.add(ss_L, gg_segment.elems);
    lcc_S.add(ss_R, hh_segment.elems);
    proof.S.elem = lcc_S.result();
//    proof.S.elem = rho * h.elem + linear_comb(ss_L, gg_segment.elems) + linear_comb(ss_R, hh_segment.elems);
    proof.S.fill_bytes();

    RistScal y = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  V.bytes,
                                                  proof.A.bytes, proof.S.bytes);


    RistScal z = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  V.bytes,
                                                  proof.A.bytes, proof.S.bytes,
                                                  y.scalar);

    RistScal tau1, tau2;
    rand_init(tau1);
    rand_init(tau2);
    RistScalVec ll_0 = aa_L - RistScalVec(n, z);
    RistScalVec ll_1 = ss_L;
    RistScalVec rr_0 = scalar_geometric_series(n, y) * (aa_R + RistScalVec(n, z))
                       + z * z * scalar_geometric_series_power_two(n);
    RistScalVec rr_1 = scalar_geometric_series(n, y) * ss_R;
    RistScal t1 = inner_prod(ll_0, rr_1) + inner_prod(ll_1, rr_0);
    RistScal t2 = inner_prod(ll_1, rr_1);
    proof.T1.elem = pedersen_commit(t1, h.elem, tau1);
    proof.T1.fill_bytes();
    proof.T2.elem = pedersen_commit(t2, h.elem, tau2);
    proof.T2.fill_bytes();

    RistScal x = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  V.bytes,
                                                  proof.A.bytes, proof.S.bytes,
                                                  y.scalar, z.scalar,
                                                  proof.T1.bytes, proof.T2.bytes);


    proof.ll = aa_L - RistScalVec(n, z) + x * ss_L;
    proof.rr = scalar_geometric_series(n, y) * (aa_R + RistScalVec(n, z) + x * ss_R)
               + z * z * scalar_geometric_series_power_two(n);
    proof.t_hat = inner_prod(proof.ll, proof.rr);
    proof.tau_x = tau2 * x * x + tau1 * x + z * z * gamma;
    proof.mu = alpha + rho * x;

    return proof;
}

bool b_verify_range_power_two(const RistP3AndBytes &h, const RistP3AndBytes &V,
                              const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                              int n,
                              const RangeProofPowerTwoP3 &proof) {
    RistP3VecAndBytes gg_segment = gg.segment(n);
    RistP3VecAndBytes hh_segment = hh.segment(n);

    RistScal y = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  V.bytes,
                                                  proof.A.bytes, proof.S.bytes);

    RistScal z = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  V.bytes,
                                                  proof.A.bytes, proof.S.bytes,
                                                  y.scalar);


    RistScal delta = inner_prod(RistScalVec(n, z - z * z), scalar_geometric_series(n, y))
                     - inner_prod(RistScalVec(n, z * z * z), scalar_geometric_series_power_two(n));

    RistScal x = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  V.bytes,
                                                  proof.A.bytes, proof.S.bytes,
                                                  y.scalar, z.scalar,
                                                  proof.T1.bytes, proof.T2.bytes);


    if (pedersen_commit(proof.t_hat, h.elem, proof.tau_x) !=
        z * z * V.elem + pedersen_commit(delta, proof.T1.elem, x) + x * x * proof.T2.elem)
        return false;
    if (proof.t_hat != inner_prod(proof.ll, proof.rr))
        return false;

    LinearCombCalculator lcc;
    lcc.sub(proof.A.elem);
    lcc.sub(x, proof.S.elem);
    lcc.add(RistScalVec(n, z) + proof.ll, gg_segment.elems);
    lcc.add(scalar_geometric_series(n, scalar_invert(y)) *
            (proof.rr - z * z * scalar_geometric_series_power_two(n))
            - RistScalVec(n, z), hh_segment.elems);
    lcc.add(proof.mu, h.elem);
    return lcc.is_0();
//    return proof.A.elem + x * proof.S.elem ==
//           linear_comb(RistScalVec(n, z) + proof.ll, gg_segment.elems)
//           +
//           linear_comb(scalar_geometric_series(n, scalar_invert(y)) *
//                       (proof.rr - z * z * scalar_geometric_series_power_two(n))
//                       - RistScalVec(n, z), hh_segment.elems)
//           +
//           proof.mu * h.elem;
}

RangeProofPowerTwoP3 generate_abs_val_range_proof_power_two(const RistP3AndBytes &h, const RistP3AndBytes &V,
                                                            const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                                            int n,
                                                            const RistScal &gamma, const RistScal &v) {
    RistScal v_shifted = v + power_of_two(n);

    RistP3AndBytes V_shifted;
    V_shifted.elem = V.elem + pedersen_zero_commit_p3(power_of_two(n));
    V_shifted.fill_bytes();

    return generate_range_proof_power_two(h, V_shifted, gg, hh, n + 1, gamma, v_shifted);
}

bool b_verify_abs_val_range_power_two(const RistP3AndBytes &h, const RistP3AndBytes &V,
                                      const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                      int n,
                                      const RangeProofPowerTwoP3 &proof) {
    RistP3AndBytes V_shifted;
    V_shifted.elem = V.elem + pedersen_zero_commit_p3(power_of_two(n));
    V_shifted.fill_bytes();

    return b_verify_range_power_two(h, V_shifted, gg, hh, n + 1, proof);
}


RangeProofP3 generate_range_proof(const RistP3AndBytes &h, const RistP3AndBytes &V,
                                  const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                  const RistScal &bound,
                                  const RistScal &gamma, const RistScal &v) {

    int n = get_power_two_bound(bound);
    RangeProofP3 proof(n);
    proof.lower = generate_range_proof_power_two(h, V, gg, hh, n, gamma, v);
    RistScal diff = power_of_two(n) - bound;
    proof.upper = generate_range_proof_power_two(h, V + pedersen_zero_commit_p3(diff), gg, hh, n, gamma, v + diff);

    return proof;
}

bool b_verify_range(const RistP3AndBytes &h, const RistP3AndBytes &V,
                    const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                    const RistScal &bound,
                    const RangeProofP3 &proof) {
    int n = get_power_two_bound(bound);
    RistScal diff = power_of_two(n) - bound;
    return b_verify_range_power_two(h, V, gg, hh, n, proof.lower) &&
           b_verify_range_power_two(h, V + pedersen_zero_commit_p3(diff), gg, hh, n, proof.upper);
}

RangeProofPowerTwoP3 generate_range_proof_power_two_agg(const RistP3AndBytes &h, const RistP3VecAndBytes &VV,
                                                        const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                                        int n, int m, // m values of range [0, 2^n)
                                                        const RistScalVec &gamma, const RistScalVec &vv) {
    assert(m == VV.size());
    assert(m == vv.size());
    assert(m == gamma.size());

    assert(gg.size() >= n * m);
    assert(hh.size() >= n * m);

    RangeProofPowerTwoP3 proof(n * m);
    RistScalVec aa_L(n * m), aa_R(n * m);
    RistP3VecAndBytes gg_segment = gg.segment(n * m);
    RistP3VecAndBytes hh_segment = hh.segment(n * m);
    for (int j = 0; j < m; j++) {
        for (int b = 0; b < n; b++) {
            aa_L[j * n + b] = b_get_bit(vv[j], b) ? c_scal_one : c_scal_zero;
            aa_R[j * n + b] = aa_L[j * n + b] - c_scal_one;
        }
    }
    RistScal alpha;
    rand_init(alpha);
    LinearCombCalculator lcc_A;
    lcc_A.add(alpha, h.elem);
    lcc_A.add(aa_L, gg_segment.elems);
    lcc_A.add(aa_R, hh_segment.elems);
    proof.A.elem = lcc_A.result();
//    proof.A.elem = alpha * h.elem + linear_comb(aa_L, gg_segment.elems) + linear_comb(aa_R, hh_segment.elems);
    proof.A.fill_bytes();

    RistScalVec ss_L(n * m), ss_R(n * m);
    RistScal rho;
    rand_init(ss_L);
    rand_init(ss_R);
    rand_init(rho);
    LinearCombCalculator lcc_S;
    lcc_S.add(rho, h.elem);
    lcc_S.add(ss_L, gg_segment.elems);
    lcc_S.add(ss_R, hh_segment.elems);
    proof.S.elem = lcc_S.result();
//    proof.S.elem = rho * h.elem + linear_comb(ss_L, gg_segment.elems) + linear_comb(ss_R, hh_segment.elems);
    proof.S.fill_bytes();

    RistScal y = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  VV.bytes,
                                                  proof.A.bytes, proof.S.bytes);


    RistScal z = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  VV.bytes,
                                                  proof.A.bytes, proof.S.bytes,
                                                  y.scalar);

    RistScal tau1, tau2;
    rand_init(tau1);
    rand_init(tau2);
    RistScalVec ll_0 = aa_L - RistScalVec(n * m, z);
    RistScalVec ll_1 = ss_L;

    RistScalVec rr_0 = scalar_geometric_series(n * m, y) * (aa_R + RistScalVec(n * m, z));
    RistScalVec power_segments;
    RistScalVec segment = z * scalar_geometric_series_power_two(n);
    for (int j = 0; j < m; j++) {
        segment = z * segment;
        power_segments.insert(power_segments.end(), segment.begin(), segment.end());
    }
    assert(power_segments.size() == n * m);
    rr_0 = rr_0 + power_segments;

    RistScalVec rr_1 = scalar_geometric_series(n * m, y) * ss_R;
    RistScal t1 = inner_prod(ll_0, rr_1) + inner_prod(ll_1, rr_0);
    RistScal t2 = inner_prod(ll_1, rr_1);
    proof.T1.elem = pedersen_commit(t1, h.elem, tau1);
    proof.T1.fill_bytes();
    proof.T2.elem = pedersen_commit(t2, h.elem, tau2);
    proof.T2.fill_bytes();

    RistScal x = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  VV.bytes,
                                                  proof.A.bytes, proof.S.bytes,
                                                  y.scalar, z.scalar,
                                                  proof.T1.bytes, proof.T2.bytes);


    proof.ll = ll_0 + x * ll_1;
    proof.rr = rr_0 + x * rr_1;
    proof.t_hat = inner_prod(proof.ll, proof.rr);

    proof.tau_x = tau2 * x * x + tau1 * x;
    RistScal z_gamma_power_sum = c_scal_zero;
    RistScal power_of_z = z;
    for (int j = 0; j < m; j++) {
        power_of_z *= z;
        z_gamma_power_sum += power_of_z * gamma[j];
    }
    proof.tau_x += z_gamma_power_sum;

    proof.mu = alpha + rho * x;

//    auto y_inv = scalar_invert(y);
//    auto y_inv_series = scalar_geometric_series(n * m, y_inv);
//    RistP3VecAndBytes hh_prime(n * m);
//    hh_prime.elems = y_inv_series * hh_segment.elems;
//    hh_prime.fill_bytes();

//    RistP3AndBytes P;
//    P.elem = proof.mu * h.elem + linear_comb(ll, gg_segment.elems) + linear_comb(rr, hh_prime.elems);
//    P.fill_bytes();
//
//    RistP3AndBytes inner_prod_commitment;
//    inner_prod_commitment.elem = P.elem + pedersen_commit(proof.t_hat, h.elem, -proof.mu);
//    inner_prod_commitment.fill_bytes();
//
//    proof.inner_prod_proof = generate_inner_prod_proof(gg_segment, hh_prime, c_p3_bytes_one, inner_prod_commitment,
//                                                       ll, rr);

    return proof;
}


bool b_verify_range_power_two_agg(const RistP3AndBytes &h, const RistP3VecAndBytes &VV,
                                  const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                  int n, int m, // m values of range [0, 2^n)
                                  const RangeProofPowerTwoP3 &proof) {
    assert(gg.size() >= n * m);
    assert(hh.size() >= n * m);

    RistP3VecAndBytes gg_segment = gg.segment(n * m);
    RistP3VecAndBytes hh_segment = hh.segment(n * m);
    assert(VV.size() == m);

    RistScal y = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  VV.bytes,
                                                  proof.A.bytes, proof.S.bytes);

    RistScal z = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  VV.bytes,
                                                  proof.A.bytes, proof.S.bytes,
                                                  y.scalar);


    RistScal delta = (z - z * z) * (power(n * m, y) - c_scal_one) / (y - c_scal_one) -
                     (power(n, RistScal(2)) - c_scal_one) * z * z * z * (power(m, z) - c_scal_one) / (z - c_scal_one);

    RistScal x = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  VV.bytes,
                                                  proof.A.bytes, proof.S.bytes,
                                                  y.scalar, z.scalar,
                                                  proof.T1.bytes, proof.T2.bytes);

    if (pedersen_commit(proof.t_hat, h.elem, proof.tau_x) !=
        z * z * linear_comb(scalar_geometric_series(m, z), VV.elems)
        +
        pedersen_commit(delta, proof.T1.elem, x) + x * x * proof.T2.elem)
        return false;
    if (proof.t_hat != inner_prod(proof.ll, proof.rr))
        return false;
    RistScalVec power_segments;
    RistScalVec segment = z * scalar_geometric_series_power_two(n);
    for (int j = 0; j < m; j++) {
        segment = z * segment;
        power_segments.insert(power_segments.end(), segment.begin(), segment.end());
    }

    LinearCombCalculator lcc;
    lcc.sub(proof.A.elem);
    lcc.sub(x, proof.S.elem);
    lcc.add(RistScalVec(n * m, z) + proof.ll, gg_segment.elems);
    lcc.add(scalar_geometric_series(n * m, scalar_invert(y)) *
            (proof.rr - power_segments)
            - RistScalVec(n * m, z), hh_segment.elems);
    lcc.add(proof.mu, h.elem);
    return lcc.is_0();


//    return proof.A.elem + x * proof.S.elem ==
//            linear_comb(RistScalVec(n * m, z) + proof.ll, gg_segment.elems)
//            +
//            linear_comb(scalar_geometric_series(n * m, scalar_invert(y)) *
//                        (proof.rr - power_segments)
//                        - RistScalVec(n * m, z), hh_segment.elems)
//            +
//            proof.mu * h.elem;
//
//
//
//    RistElemP3 P;
//    RistScalVec power_segments;
//    RistScalVec segment = z * scalar_geometric_series_power_two(n);
//
//    for (int j = 0; j < m; j++) {
//        segment = z * segment;
//        power_segments.insert(power_segments.end(), segment.begin(), segment.end());
//    }
//    assert(power_segments.size() == n * m);
//    RistScalVec coeffs_temp = scalar_geometric_series(n * m, scalar_invert(y)) * power_segments + RistScalVec(n * m, z);
//
//    P = proof.A.elem + x * proof.S.elem - z * sum_single_thread(gg_segment.elems) +
//        linear_comb(coeffs_temp, hh_segment.elems);
//
//    RistElemP3 inner_prod_commitment;
//    inner_prod_commitment = P + pedersen_commit(proof.t_hat, h.elem, -proof.mu);
////    inner_prod_commitment.fill_bytes();
//
//    auto y_inv = scalar_invert(y);
//    auto y_inv_series = scalar_geometric_series(n * m, y_inv);
////    RistP3VecAndBytes hh_prime(n * m);
////    hh_prime.elems = y_inv_series * hh_segment.elems;
////    hh_prime.fill_bytes();
//
//    return inner_prod_commitment == linear_comb(proof.ll, gg_segment.elems) + linear_comb(proof.rr * y_inv_series, hh_segment.elems) +
//                                    proof.t_hat * c_p3_bytes_one.elem;
//
////    return b_verify_inner_prod_proof(gg_segment, hh_prime, c_p3_bytes_one, inner_prod_commitment,
////                                     proof.inner_prod_proof);
}

RangeProofPowerTwoP3 generate_abs_val_range_proof_power_two_agg(const RistP3AndBytes &h, const RistP3VecAndBytes &VV,
                                                                const RistP3VecAndBytes &gg,
                                                                const RistP3VecAndBytes &hh,
                                                                int n, int m, // m values of range [-2^(n-1), 2^(n-1))
                                                                const RistScalVec &gamma, const RistScalVec &vv) {
    RistScalVec vv_shifted = vv + RistScalVec(m, power_of_two(n - 1));

    RistP3VecAndBytes VV_shifted(m);
    VV_shifted.elems = VV.elems + RistElemP3Vec(m, pedersen_zero_commit_p3(power_of_two(n - 1)));
    VV_shifted.fill_bytes();

    return generate_range_proof_power_two_agg(h, VV_shifted, gg, hh, n, m, gamma, vv_shifted);
}


bool b_verify_abs_val_range_power_two_agg(const RistP3AndBytes &h, const RistP3VecAndBytes &VV,
                                          const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                          int n, int m, // m values of range [-2^(n-1), 2^(n-1))
                                          const RangeProofPowerTwoP3 &proof) {
    RistP3VecAndBytes VV_shifted(m);
    VV_shifted.elems = VV.elems + RistElemP3Vec(m, pedersen_zero_commit_p3(power_of_two(n - 1)));

    VV_shifted.fill_bytes();

    return b_verify_range_power_two_agg(h, VV_shifted, gg, hh, n, m, proof);
}

RistScalVec segment(const RistScalVec &rr, int start, int length) {
    return RistScalVec(rr.begin() + start, rr.begin() + start + length);
}

InnerProdProof generate_inner_prod_proof_hh_coeffs(const RistP3VecAndBytes &gg,
                                                   const RistP3VecAndBytes &hh,
                                                   const RistScalVec &hh_coeffs_orig,
                                                   const RistP3AndBytes &u,
                                                   const RistP3AndBytes &P,
                                                   const RistScalVec &aa, const RistScalVec &bb) {
    int n = gg.size();
    assert(n == hh.size());
    assert(n == hh_coeffs_orig.size());
    assert(n == aa.size());
    assert(n == bb.size());

    int k = get_bit_length(n - 1);

    InnerProdProof proof(n);
    if (n == 1) {
        proof.a = aa[0];
        proof.b = bb[0];
        return proof;
    }

    RistScalVec aa_ext = aa;
    aa_ext.reserve(1 << k);
    RistScalVec bb_ext = bb;
    bb_ext.reserve(1 << k);
    for (int i = n; i < (1 << k); i++) {
        aa_ext.emplace_back(c_scal_zero);
        bb_ext.emplace_back(c_scal_zero);
    }

    RistScalVec aa_temp(aa_ext), bb_temp(bb_ext);
    RistScalVec gg_coeffs(n, c_scal_one), hh_coeffs(hh_coeffs_orig);
    RistScalVec xx(k);
    for (int i = 0; i < k; i++) {
        int half_step = 1 << (k - 1 - i);
        auto aa_left = segment(aa_temp, 0, half_step);
        auto aa_right = segment(aa_temp, half_step, half_step);
        auto bb_left = segment(bb_temp, 0, half_step);
        auto bb_right = segment(bb_temp, half_step, half_step);

        auto c_L = inner_prod(aa_left, bb_right);
        auto c_R = inner_prod(aa_right, bb_left);

        proof.LL.elems[i] = linear_comb_trim(block_mult_half(aa_left, gg_coeffs, 1), gg.elems) +
                            linear_comb_trim(block_mult_half(bb_right, hh_coeffs, 0), hh.elems) +
                            c_L * u.elem;
        proof.LL.fill_bytes(i);

        proof.RR.elems[i] = linear_comb_trim(block_mult_half(aa_right, gg_coeffs, 0), gg.elems) +
                            linear_comb_trim(block_mult_half(bb_left, hh_coeffs, 1), hh.elems) +
                            c_R * u.elem;
        proof.RR.fill_bytes(i);

        xx[i] = rist_scalar_from_hash_from_bytes(gg.bytes, hh.bytes, u.bytes, P.bytes,
                                                 proof.LL.segment(i + 1).bytes,
                                                 proof.RR.segment(i + 1).bytes);
        auto x_inv = scalar_invert(xx[i]);

        if (i < k - 1) {
            gg_coeffs = block_mult(x_inv, xx[i], gg_coeffs, half_step);
            hh_coeffs = block_mult(xx[i], x_inv, hh_coeffs, half_step);
        }

        aa_temp = xx[i] * aa_left + x_inv * aa_right;
        bb_temp = x_inv * bb_left + xx[i] * bb_right;
    }
    proof.a = aa_temp[0];
    proof.b = bb_temp[0];

    return proof;
}


InnerProdProof generate_inner_prod_proof(const RistP3VecAndBytes &gg,
                                         const RistP3VecAndBytes &hh,
                                         const RistP3AndBytes &u,
                                         const RistP3AndBytes &P,
                                         const RistScalVec &aa, const RistScalVec &bb) {
    return generate_inner_prod_proof_hh_coeffs(gg, hh, RistScalVec(hh.size(), c_scal_one),
                                               u, P, aa, bb);
}

// returns ret, such that:
// ret[0] is the coeff of gg
// ret[1] is the coeff of hh
// ret[2] is the coeff of proof.LL
// ret[3] is the coeff of proof.RR
// to decide if ret[0] * gg + ret[1] * hh + proof.a * proof.b * u == P + ret[2] * proof.LL + ret[3] * proof.RR
RistScalMat compute_coeffs_inner_prod_proof(const RistP3VecAndBytes &gg,
                                            const RistP3VecAndBytes &hh,
                                            const RistP3AndBytes &u,
                                            const RistP3AndBytes &P,
                                            const InnerProdProof &proof) {
    int n = gg.size();
    assert(n == hh.size());
    int k = get_bit_length(n - 1);
    assert(k == proof.LL.size());

    RistScalMat ret(4);
    if (n == 1) {
        ret[0] = RistScalVec(1, proof.a);
        ret[1] = RistScalVec(1, proof.b);
        return ret;
    }

    RistScalVec xx(k), xx_sq(k), xx_sq_inv(k);
    for (int i = 0; i < k; i++) {
        xx[i] = rist_scalar_from_hash_from_bytes(gg.bytes, hh.bytes, u.bytes, P.bytes,
                                                 proof.LL.segment(i + 1).bytes,
                                                 proof.RR.segment(i + 1).bytes);
        xx_sq[i] = xx[i] * xx[i];
        xx_sq_inv[i] = scalar_invert(xx_sq[i]);
    }
    RistScalVec ss(n), ss_inv(n);
    RistScal temp = xx[0];
    for (int i = 1; i < k; i++) {
        temp *= xx[i];
    }
    ss[0] = scalar_invert(temp);
    for (int i = 1; i < n; i++) {
        int step = get_bit_length(i) - 1;
        ss[i] = ss[i - (1 << step)] * xx_sq[k - 1 - step];
    }
    for (int i = 0; i < n; i++) {
        ss_inv[i] = scalar_invert(ss[i]);
    }
    ret[0] = proof.a * ss;
    ret[1] = proof.b * ss_inv;
    ret[2] = xx_sq;
    ret[3] = xx_sq_inv;
    return ret;
}

bool b_verify_inner_prod_proof(const RistP3VecAndBytes &gg,
                               const RistP3VecAndBytes &hh,
                               const RistP3AndBytes &u,
                               const RistP3AndBytes &P,
                               const InnerProdProof &proof) {
    auto coeffs = compute_coeffs_inner_prod_proof(gg, hh, u, P, proof);
    return linear_comb(coeffs[0], gg.elems) + linear_comb(coeffs[1], hh.elems) +
           proof.a * proof.b * u.elem
           == P.elem + linear_comb(coeffs[2], proof.LL.elems) + linear_comb(coeffs[3], proof.RR.elems);
}

RangeProofPowerTwoP3Log generate_range_proof_power_two_log(const RistP3AndBytes &h, const RistP3AndBytes &V,
                                                           const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                                           int n,
                                                           const RistScal &gamma, const RistScal &v) {
    assert(gg.size() >= n);
    assert(hh.size() >= n);

    RangeProofPowerTwoP3Log proof;

    RistScalVec aa_L(n), aa_R(n);
    RistP3VecAndBytes gg_segment = gg.segment(n);
    RistP3VecAndBytes hh_segment = hh.segment(n);
    for (int b = 0; b < n; b++) {
        aa_L[b] = b_get_bit(v, b) ? c_scal_one : c_scal_zero;
        aa_R[b] = aa_L[b] - c_scal_one;
    }

    RistScal alpha;
    rand_init(alpha);
    LinearCombCalculator lcc_A;
    lcc_A.add(alpha, h.elem);
    lcc_A.add(aa_L, gg_segment.elems);
    lcc_A.add(aa_R, hh_segment.elems);
    proof.A.elem = lcc_A.result();
//    proof.A.elem = alpha * h.elem + linear_comb(aa_L, gg_segment.elems) + linear_comb(aa_R, hh_segment.elems);
    proof.A.fill_bytes();

    RistScalVec ss_L(n), ss_R(n);
    RistScal rho;
    rand_init(ss_L);
    rand_init(ss_R);
    rand_init(rho);
    LinearCombCalculator lcc_S;
    lcc_S.add(rho, h.elem);
    lcc_S.add(ss_L, gg_segment.elems);
    lcc_S.add(ss_R, hh_segment.elems);
    proof.S.elem = lcc_S.result();
//    proof.S.elem = rho * h.elem + linear_comb(ss_L, gg_segment.elems) + linear_comb(ss_R, hh_segment.elems);
    proof.S.fill_bytes();

    RistScal y = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  V.bytes,
                                                  proof.A.bytes, proof.S.bytes);


    RistScal z = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  V.bytes,
                                                  proof.A.bytes, proof.S.bytes,
                                                  y.scalar);

    RistScal tau1, tau2;
    rand_init(tau1);
    rand_init(tau2);
    RistScalVec ll_0 = aa_L - RistScalVec(n, z);
    RistScalVec ll_1 = ss_L;
    RistScalVec rr_0 = scalar_geometric_series(n, y) * (aa_R + RistScalVec(n, z))
                       + z * z * scalar_geometric_series_power_two(n);
    RistScalVec rr_1 = scalar_geometric_series(n, y) * ss_R;
    RistScal t1 = inner_prod(ll_0, rr_1) + inner_prod(ll_1, rr_0);
    RistScal t2 = inner_prod(ll_1, rr_1);
    proof.T1.elem = pedersen_commit(t1, h.elem, tau1);
    proof.T1.fill_bytes();
    proof.T2.elem = pedersen_commit(t2, h.elem, tau2);
    proof.T2.fill_bytes();

    RistScal x = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  V.bytes,
                                                  proof.A.bytes, proof.S.bytes,
                                                  y.scalar, z.scalar,
                                                  proof.T1.bytes, proof.T2.bytes);

    auto ll = ll_0 + x * ll_1;
    auto rr = rr_0 + x * rr_1;
    proof.t_hat = inner_prod(ll, rr);
    proof.tau_x = tau2 * x * x + tau1 * x + z * z * gamma;
    proof.mu = alpha + rho * x;

    auto y_inv = scalar_invert(y);
    auto y_inv_series = scalar_geometric_series(n, y_inv);

//    RistP3AndBytes P;
//    RistElemP3 P;
//    P = proof.mu * h.elem + linear_comb(ll, gg_segment.elems) + linear_comb(rr * y_inv_series, hh_segment.elems);
//    P.fill_bytes();

    RistP3AndBytes inner_prod_commitment;
    LinearCombCalculator lcc_ipc;
    lcc_ipc.add(proof.mu, h.elem);
    lcc_ipc.add(ll, gg_segment.elems);
    lcc_ipc.add(rr * y_inv_series, hh_segment.elems);
    lcc_ipc.add_base(proof.t_hat);
    lcc_ipc.sub(proof.mu, h.elem);
    inner_prod_commitment.elem = lcc_ipc.result();
//    inner_prod_commitment.elem = P + pedersen_commit(proof.t_hat, h.elem, -proof.mu);
    inner_prod_commitment.fill_bytes();

    proof.inner_prod_proof = generate_inner_prod_proof_hh_coeffs(gg_segment, hh_segment, y_inv_series,
                                                                 c_p3_bytes_one, inner_prod_commitment,
                                                                 ll, rr);

    return proof;
}

bool b_verify_range_power_two_log(const RistP3AndBytes &h, const RistP3AndBytes &V,
                                  const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                  int n,
                                  const RangeProofPowerTwoP3Log &proof) {
    assert(gg.size() >= n);
    assert(hh.size() >= n);
    RistP3VecAndBytes gg_segment = gg.segment(n);
    RistP3VecAndBytes hh_segment = hh.segment(n);

    RistScal y = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  V.bytes,
                                                  proof.A.bytes, proof.S.bytes);

    RistScal z = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  V.bytes,
                                                  proof.A.bytes, proof.S.bytes,
                                                  y.scalar);


    RistScal delta = (z - z * z) * (power(n, y) - c_scal_one) / (y - c_scal_one) -
                     (power(n, RistScal(2)) - c_scal_one) * z * z * z;

    RistScal x = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  V.bytes,
                                                  proof.A.bytes, proof.S.bytes,
                                                  y.scalar, z.scalar,
                                                  proof.T1.bytes, proof.T2.bytes);


    if (pedersen_commit(proof.t_hat, h.elem, proof.tau_x) !=
        z * z * V.elem + pedersen_commit(delta, proof.T1.elem, x) + x * x * proof.T2.elem)
        return false;


//    RistElemP3 P;
//    P = proof.A.elem + x * proof.S.elem - z * sum_single_thread(gg_segment.elems) +
//        linear_comb(scalar_geometric_series(n, scalar_invert(y)) *
//                    (z * z * scalar_geometric_series_power_two(n))
//                    + RistScalVec(n, z), hh_segment.elems);
    RistP3AndBytes inner_prod_commitment;
    LinearCombCalculator lcc_ipc;
    lcc_ipc.add(proof.A.elem);
    lcc_ipc.add(x, proof.S.elem);
    lcc_ipc.sub(z, sum_single_thread(gg_segment.elems));
    lcc_ipc.add(scalar_geometric_series(n, scalar_invert(y)) *
                (z * z * scalar_geometric_series_power_two(n))
                + RistScalVec(n, z),
                hh_segment.elems);
    lcc_ipc.add_base(proof.t_hat);
    lcc_ipc.sub(proof.mu, h.elem);
    inner_prod_commitment.elem = lcc_ipc.result();
//    inner_prod_commitment.elem = P + pedersen_commit(proof.t_hat, h.elem, -proof.mu);
    inner_prod_commitment.fill_bytes();

    auto y_inv = scalar_invert(y);
    auto y_inv_series = scalar_geometric_series(n, y_inv);

    auto coeffs = compute_coeffs_inner_prod_proof(gg_segment, hh_segment, c_p3_bytes_one, inner_prod_commitment,
                                                  proof.inner_prod_proof);

    LinearCombCalculator lcc;
    lcc.add(coeffs[0], gg_segment.elems);
    lcc.add(coeffs[1] * y_inv_series, hh_segment.elems);
    lcc.add_base(proof.inner_prod_proof.a * proof.inner_prod_proof.b);

//    lcc.sub(proof.A.elem);
//    lcc.sub(x, proof.S.elem);
//    lcc.add(z, sum_single_thread(gg_segment.elems));
//    lcc.sub(scalar_geometric_series(n, scalar_invert(y)) *
//                (z * z * scalar_geometric_series_power_two(n))
//                + RistScalVec(n, z),
//                hh_segment.elems);
//    lcc.sub_base(proof.t_hat);
//    lcc.add(proof.mu, h.elem);
    lcc.sub(inner_prod_commitment.elem);

    lcc.sub(coeffs[2], proof.inner_prod_proof.LL.elems);
    lcc.sub(coeffs[3], proof.inner_prod_proof.RR.elems);
    return lcc.is_0();

//    return linear_comb(coeffs[0], gg_segment.elems) + linear_comb(coeffs[1] * y_inv_series, hh_segment.elems) +
//           proof.inner_prod_proof.a * proof.inner_prod_proof.b * c_p3_bytes_one.elem
//           == inner_prod_commitment.elem + linear_comb(coeffs[2], proof.inner_prod_proof.LL.elems) + linear_comb(coeffs[3], proof.inner_prod_proof.RR.elems);
}

RangeProofPowerTwoP3Log generate_range_proof_power_two_agg_log(const RistP3AndBytes &h, const RistP3VecAndBytes &VV,
                                                               const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                                               int n, int m, // m values of range [0, 2^n)
                                                               const RistScalVec &gamma, const RistScalVec &vv) {
    assert(m == VV.size());
    assert(m == vv.size());
    assert(m == gamma.size());

    assert(gg.size() >= n * m);
    assert(hh.size() >= n * m);

    RangeProofPowerTwoP3Log proof;
    RistScalVec aa_L(n * m), aa_R(n * m);
    RistP3VecAndBytes gg_segment = gg.segment(n * m);
    RistP3VecAndBytes hh_segment = hh.segment(n * m);
    for (int j = 0; j < m; j++) {
        for (int b = 0; b < n; b++) {
            aa_L[j * n + b] = b_get_bit(vv[j], b) ? c_scal_one : c_scal_zero;
            aa_R[j * n + b] = aa_L[j * n + b] - c_scal_one;
        }
    }
    RistScal alpha;
    rand_init(alpha);
    LinearCombCalculator lcc_A;
    lcc_A.add(alpha, h.elem);
    lcc_A.add(aa_L, gg_segment.elems);
    lcc_A.add(aa_R, hh_segment.elems);
    proof.A.elem = lcc_A.result();
//    proof.A.elem = alpha * h.elem + linear_comb(aa_L, gg_segment.elems) + linear_comb(aa_R, hh_segment.elems);
    proof.A.fill_bytes();

    RistScalVec ss_L(n * m), ss_R(n * m);
    RistScal rho;
    rand_init(ss_L);
    rand_init(ss_R);
    rand_init(rho);
    LinearCombCalculator lcc_S;
    lcc_S.add(rho, h.elem);
    lcc_S.add(ss_L, gg_segment.elems);
    lcc_S.add(ss_R, hh_segment.elems);
    proof.S.elem = lcc_S.result();
//    proof.S.elem = rho * h.elem + linear_comb(ss_L, gg_segment.elems) + linear_comb(ss_R, hh_segment.elems);
    proof.S.fill_bytes();

    RistScal y = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  VV.bytes,
                                                  proof.A.bytes, proof.S.bytes);


    RistScal z = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  VV.bytes,
                                                  proof.A.bytes, proof.S.bytes,
                                                  y.scalar);

    RistScal tau1, tau2;
    rand_init(tau1);
    rand_init(tau2);
    RistScalVec ll_0 = aa_L - RistScalVec(n * m, z);
    RistScalVec ll_1 = ss_L;

    RistScalVec rr_0 = scalar_geometric_series(n * m, y) * (aa_R + RistScalVec(n * m, z));
    RistScalVec power_segments;
    RistScalVec segment = z * scalar_geometric_series_power_two(n);
    for (int j = 0; j < m; j++) {
        segment = z * segment;
        power_segments.insert(power_segments.end(), segment.begin(), segment.end());
    }
    assert(power_segments.size() == n * m);
    rr_0 = rr_0 + power_segments;

    RistScalVec rr_1 = scalar_geometric_series(n * m, y) * ss_R;
    RistScal t1 = inner_prod(ll_0, rr_1) + inner_prod(ll_1, rr_0);
    RistScal t2 = inner_prod(ll_1, rr_1);
    proof.T1.elem = pedersen_commit(t1, h.elem, tau1);
    proof.T1.fill_bytes();
    proof.T2.elem = pedersen_commit(t2, h.elem, tau2);
    proof.T2.fill_bytes();

    RistScal x = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  VV.bytes,
                                                  proof.A.bytes, proof.S.bytes,
                                                  y.scalar, z.scalar,
                                                  proof.T1.bytes, proof.T2.bytes);


    auto ll = ll_0 + x * ll_1;
    auto rr = rr_0 + x * rr_1;
    proof.t_hat = inner_prod(ll, rr);

    proof.tau_x = tau2 * x * x + tau1 * x;
    RistScal z_gamma_power_sum = c_scal_zero;
    RistScal power_of_z = z;
    for (int j = 0; j < m; j++) {
        power_of_z *= z;
        z_gamma_power_sum += power_of_z * gamma[j];
    }
    proof.tau_x += z_gamma_power_sum;

    proof.mu = alpha + rho * x;

    auto y_inv = scalar_invert(y);
    auto y_inv_series = scalar_geometric_series(n * m, y_inv);

//    RistP3AndBytes P;
//    P.elem = proof.mu * h.elem + linear_comb(ll, gg_segment.elems) + linear_comb(rr * y_inv_series, hh_segment.elems);
//    P.fill_bytes();

    RistP3AndBytes inner_prod_commitment;
    LinearCombCalculator lcc_ipc;
    lcc_ipc.add(proof.mu, h.elem);
    lcc_ipc.add(ll, gg_segment.elems);
    lcc_ipc.add(rr * y_inv_series, hh_segment.elems);
    lcc_ipc.add_base(proof.t_hat);
    lcc_ipc.sub(proof.mu, h.elem);
    inner_prod_commitment.elem = lcc_ipc.result();
//    inner_prod_commitment.elem = P.elem + pedersen_commit(proof.t_hat, h.elem, -proof.mu);
    inner_prod_commitment.fill_bytes();

    proof.inner_prod_proof = generate_inner_prod_proof_hh_coeffs(gg_segment, hh_segment, y_inv_series,
                                                                 c_p3_bytes_one, inner_prod_commitment,
                                                                 ll, rr);

    return proof;
}


bool b_verify_range_power_two_agg_log(const RistP3AndBytes &h, const RistP3VecAndBytes &VV,
                                      const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                      int n, int m, // m values of range [0, 2^n)
                                      const RangeProofPowerTwoP3Log &proof) {
    assert(gg.size() >= n * m);
    assert(hh.size() >= n * m);

    RistP3VecAndBytes gg_segment = gg.segment(n * m);
    RistP3VecAndBytes hh_segment = hh.segment(n * m);
    assert(VV.size() == m);

    RistScal y = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  VV.bytes,
                                                  proof.A.bytes, proof.S.bytes);

    RistScal z = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  VV.bytes,
                                                  proof.A.bytes, proof.S.bytes,
                                                  y.scalar);


    RistScal delta = (z - z * z) * (power(n * m, y) - c_scal_one) / (y - c_scal_one) -
                     (power(n, RistScal(2)) - c_scal_one) * z * z * z * (power(m, z) - c_scal_one) / (z - c_scal_one);

    RistScal x = rist_scalar_from_hash_from_bytes(h.bytes, gg_segment.bytes, hh_segment.bytes,
                                                  VV.bytes,
                                                  proof.A.bytes, proof.S.bytes,
                                                  y.scalar, z.scalar,
                                                  proof.T1.bytes, proof.T2.bytes);

    if (pedersen_commit(proof.t_hat, h.elem, proof.tau_x) !=
        z * z * linear_comb(scalar_geometric_series(m, z), VV.elems)
        +
        pedersen_commit(delta, proof.T1.elem, x) + x * x * proof.T2.elem)
        return false;

    RistElemP3 P;
    RistScalVec power_segments;
    RistScalVec segment = z * scalar_geometric_series_power_two(n);

    for (int j = 0; j < m; j++) {
        segment = z * segment;
        power_segments.insert(power_segments.end(), segment.begin(), segment.end());
    }
    assert(power_segments.size() == n * m);
    RistScalVec coeffs_temp = scalar_geometric_series(n * m, scalar_invert(y)) * power_segments + RistScalVec(n * m, z);

    P = proof.A.elem + x * proof.S.elem - z * sum_single_thread(gg_segment.elems) +
        linear_comb(coeffs_temp, hh_segment.elems);

    RistP3AndBytes inner_prod_commitment;
    LinearCombCalculator lcc_ipc;
    lcc_ipc.add(proof.A.elem);
    lcc_ipc.add(x, proof.S.elem);
    lcc_ipc.sub(z, sum_single_thread(gg_segment.elems));
    lcc_ipc.add(coeffs_temp,
                hh_segment.elems);
    lcc_ipc.add_base(proof.t_hat);
    lcc_ipc.sub(proof.mu, h.elem);
    inner_prod_commitment.elem = lcc_ipc.result();

//    inner_prod_commitment.elem = P + pedersen_commit(proof.t_hat, h.elem, -proof.mu);
    inner_prod_commitment.fill_bytes();

    auto y_inv = scalar_invert(y);
    auto y_inv_series = scalar_geometric_series(n * m, y_inv);

    auto coeffs = compute_coeffs_inner_prod_proof(gg_segment, hh_segment, c_p3_bytes_one, inner_prod_commitment,
                                                  proof.inner_prod_proof);

    LinearCombCalculator lcc;
    lcc.add(coeffs[0], gg_segment.elems);
    lcc.add(coeffs[1] * y_inv_series, hh_segment.elems);
    lcc.add_base(proof.inner_prod_proof.a * proof.inner_prod_proof.b);
    lcc.sub(inner_prod_commitment.elem);
    lcc.sub(coeffs[2], proof.inner_prod_proof.LL.elems);
    lcc.sub(coeffs[3], proof.inner_prod_proof.RR.elems);

    return lcc.is_0();

//    return linear_comb(coeffs[0], gg_segment.elems) + linear_comb(coeffs[1] * y_inv_series, hh_segment.elems) +
//           proof.inner_prod_proof.a * proof.inner_prod_proof.b * c_p3_bytes_one.elem
//           == inner_prod_commitment.elem + linear_comb(coeffs[2], proof.inner_prod_proof.LL.elems) + linear_comb(coeffs[3], proof.inner_prod_proof.RR.elems);
}

RangeProofPowerTwoP3Log
generate_abs_val_range_proof_power_two_agg_log(const RistP3AndBytes &h, const RistP3VecAndBytes &VV,
                                               const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                               int n, int m, // m values of range [-2^(n-1), 2^(n-1))
                                               const RistScalVec &gamma, const RistScalVec &vv) {
    RistScalVec vv_shifted = vv + RistScalVec(m, power_of_two(n - 1));

    RistP3VecAndBytes VV_shifted(m);
    VV_shifted.elems = VV.elems + RistElemP3Vec(m, pedersen_zero_commit_p3(power_of_two(n - 1)));
    VV_shifted.fill_bytes();

    return generate_range_proof_power_two_agg_log(h, VV_shifted, gg, hh, n, m, gamma, vv_shifted);
}


bool b_verify_abs_val_range_power_two_agg_log(const RistP3AndBytes &h, const RistP3VecAndBytes &VV,
                                              const RistP3VecAndBytes &gg, const RistP3VecAndBytes &hh,
                                              int n, int m, // m values of range [-2^(n-1), 2^(n-1))
                                              const RangeProofPowerTwoP3Log &proof) {
    RistP3VecAndBytes VV_shifted(m);
    VV_shifted.elems = VV.elems + RistElemP3Vec(m, pedersen_zero_commit_p3(power_of_two(n - 1)));
    VV_shifted.fill_bytes();
    return b_verify_range_power_two_agg_log(h, VV_shifted, gg, hh, n, m, proof);
}