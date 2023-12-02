//
// Created by yizheng on 7/3/23.
//

#ifndef DI_ZKP_CRYPTO_RISTRETTO_VECTOR_H
#define DI_ZKP_CRYPTO_RISTRETTO_VECTOR_H

#include <vector>
#include <array>
#include <sodium.h>
#include <algorithm>
#include <execution>
#include <thread>
#include <chrono>
#include "ristretto.h"

//using RistElemVec = std::vector<RistElem>;
//using RistScalVec = std::vector<RistScal>;
typedef std::vector<RistElem> RistElemVec;
typedef std::vector<RistElemVec> RistElemMat;
typedef std::vector<RistScal> RistScalVec;

void export_to_bytestream(const RistScalVec &rr, std::vector<unsigned char>::iterator &it);

void import_from_bytestream(RistScalVec &rr, std::vector<unsigned char>::const_iterator &it);

typedef std::vector<RistScalVec> RistScalMat;

void export_to_bytestream(const RistScalMat &rr, std::vector<unsigned char>::iterator &it);

void import_from_bytestream(RistScalMat &rr, std::vector<unsigned char>::const_iterator &it);

typedef std::vector<RistHashbytes> RistHashbytesVec;
typedef std::vector<RistNonreducedScalar> RistNonreducedScalarVec;

RistScalVec rist_scal_vec_from_int_vec(const std::vector<int> &s);
RistScalVec rist_scal_vec_from_long_vec(const std::vector<long> &s);

std::vector<unsigned char> bytes_from_rist_elem_vec(const RistElemVec &gg);

std::vector<unsigned char> bytes_from_rist_scal_vec(const RistScalVec &gg);

RistElemVec rist_elem_vec_from_bytes(const std::vector<unsigned char> &vv);

RistScalVec rist_scal_vec_from_bytes(const std::vector<unsigned char> &vv);


void rand_init(RistElemVec &pp);

void hash_init(RistElemVec &pp, const RistHashbytesVec &rr);

void add(RistElemVec &rr, const RistElemVec &pp, const RistElemVec &qq);


void scalar_mult(RistElemVec &qq, const RistScalVec &nn, const RistElemVec pp);

RistElemVec operator*(const RistScalVec &nn, const RistElemVec pp);

void scalar_mult_base(RistElemVec &qq, const RistScalVec &nn);

void rand_init(RistScalVec &pp);

void scalar_reduce(RistScalVec &rr, const RistNonreducedScalarVec &ss);

void scalar_reduce(RistScalVec &rr, const RistHashbytesVec &ss);

void scalar_add(RistScalVec &zz, const RistScalVec &xx, const RistScalVec &yy);

RistScalVec operator+(const RistScalVec &xx, const RistScalVec &yy);

void scalar_sub(RistScalVec &zz, const RistScalVec &xx, const RistScalVec &yy);

RistScalVec operator-(const RistScalVec &xx, const RistScalVec &yy);

RistScalVec operator-(const RistScalVec &xx);

void scalar_mul(RistScalVec &zz, const RistScalVec &xx, const RistScalVec &yy);

RistScalVec operator*(const RistScalVec &xx, const RistScalVec &yy);

RistScalVec operator*(const RistScal &x, const RistScalVec &yy);

RistElem sum(const RistElemVec &pp);

RistScal sum(const RistScalVec &pp);

// possible optimizations of linear combination:
// Pippenger'scalar_seed approach "On the evaluation of powers and related problems"
//      http://ieeexplore.ieee.org/document/4567910/
// and "On the evaluation of powers and monomials"
//      https://scholarship.claremont.edu/cgi/viewcontent.cgi?article=1141&context=hmc_fac_pub
//  which are mentioned in "Faster Batch Forgery Identification"
//      https://cr.yp.to/badbatch/badbatch-20120919.pdf

// These are implemented in rist_fast_computations.h, rist_fast_computations.cpp

RistElem linear_comb(const RistScalVec &nn, const RistElemVec &pp);


RistScalMat operator*(const RistScal &r, const RistScalMat &B);

RistScalMat operator+(const RistScalMat &A, const RistScalMat &B);

RistScalMat operator-(const RistScalMat &A);


#endif //DI_ZKP_CRYPTO_RISTRETTO_VECTOR_H
