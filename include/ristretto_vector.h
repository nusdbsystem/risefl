//
// Created by yizheng on 7/3/23.
//

// ristretto_vector.h
// 
// Vectorized ristretto group operations.

#ifndef RISEFL_CRYPTO_RISTRETTO_VECTOR_H
#define RISEFL_CRYPTO_RISTRETTO_VECTOR_H

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

// convert the vector of scalar to bytes, write the bytes to a vector starting from the place pointed by the iterator it, and move the iterator it to one place after the end
void export_to_bytestream(const RistScalVec &rr, std::vector<unsigned char>::iterator &it);

// read bytes from a vector starting from the place pointed by the iterator it, convert the bytes to the scalar of bytes, and move the iterator it to one place after the end
void import_from_bytestream(RistScalVec &rr, std::vector<unsigned char>::const_iterator &it);

typedef std::vector<RistScalVec> RistScalMat;

// convert the matrix (i.e. double vector) of scalar to bytes, write the bytes to a vector starting from the place pointed by the iterator it, and move the iterator it to one place after the end
void export_to_bytestream(const RistScalMat &rr, std::vector<unsigned char>::iterator &it);

// read bytes from a vector starting from the place pointed by the iterator it, convert the bytes to the matrix (i.e. double vector) of bytes, and move the iterator it to one place after the end
void import_from_bytestream(RistScalMat &rr, std::vector<unsigned char>::const_iterator &it);

typedef std::vector<RistHashbytes> RistHashbytesVec;
typedef std::vector<RistNonreducedScalar> RistNonreducedScalarVec;

// convert a vector of integers to a vector of scalars
RistScalVec rist_scal_vec_from_int_vec(const std::vector<int> &s);
RistScalVec rist_scal_vec_from_long_vec(const std::vector<long> &s);

// convert a vector of group elements to a vector of bytes
std::vector<unsigned char> bytes_from_rist_elem_vec(const RistElemVec &gg);

// convert a vector of scalars to a vector of bytes
std::vector<unsigned char> bytes_from_rist_scal_vec(const RistScalVec &gg);

// convert a vector of bytes to a vector of group elements
RistElemVec rist_elem_vec_from_bytes(const std::vector<unsigned char> &vv);

// convert a vector of bytes to a vector of scalars
RistScalVec rist_scal_vec_from_bytes(const std::vector<unsigned char> &vv);

// randomly initialize every element of a vector of group elements
void rand_init(RistElemVec &pp);

// initialize initialize every element of a vector of group elements from the corresponding member in the vector of hashes
void hash_init(RistElemVec &pp, const RistHashbytesVec &rr);

// pointwise add pp and qq, write the result into rr
void add(RistElemVec &rr, const RistElemVec &pp, const RistElemVec &qq);

// pointwise multiply nn by pp, write the result into qq
void scalar_mult(RistElemVec &qq, const RistScalVec &nn, const RistElemVec pp);

// pointwise multiply
RistElemVec operator*(const RistScalVec &nn, const RistElemVec pp);

// multiply the fixed group generator by every element of nn, write them into qq
void scalar_mult_base(RistElemVec &qq, const RistScalVec &nn);

// randomly initialize every element of a vector of scalars
void rand_init(RistScalVec &pp);

// pointwise reduce
void scalar_reduce(RistScalVec &rr, const RistNonreducedScalarVec &ss);

// pointwise reduce
void scalar_reduce(RistScalVec &rr, const RistHashbytesVec &ss);

// pointwise add xx and yy, write the result into zz
void scalar_add(RistScalVec &zz, const RistScalVec &xx, const RistScalVec &yy);

// pointwise add
RistScalVec operator+(const RistScalVec &xx, const RistScalVec &yy);

// pointwise substract xx by yy, write the result into zz
void scalar_sub(RistScalVec &zz, const RistScalVec &xx, const RistScalVec &yy);

// pointwise substract
RistScalVec operator-(const RistScalVec &xx, const RistScalVec &yy);

// pointwise take the additive inverse
RistScalVec operator-(const RistScalVec &xx);

// pointwise multiply xx by yy, write the result into zz
void scalar_mul(RistScalVec &zz, const RistScalVec &xx, const RistScalVec &yy);

// pointwise multiply
RistScalVec operator*(const RistScalVec &xx, const RistScalVec &yy);

// multiply x by every every element of yy (i.e. broadcasting)
RistScalVec operator*(const RistScal &x, const RistScalVec &yy);

// compute the sum of every group element in pp
RistElem sum(const RistElemVec &pp);

// compute the sum of every scalar in pp
RistScal sum(const RistScalVec &pp);

// possible optimizations of linear combination:
// Pippenger'scalar_seed approach "On the evaluation of powers and related problems"
//      http://ieeexplore.ieee.org/document/4567910/
// and "On the evaluation of powers and monomials"
//      https://scholarship.claremont.edu/cgi/viewcontent.cgi?article=1141&context=hmc_fac_pub
//  which are mentioned in "Faster Batch Forgery Identification"
//      https://cr.yp.to/badbatch/badbatch-20120919.pdf

// These are implemented in rist_fast_computations.h, rist_fast_computations.cpp

// nn[0] * pp[0] + ... + nn[nn.size() - 1] * pp[nn.size() - 1]
// requirement: nn.size() == pp.size()
RistElem linear_comb(const RistScalVec &nn, const RistElemVec &pp);

// multiply r by every every element of B (i.e. broadcasting)
RistScalMat operator*(const RistScal &r, const RistScalMat &B);

// pointwise add
RistScalMat operator+(const RistScalMat &A, const RistScalMat &B);

// pointwise take additive inverse
RistScalMat operator-(const RistScalMat &A);


#endif //RISEFL_CRYPTO_RISTRETTO_VECTOR_H
