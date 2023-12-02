//
// Created by yizheng on 10/5/23.
//
#include <stdio.h>
#include <iostream>
#include <cmath>
#include "mersenne-twister-avx2.h"
#include "immintrin.h"
#include "../vectorclass/version2-master/vectormath_exp.h"
#include "../vectorclass/version2-master/vectormath_trig.h"


void mt_avx2::generate_numbers() {
    size_t i = 0;
    __m256i y;

    while (i < DIFF) {
        y = _mm256_or_si256(_mm256_and_si256(m256_0x80000000, MT[i]), _mm256_and_si256(m256_0x7FFFFFFF, MT[i + 1])); \
  MT[i] = _mm256_xor_si256(_mm256_xor_si256(MT[i + PERIOD], _mm256_srli_epi32(y, 1)), \
    _mm256_mullo_epi32(m256_MAGIC, _mm256_and_si256(y, m256_0x1))); \
  ++i;
    }

    while (i < SIZE - 1) {
        y = _mm256_or_si256(_mm256_and_si256(m256_0x80000000, MT[i]), _mm256_and_si256(m256_0x7FFFFFFF, MT[i + 1])); \
  MT[i] = _mm256_xor_si256(_mm256_xor_si256(MT[i - DIFF], _mm256_srli_epi32(y, 1)), \
    _mm256_mullo_epi32(m256_MAGIC, _mm256_and_si256(y, m256_0x1))); \
  ++i;
    }
    {
        y = _mm256_or_si256(_mm256_and_si256(m256_0x80000000, MT[SIZE - 1]), _mm256_and_si256(m256_0x7FFFFFFF, MT[0])); \
  MT[SIZE - 1] = _mm256_xor_si256(_mm256_xor_si256(MT[PERIOD - 1], _mm256_srli_epi32(y, 1)), \
    _mm256_mullo_epi32(m256_MAGIC, _mm256_and_si256(y, m256_0x1))); \

    }

    for (size_t i = 0; i < SIZE; ++i) {
        y = MT[i];
        y = _mm256_xor_si256(y, _mm256_srli_epi32(y, 11));
        y = _mm256_xor_si256(y, _mm256_and_si256(_mm256_slli_epi32(y, 7), m256_0x9d2c5680));
        y = _mm256_xor_si256(y, _mm256_and_si256(_mm256_slli_epi32(y, 15), m256_0xefc60000));
        y = _mm256_xor_si256(y, _mm256_srli_epi32(y, 18));
        MT_TEMPERED[i] = y;
    }

    index = 0;
}

//    extern "C"
void mt_avx2::seed_avx2(const uint32_t *seed_value_pointer) {
    __m256i v = _mm256_loadu_si256((const __m256i *) seed_value_pointer);
    _mm256_storeu_si256(&MT[0], v);
//        state.MT[0] = _mm256_loadu_si256((const __m256i *) seed_value_pointer);
    index = SIZE;

    for (uint_fast32_t i = 1; i < SIZE; ++i) {
        _mm256_storeu_si256(&MT[i],
                            _mm256_add_epi32(_mm256_mullo_epi32(m256_0x6c078965,
                                                                _mm256_xor_si256(MT[i - 1],
                                                                                 _mm256_srli_epi32(MT[i - 1],
                                                                                                   30))),
                                             _mm256_set1_epi32(i)));
    }

//        uint32_t temp_int[8];
//        _mm256_storeu_si256((__m256i *)temp_int, MT[0]);
//        std::cout << "avx seed finished, MT[0] = " << temp_int[0] << std::endl;
//        _mm256_storeu_si256((__m256i *)temp_int, MT[1]);
//        std::cout << "avx seed finished, MT[1] = " << temp_int[0] << std::endl;
//        _mm256_storeu_si256((__m256i *)temp_int, MT[10]);
//        std::cout << "avx seed finished, MT[10] = " << temp_int[0] << std::endl;

}


//extern "C"
__m256i mt_avx2::rand_u32_avx2() {
    if (index == SIZE) {
        generate_numbers();
        index = 0;
    }
    return MT_TEMPERED[index++];
}

//extern "C"
void mt_avx2::rand_u32_avx2(uint32_t *dest) {
    _mm256_storeu_si256((__m256i *) dest, rand_u32_avx2());
}

__m256 mt_avx2::rand_float_avx2() {
    if (index == SIZE) {
        generate_numbers();
        index = 0;
    }
    __m256i temp_uint = _mm256_srli_epi32(MT_TEMPERED[index++], 1); // shift right by 1 to remove sign
    __m256 temp_float = _mm256_cvtepi32_ps(temp_uint);
    temp_float = _mm256_mul_ps(temp_float, m256_float_multiplier); // then multiply by 2^(-31) to convert to [0,1)
    return temp_float;
}


void mt_avx2::rand_float64_avx2(__m256d *dest1, __m256d *dest2){
    if (index == SIZE) {
        generate_numbers();
        index = 0;
    }
    __m256i temp_uint = MT_TEMPERED[index++];
    __m128i* temp_uint1_pt = (__m128i *) &temp_uint;
    __m128i* temp_uint2_pt = temp_uint1_pt + 1;
    __m256d temp_double_1 = _mm256_cvtepi32_pd(*temp_uint1_pt);
    __m256d temp_double_2 = _mm256_cvtepi32_pd(*temp_uint2_pt);

    const __m256d half = _mm256_set1_pd(0.5);

    (*dest1) = _mm256_fmadd_pd(temp_double_1, m256d_float64_multiplier, half);
    (*dest2) = _mm256_fmadd_pd(temp_double_2, m256d_float64_multiplier, half);
//    (*dest1) = _mm256_add_pd(_mm256_mul_pd(temp_double_1, m256d_float64_multiplier), half);
//    (*dest2) = _mm256_add_pd(_mm256_mul_pd(temp_double_2, m256d_float64_multiplier), half);
}


void mt_avx2::rand_float_avx2(float *dest) {
    _mm256_storeu_ps(dest, rand_float_avx2());
}

void mt_avx2::rand_float64_avx2(double *dest) {
    __m256d num1, num2;
    rand_float64_avx2(&num1, &num2);
    _mm256_storeu_pd(dest, num1);
    _mm256_storeu_pd(dest + 4, num2);
}

void mt_avx2::rand_normal_avx2(float *dest) {
    const __m256 twopi = _mm256_set1_ps(2.0f * 3.14159265358979323846f);
    const __m256 one = _mm256_set1_ps(1.0f);
    const __m256 minustwo = _mm256_set1_ps(-2.0f);

    __m256 u1 = _mm256_sub_ps(one, rand_float_avx2());
    __m256 u2 = rand_float_avx2();
    __m256 logu1 = __m256(log(Vec8f(u1)));
    __m256 radius = _mm256_sqrt_ps(_mm256_mul_ps(minustwo, logu1));
//    __m256 radius = _mm256_sqrt_ps(_mm256_mul_ps(minustwo, log256_ps(u1)));
    __m256 theta = _mm256_mul_ps(twopi, u2);
    __m256 sintheta, costheta;
//    sincos256_ps(theta, &sintheta, &costheta);
    Vec8f costheta_vec8f;
    Vec8f sintheta_vec8f = sincos(&costheta_vec8f, Vec8f(theta));
    sintheta = __m256(sintheta_vec8f);
    costheta = __m256(costheta_vec8f);

    _mm256_storeu_ps(dest, _mm256_mul_ps(radius, costheta));
    _mm256_storeu_ps(dest + 8, _mm256_mul_ps(radius, sintheta));
}

void mt_avx2::rand_normal_avx2(double *dest) {
    const __m256d twopi = _mm256_set1_pd(2.0 * 3.14159265358979323846);
    const __m256d one = _mm256_set1_pd(1.0);
    const __m256d minustwo = _mm256_set1_pd(-2.0);

    __m256d u1, u2;
    rand_float64_avx2(&u1, &u2);
    u1 = _mm256_sub_pd(one, u1);

    __m256d logu1 = __m256d(log(Vec4d(u1)));
    __m256d radius = _mm256_sqrt_pd(_mm256_mul_pd(minustwo, logu1));
    __m256d theta = _mm256_mul_pd(twopi, u2);
    Vec4d costheta_vec4d;
    Vec4d sintheta_vec4d = sincos(&costheta_vec4d, Vec4d(theta));
    __m256d sintheta, costheta;
    sintheta = __m256d(sintheta_vec4d);
    costheta = __m256d(costheta_vec4d);
    _mm256_storeu_pd(dest, _mm256_mul_pd(radius, costheta));
    _mm256_storeu_pd(dest + 4, _mm256_mul_pd(radius, sintheta));

}

    void mt_avx2::rand_normal_avx2(float *dest, int count) {
    int i = 0;
    for (; i + 16 < count; i += 16) {
        rand_normal_avx2(dest + i);
    }
    float temp[16];
    rand_normal_avx2(temp);
    for (int j = i; j < count; j++)
        dest[j] = temp[j - i];
}

void mt_avx2::rand_normal_avx2(double *dest, int count) {
    int i = 0;
    for (; i + 8 < count; i += 8) {
        rand_normal_avx2(dest + i);
    }
    double temp[8];
    rand_normal_avx2(temp);
    for (int j = i; j < count; j++)
        dest[j] = temp[j - i];
}

void mt_avx2::rand_normal_avx2_shifted(int *dest, int count, int bit_shifter) {
    if (bit_shifter <= 20) {
        float temp[count];
        rand_normal_avx2(temp, count);
        for (int i = 0; i < count; i++) {
            dest[i] = std::lround(temp[i] * (1 << bit_shifter));
        }
    }
    else {
        double temp[count];
        rand_normal_avx2(temp, count);
        for (int i = 0; i < count; i++) {
            dest[i] = std::lround(temp[i] * (1L << bit_shifter));
        }
    }
}