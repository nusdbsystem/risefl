//
// Created by yizheng on 10/5/23.
//

#ifndef RANDOM_AVX_MERSENNE_TWISTER_AVX2_H
#define RANDOM_AVX_MERSENNE_TWISTER_AVX2_H

#include <stdint.h>
#include "immintrin.h"
//
//#ifdef __cplusplus
//extern "C" {
//#endif


class mt_avx2 {
private:
    static const size_t SIZE = 624;
    static const size_t PERIOD = 397;
    static const size_t DIFF = SIZE - PERIOD;

    inline static const __m256i m256_MAGIC = _mm256_set1_epi32(0x9908b0df);
    inline static const __m256i m256_0x6c078965 = _mm256_set1_epi32(0x6c078965);
    inline static const __m256i m256_0x9d2c5680 = _mm256_set1_epi32(0x9d2c5680);
    inline static const __m256i m256_0xefc60000 = _mm256_set1_epi32(0xefc60000);
    inline static const __m256i m256_0x80000000 = _mm256_set1_epi32(0x80000000);
    inline static const __m256i m256_0x7FFFFFFF = _mm256_set1_epi32(0x7FFFFFFF);
    inline static const __m256i m256_0x1 = _mm256_set1_epi32(0x1);

    constexpr static float float_multiplier = 1.0f / (1L << 31);
    inline static const __m256 m256_float_multiplier = _mm256_broadcast_ss(&float_multiplier);
    constexpr static double float64_multiplier = 1.0 / (1L << 32);
    inline static const __m256d m256d_float64_multiplier = _mm256_broadcast_sd(&float64_multiplier);


// State for a singleton Mersenne Twister. If you want to make this into a
// class, these are what you need to isolate.

    __m256i MT[SIZE];
    __m256i MT_TEMPERED[SIZE];
    size_t index = SIZE;

    void generate_numbers();

    // Extract 8 uint32 values in [0,1)
    __m256i rand_u32_avx2();

    // Extract 8 float32 values in [0,1)
    __m256 rand_float_avx2();

    // Extract 8 float64 values in [0,1)
    void rand_float64_avx2(__m256d *dest1, __m256d *dest2);

public:
    mt_avx2(const uint32_t *seed_value_pointer) {
        seed_avx2(seed_value_pointer);
    }

    mt_avx2(const uint32_t seed_value) {
        uint32_t seed_value_pointer[8];
        for (int i = 0; i < 8; i++) {
            seed_value_pointer[i] = seed_value + i;
        }
        seed_avx2(seed_value_pointer);
    }

/*
 * Initialize Mersenne Twister with given seed_value_pointer value.
 */
    void seed_avx2(const uint32_t *seed_value_pointer);

    /*
* Extract 8 pseudo-random unsigned 32-bit integers starting with dest in the range 0 ... UINT32_MAX
*/
    void rand_u32_avx2(uint32_t *dest);

    // Extract 8 float32 values in [0,1)
    void rand_float_avx2(float *dest);

    // Extract 8 float64 values in [0,1)
    void rand_float64_avx2(double *dest);

    // Extract 16 float32 pseudo-random normal samples from N(0,1)
    void rand_normal_avx2(float *dest);

    // Extract 8 float64 pseudo-random normal samples from N(0,1)
    void rand_normal_avx2(double *dest);

    // Extract count float32 pseudo-random normal samples from N(0,1)
    void rand_normal_avx2(float *dest, int count);

    // Extract count float64 pseudo-random normal samples from N(0,1)
    void rand_normal_avx2(double *dest, int count);

    void rand_normal_avx2_shifted(int *dest, int count, int bit_shifter);

};

//#ifdef __cplusplus
//} // extern "C"
//#endif


#endif //RANDOM_AVX_MERSENNE_TWISTER_AVX2_H
