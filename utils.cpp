//
// Created by yizheng on 10/4/23.
//

#include "include/utils.h"

#include <functional>
#include <chrono>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <execution>
#include <random>
#include "immintrin.h"

double measure_time(std::function<void(void)> f) {
    typedef std::chrono::high_resolution_clock Clock;
    auto t1 = Clock::now();
    f();
    auto t2 = Clock::now();
    long double mytime;
    return std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / double(1000000000);
}


void measure_time_and_add_to_bench(std::function<void(void)> f, std::shared_ptr<double> cost) {
    if (cost != nullptr) {
        *cost += measure_time(f);
    } else {
        f();
    }
}

std::vector<unsigned char> bytes_from_int(int c) {
    std::vector<unsigned char> b(sizeof c);
    memcpy(b.data(), &c, sizeof c);
    return b;
}

bool b_get_bit(const RistScal &p, int n) {
    if (n >= 255)
        return false;
    return (p.scalar[n / 8] >> (n % 8)) & (unsigned char) 1;
}

int get_power_two_bound(const RistScal &bound) {
    int n = 255; // largest bit of bound is ignored. see the end of https://doc.libsodium.org/advanced/point-arithmetic/ristretto
    assert (bound != c_scal_zero);
    while (!b_get_bit(bound, n - 1)) {
        n--;
    }
    return n;
}


RistScalVec scalar_geometric_series(int n, const RistScal &p) {
    RistScalVec series(n);
    series[0] = c_scal_one;
    for (int i = 1; i < n; i++) {
        series[i] = series[i - 1] * p;
    }
    return series;
}

RistScal inner_prod(const std::vector<long> &aa, const std::vector<long> &bb) {
    // long ret = 0;
    // assert(aa.size() == bb.size());
    // for (int i = 0; i < aa.size(); i++) {
    //     ret += aa[i] * bb[i];
    // }
    // return ret;


    MyMpz ret_mpz;
    assert(aa.size() == bb.size());
    for (int i = 0; i < aa.size(); i++) {
        ret_mpz += MyMpz(aa[i]) * MyMpz(bb[i]);
    }
    RistScal ret;
    rist_from_mpz(ret, ret_mpz);
    return ret;
}

RistScal inner_prod(const std::vector<int> &aa, const std::vector<long> &bb) {
    // long ret = 0;
    // assert(aa.size() == bb.size());
    // for (int i = 0; i < aa.size(); i++) {
    //     ret += aa[i] * bb[i];
    // }
    // return ret;


    MyMpz ret_mpz;
    assert(aa.size() == bb.size());
    for (int i = 0; i < aa.size(); i++) {
        ret_mpz += MyMpz(aa[i]) * MyMpz(bb[i]);
    }
    RistScal ret;
    rist_from_mpz(ret, ret_mpz);
    return ret;
}

RistScal inner_prod(const std::vector<long> &aa, const std::vector<int> &bb) {
    // long ret = 0;
    // assert(aa.size() == bb.size());
    // for (int i = 0; i < aa.size(); i++) {
    //     ret += aa[i] * bb[i];
    // }
    // return ret;


    MyMpz ret_mpz;
    assert(aa.size() == bb.size());
    for (int i = 0; i < aa.size(); i++) {
        ret_mpz += MyMpz(aa[i]) * MyMpz(bb[i]);
    }
    RistScal ret;
    rist_from_mpz(ret, ret_mpz);
    return ret;
}


double inner_prod(const std::vector<float> &aa, const std::vector<float> &bb) {
    double ret = 0;
    assert(aa.size() == bb.size());
    for (int i = 0; i < aa.size(); i++) {
        ret += static_cast<double>(aa[i]) * static_cast<double>(bb[i]);
    }
    return ret;
}


RistScal inner_prod(const std::vector<int> &aa, const RistScalVec &bb) {
//    MyMpz ret_mpz;
    assert(aa.size() == bb.size());
//    for (int i = 0; i < aa.size(); i++) {
//        ret_mpz += MyMpz(aa[i]) * MyMpz(bb[i]);
//    }
    RistScal ret;
//    rist_from_mpz(ret, ret_mpz);
//    return ret;

    IntCombRistScal ret_ttmath = 0;
    for (int i = 0; i < aa.size(); i++) {
        if (aa[i] >= 0)
            ret_ttmath += multiply_nonnegint(aa[i], bb[i]);
        else
            ret_ttmath -= multiply_nonnegint(-aa[i], bb[i]);
    }
    ttmath_to_rist(ret, ret_ttmath);
    return ret;

}

RistScal inner_prod(const std::vector<long> &aa, const RistScalVec &bb) {
    RistScalVec aa_rist(aa.size());
    for (int i = 0; i < aa.size(); i++) {
        aa_rist[i] = RistScal(aa[i]);
    }
    return inner_prod(aa_rist, bb);
}


RistScal inner_prod(const RistScalVec &aa, const RistScalVec &bb) {
    return sum(aa * bb);
//
//    MyMpz ret_mpz;
//    assert(aa.size() == bb.size());
//    for (int i = 0; i < aa.size(); i++) {
//        ret_mpz += MyMpz(aa[i]) * MyMpz(bb[i]);
//    }
//    RistScal ret;
//    rist_from_mpz(ret, ret_mpz);
//    return ret;
}

// aa: (m + 1)
// V_top: (n)
// V: (m, n)
// V_top stacked on top of V: (m + 1, n)
// computes aa * (V_top stacked on top of V)
RistScalVec vec_mult_mat(const RistScalVec &aa, const RistScalVec &V_top, const std::vector<std::vector<int>> &V) {
    int m = aa.size() - 1;
    int n = V_top.size();
    assert(V.size() == m);
    assert(V[0].size() == n);

    MpzVec ret_mpz(n);

    MpzVec aa_mpz(m + 1);
    mpz_vec_from_rist(aa_mpz, aa);

    std::for_each(std::execution::par_unseq, ret_mpz.begin(), ret_mpz.end(),
                  [&ret_mpz, &aa_mpz, &V_top, &V, &m](auto &&r) {
                      auto j = &r - ret_mpz.data();
                      for (int i = 0; i < m; i++) {
                          r += aa_mpz[i + 1] * MyMpz(V[i][j]);
                      }
                      r += aa_mpz[0] * MyMpz(V_top[j]);
                  });

    RistScalVec ret(n);
    rist_from_mpz_vec(ret, ret_mpz);
    return ret;
}

// aa: (m + 1, k)
// V_top: (n)
// V: (m, n)
// V_top stacked on top of V: (m + 1, n)
// extends aa horizontally to (m + 1, n) by copying the entire matrix (n / k) times, and possibly partially if remainder > 0
// computes the pointwise product, and sums up every (m + 1, k) block
// result: (n / k) if k divides n, (n / k + 1) otherwise
RistScalVec
vec_block_mult_mat(const RistScalMat &aa, const RistScalVec &V_top, const std::vector<std::vector<int>> &V) {
    int m = aa.size() - 1;
    int k = aa[0].size();
    int n = V_top.size();
    assert(V.size() == m);
    assert(V[0].size() == n);

    int ret_size = (n % k == 0) ? n / k : n / k + 1;
    MpzVec ret_mpz(ret_size);

    MpzMat aa_mpz(m + 1, MpzVec(k));
    mpz_mat_from_rist(aa_mpz, aa);

    std::for_each(std::execution::par_unseq, ret_mpz.begin(), ret_mpz.end(),
                  [&ret_mpz, &aa_mpz, &V_top, &V, &m, &k, &n](auto &&r) {
                      auto block_id = &r - ret_mpz.data();
                      for (int i = 0; i < m; i++) {
                          for (int j = k * block_id; j < k * (block_id + 1) && j < n; j++) {
                              r += aa_mpz[i + 1][j - k * block_id] * MyMpz(V[i][j]);

                          }
                      }
                      for (int j = k * block_id; j < k * (block_id + 1) && j < n; j++)
                          r += aa_mpz[0][j - k * block_id] * MyMpz(V_top[j]);
                  });

    RistScalVec ret(ret_size);
    rist_from_mpz_vec(ret, ret_mpz);

    return ret;
}

RistScalVec
mat_mult_vec(const RistScalVec &V_top, const std::vector<std::vector<int>> &V, const std::vector<long> &aa) {
    int n = V_top.size();
    int m = V.size();
    assert(V[0].size() == n);
    assert(aa.size() == n);

    RistScalVec ret(m + 1);
    std::for_each(std::execution::par_unseq, ret.begin(), ret.end(),
                  [&ret, &aa, &V_top, &V, &m, &n](auto &&r) {
                      auto i = &r - ret.data();
                      if (i == 0) {
                          r = inner_prod(aa, V_top);
                      } else {
                          r = RistScal(inner_prod(aa, V[i - 1]));
                      }
                  });
    return ret;
}

RistScal power_of_two(int n) {
    RistScal power = c_scal_zero;
    power.scalar[n / 8] = 1 << (n % 8);
    return power;
}

RistScal ristscal_from_positive_float(float a, int bit_shifter) {
    if (bit_shifter >= 32) {
        long a_shifted_32 = static_cast<long>(static_cast<double>(a) * (1L << 32));
        return RistScal(a_shifted_32) * power_of_two(bit_shifter - 32);
    }
    else {
        long a_shifted = static_cast<long>(static_cast<double>(a) * (1L << bit_shifter));
        return RistScal(a_shifted);
    }
}

RistScal ristscal_from_float(float a, int bit_shifter) {
    if (a >= 0) {
        return ristscal_from_positive_float(a, bit_shifter);
    }
    else {
        return -ristscal_from_positive_float(-a, bit_shifter);
    }
}

long int_from_ristscal(const RistScal &r) {
    long a;
    memcpy(&a, r.scalar.data(), sizeof a);
    return a;
}

float positive_float_from_ristscal(const RistScal &r, int bit_shifter, int bound) {
    if (bound <= 32) {
//        std::cout << "positive float from ristscal small bound <= 32: " << int_from_ristscal(r) << "\n";
//        std::cout << "return val: " << static_cast<double>(int_from_ristscal(r)) / (1L << bit_shifter) << "\n";
        return static_cast<double>(int_from_ristscal(r)) / (1L << bit_shifter);
    }
    long a = 0;
    int byte_count = bound / 8 + 1;
    memcpy(&a, r.scalar.data() + byte_count - 5, 5);

//    std::cout << "positive float from ristscal: " << a << "\n";

    return static_cast<double>(a) / (1L << (bit_shifter - 8 * (byte_count - 5)));
}


float float_from_ristscal(const RistScal &r, int bit_shifter, int max_bit_bound) {
    // assume: bit_shifter < max_bit_bound

    std::cout << "\n\tfloat from ristscal: " << "\n"
        << "\tbit shifter: " << bit_shifter << "\n"
        << "\tmax_bit_bound: " << max_bit_bound << "\n";

    int bound = get_power_two_bound(r);

    std::cout << "\tbound: " << bound << "\n";

    if (bound < max_bit_bound) {
        return positive_float_from_ristscal(r, bit_shifter, bound);
    }
    else {
        return -positive_float_from_ristscal(-r, bit_shifter, bound);
    }
}

IntCombRistScal rist_to_ttmath(const RistScal &r) {
    IntCombRistScal n;
    rist_to_ttmath(n, r);
    return n;
}

void rist_to_ttmath(IntCombRistScal &n, const RistScal &r) {
    n = 0;
    memcpy(n.table, r.scalar.data(), RISTSCALBYTES);
}

void rist_to_ttmath(std::vector<IntCombRistScal> &nn, const RistScalVec &rr) {
    nn.resize(rr.size());
    for (int i = 0; i < rr.size(); i++) {
        rist_to_ttmath(nn[i], rr[i]);
    }
}

void rist_to_ttmath(std::vector<std::vector<IntCombRistScal>> &N, const RistScalMat &R) {
    N.resize(R.size());
    for (int i = 0; i < R.size(); i++) {
        rist_to_ttmath(N[i], R[i]);
    }
}


RistScal ttmath_to_rist(const IntCombRistScal &n) {
    RistScal r;
    ttmath_to_rist(r, n);
    return r;
}

void ttmath_to_rist(RistScal &r, const IntCombRistScal &n) {
    static const IntCombRistScal c_rist_order_ttmath = "7237005577332262213973186563042994240857116359379907606001950938285454250989";
    auto rem = n % c_rist_order_ttmath;
    if (rem < 0) {
        rem += c_rist_order_ttmath;
    }
    memcpy(r.scalar.data(), rem.table, RISTSCALBYTES);
}

void ttmath_to_rist(RistScalVec &rr, const std::vector<IntCombRistScal> &nn) {
    rr.resize(nn.size());
    for (int i = 0; i < nn.size(); i++) {
        ttmath_to_rist(rr[i], nn[i]);
    }
}

// assumption: a >= 0 (not checked)
IntCombRistScal multiply_nonnegint(int a, const RistScal &r) {
    __m256i a_broadcast = _mm256_set1_epi64x(a);
    __m256i r_256 = _mm256_loadu_si256((__m256i *) r.scalar.data());
    __m256i r_hi = _mm256_srli_epi64(r_256, 32);
    __m256i r_lo = _mm256_srli_epi64(_mm256_slli_epi64(r_256, 32), 32);
    __m256i prod_hi = _mm256_mul_epu32(a_broadcast, r_hi);
    __m256i prod_lo = _mm256_mul_epu32(a_broadcast, r_lo);
    IntCombRistScal ret_lo, ret_hi;
    memcpy((unsigned char *) ret_lo.table, &prod_lo, 256 / 8);
    memset((unsigned char *) ret_lo.table + 256 / 8, 0, 64 / 8);
    memset((unsigned char *) ret_hi.table, 0, 32 / 8);
    memcpy((unsigned char *) ret_hi.table + 32 / 8, &prod_hi, 256 / 8);
    memset((unsigned char *) ret_hi.table + (256 + 32) / 8, 0, 32 / 8);
    return ret_lo + ret_hi;
}

std::vector<std::uint32_t> generate_multiple_seeds(const RistHashbytes &seed_source, int count) {
    int hashbytes_in_int_length = seed_source.hashbytes.size() / sizeof(unsigned int); // this value is 16 = 64 / 4.
    std::vector<std::uint32_t> seed_source_uint(hashbytes_in_int_length);
    std::copy(seed_source.hashbytes.begin(), seed_source.hashbytes.end(),
              (unsigned char *) seed_source_uint.data());
    std::seed_seq seeds_converter(seed_source_uint.begin(),
                                  seed_source_uint.end());
    std::vector<std::uint32_t> seeds(count);
    seeds_converter.generate(seeds.begin(), seeds.end());
    return seeds;
}

int get_bit_length(int n) {
    assert(n >= 0);
    int b = 0;
    int cur = n;
    while (cur) {
        b++;
        cur >>= 1;
    }
    return b;
}

RistScalVec left_half_power_two(const RistScalVec &rr) {
    int n = rr.size();
    return RistScalVec(rr.begin(), rr.begin() + n / 2);
}

RistScalVec right_half_power_two(const RistScalVec &rr) {
    int n = rr.size();
    return RistScalVec(rr.begin() + n / 2, rr.end());
}

// p=0: left block
// p=1: right block
RistScalVec block_mult_half(const RistScalVec &aa, const RistScalVec &rr, int p) {
    RistScalVec ret(rr.size());
    for (int i = 0; i < rr.size(); i++) {
        if (((i / aa.size()) & 1) == p) {
            ret[i] = aa[i % aa.size()] * rr[i];
        } else {
            ret[i] = c_scal_zero;
        }
    }
    return ret;
}

RistScalVec block_mult(const RistScal &a, const RistScal &b, const RistScalVec &rr, int half_step) {
    RistScalVec ret(rr.size());
    for (int i = 0; i < rr.size(); i++) {
        if (((i / half_step) & 1) == 0) {
            ret[i] = a * rr[i];
        } else {
            ret[i] = b * rr[i];
        }
    }
    return ret;
}

RistScal power(int n, const RistScal &r) {
    RistScal ret = c_scal_one;
    for (int i = 0; i < n; i++) {
        ret *= r;
    }
    return ret;
}