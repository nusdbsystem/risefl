//
// Created by yizheng on 4/4/23.
//
#include <vector>
#include <random>
#include <iostream>
#include <algorithm>
#include <execution>
#include "simulate.h"
#include <random>

void bench_random(int k, int d) {
    static const double twopi = 2.0 * 3.14159265358979323846;

    // generate seed for each a_i
//    int hashbytes_in_int_length = seed.hashbytes.size() / sizeof(unsigned int); // this value is 16 = 64 / 4.
//    std::vector<std::uint32_t> seed_of_a_i_source_uint(hashbytes_in_int_length);
//    std::copy(seed.hashbytes.begin(), seed.hashbytes.end(),
//              (unsigned char *) seed_of_a_i_source_uint.data());
//    std::seed_seq seed_at_a_i_converter(seed_of_a_i_source_uint.begin(),
//                                        seed_of_a_i_source_uint.end());
//    std::vector<std::uint32_t> seeds_at_a_m(num_samples + 1);
//    seed_at_a_i_converter.generate(seeds_at_a_m.begin(), seeds_at_a_m.end());

    std::mt19937_64 generator;
    generator.seed(111);


//    for (auto &&a: aa_0) {
//        std::uniform_int_distribution<unsigned long> unif_unsigned_long;
//        RistHashbytes hash;
//        for (int j = 0; j < hash.hashbytes.size(); j += sizeof(unsigned long)) {
//            auto u = unif_unsigned_long(generator);
//            memcpy(hash.hashbytes.data() + j, &u, sizeof(unsigned long));
//        }
//        scalar_reduce(a, hash);
//    }

//    int k = 100;
//    int d = 1e6;
    std::vector<std::vector<int>> aa_pos(k, std::vector<int>(d));

    std::for_each(
            std::execution::par_unseq, // or std::execution::par_unseq if no conflicts
            aa_pos.begin(),
            aa_pos.end(),
            [&](auto &&a_m) {
                int m = &a_m - aa_pos.data();

                std::mt19937_64 generator;
                generator.seed(111 + m);

                // Box-Muller transform
                // https://github.com/miloyip/normaldist-benchmark/blob/master/src/boxmuller.cpp
                for (int i = 0; i < d; i += 2) {
                    double u1 = double(1.0) - std::generate_canonical<double, 64>(generator);
                    double u2 = std::generate_canonical<double, 64>(generator);
                    double radius = std::sqrt(-2 * std::log(u1));
                    double theta = twopi * u2;
                    double r1 = radius * std::cos(theta);
                    double r2 = radius * std::cos(theta);
                    aa_pos[m][i] = std::lround(r1 * (1 << 16));
                    if (i + 1 < d)
                        aa_pos[m][i + 1] = std::lround(r2 * (1 << 16));
                }
            });

}

void bench_random() {
    int k = 100;
    int d = 1e6;
    std::cout << "k = " << k
    << "\n d = " << d
    << "\n " << measure_time([&](){ bench_random(k, d); }) << "\n";
}