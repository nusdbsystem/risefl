//
// Created by yizheng on 15/3/23.
//

#include <vector>
#include <cstring>
#include <algorithm>
#include <execution>
#include <cmath>
#include <random>
#include "include/risefl_interface_common.h"
#include "include/zkp_hash.h"
#include "include/ristretto_vector.h"
#include "include/mpz_conversion.h"

std::vector<unsigned char> bytes_from_pub_vec(const std::vector<PublicKey> &p) {
    std::vector<unsigned char> b;
    b.reserve(p.size());
    for (auto &&pub: p) {
        b.insert(b.end(), pub.pub.begin(), pub.pub.end());
    }
    return b;
}

std::vector<unsigned char> int_vec_to_bytes(const std::vector<int> &p) {
    std::vector<unsigned char> b;
    b.reserve(sizeof(int) * p.size());
    for (auto &&i: p) {
        auto temp = bytes_from_int(i);
        b.insert(b.end(), temp.begin(), temp.end());
    }
    return b;
}

int round_to_integer(double d, int bit_length) {
    return std::lround(d * (1 << bit_length));
}

int extra_inner_prod(CHECK_TYPE c) {
    if (c == CHECK_TYPE::L2NORM || c == CHECK_TYPE::SPHERE) {
        return 0;
    }
    else {
        return 1;
    }
}

int extra_square(CHECK_TYPE c) {
    if (c == CHECK_TYPE::COSINE_SIM) {
        return 1;
    }
    else {
        return 0;
    }
}


CheckParam
CommonInterface::set_normalizing_factor_and_compute_check_param_from_gamma(CheckParamFloat check_param_float) {
    float gamma = gammas_from_num_samples[predicate.num_samples - 1];
    double temp = sqrt(gamma) + sqrt(predicate.num_samples * dim) / (1 << (predicate.random_normal_bit_shifter + 1));
    float gamma_calibrated = temp * temp;
    int weight_bit_shifter = weight_bits - 1;
    int bit_shifter = 2 * (weight_bits - 1 + predicate.random_normal_bit_shifter);

    auto check_type = check_param_float.check_type;
    CheckParam ret(check_type);
    if (check_type == CHECK_TYPE::L2NORM) {
        normalizing_factor = check_param_float.l2_param.bound;
        ret.l2_param.bound_sq = ristscal_from_positive_float(gamma_calibrated, bit_shifter);
    }
    if (check_type == CHECK_TYPE::SPHERE) {
        double center_norm = sqrt(inner_prod(check_param_float.sphere_param.center, check_param_float.sphere_param.center));
        normalizing_factor = center_norm + check_param_float.sphere_param.bound;
        ret.sphere_param.bound_sq = ristscal_from_positive_float(gamma_calibrated / (normalizing_factor * normalizing_factor), bit_shifter);
        ret.sphere_param.center.resize(dim);
        for (int i = 0; i < dim; i++) {
            ret.sphere_param.center[i] = ristscal_from_float(check_param_float.sphere_param.center[i] / normalizing_factor, weight_bit_shifter);
        }
    }
    if (check_type == CHECK_TYPE::COSINE_SIM) {
        auto pivot = check_param_float.cosine_param.pivot;
        double pivot_norm = sqrt(inner_prod(pivot, pivot));
        normalizing_factor = check_param_float.cosine_param.bound;
        ret.cosine_param.bound_sq = ristscal_from_positive_float(gamma_calibrated, bit_shifter);
        ret.cosine_param.pivot.resize(dim);
        for (int i = 0; i < dim; i++) {
            ret.cosine_param.pivot[i] = ristscal_from_float(pivot[i] / pivot_norm, weight_bit_shifter);
        }

        auto cosine_bound = check_param_float.cosine_param.cosine_bound;
        ret.cosine_param.l2_sq_multiplier = c_scal_one;
        ret.cosine_param.inner_prod_sq_multiplier = ristscal_from_float(gamma_calibrated / (cosine_bound * cosine_bound),   2 * predicate.random_normal_bit_shifter);
    }
    return ret;
}


void Predicate::initialize_from_seed(const std::vector<unsigned char> &seed) {
    std::for_each(std::execution::par_unseq, hh.begin(), hh.end(),
                  [this, &seed](auto &&h) {
                      auto i = &h - hh.data();
                      h = p3_from_hash_from_bytes(seed, bytes_from_int(0), bytes_from_int(i));
                      if (b_precomp)
                          generate_precomp_table(hh_precomp[i], h);
                  });
//    for (int i = 0; i < hh.size(); i++) {
//        hh[i] = p3_from_hash_from_bytes(seed, bytes_from_int(i));
//        generate_precomp_table(hh_precomp[i], hh[i]);
//    }
    for (int i = 0; i < bound_elem_keys_1.size(); i++) {
        bound_elem_keys_1.elems[i] = p3_from_hash_from_bytes(seed, bytes_from_int(1), bytes_from_int(i));
        bound_elem_keys_2.elems[i] = p3_from_hash_from_bytes(seed, bytes_from_int(2), bytes_from_int(i));
    }
    bound_elem_keys_1.fill_bytes();
    bound_elem_keys_2.fill_bytes();

    square_key.elem = p3_from_hash_from_bytes(seed, bytes_from_int(3));
    square_key.fill_bytes();
}


void import_flags_start_from_1_from_bytestream(std::vector<int> &flags,
                                               std::vector<unsigned char>::const_iterator &it) {
    for (int i = 1; i < flags.size(); i++) {
        flags[i] = *it;
        ++it;
    }
}

void export_flags_start_from_1_to_bytestream(const std::vector<int> &flags, std::vector<unsigned char>::iterator &it) {
    for (int i = 1; i < flags.size(); i++) {
        *it = flags[i];
        ++it;
    }
}
