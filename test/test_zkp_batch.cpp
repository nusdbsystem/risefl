//
// Created by yizheng on 13/4/23.
//
#include <iostream>
#include <chrono>
#include <cassert>
#include <vector>
#include <string>
#include <algorithm>
#include <execution>
#include <thread>
#include <sodium.h>
#include "../include/ristretto.h"
#include "../include/ristretto_vector.h"
#include "../include/zkp.h"

//void test_zkp_batch() {
//    RistElem h_orig, f_orig;
//    rand_init(h_orig);
//    rand_init(f_orig);
//
//    RistP3MatAndBytes hh(1, 1);
//    p3_from_bytes(hh.elems[0][0], h_orig.element);
//    hh.fill_bytes();
//
//    RistP3AndBytes f;
//    p3_from_bytes(f.elem, f_orig.element);
//    f.fill_bytes();
//
//    RistP3VecAndBytes yy(1), zz(1), yy_prime(1);
//    RistScalVec xx(1), qq(1), rr(1);
//    rand_init(xx);
//    rand_init(qq);
//    rand_init(rr);
//
//    pedersen_commit(yy.elems[0], xx[0], hh.elems[0][0], rr[0]);
//    pedersen_commit(yy_prime.elems[0], xx[0], f.elem, qq[0]);
//    pedersen_zero_commit(zz.elems[0], rr[0]);
//
//    auto proof = generate_pedersen_batch_with_zero_proof(hh, zz, yy,
//                                                         f,
//                                                         yy_prime,
//                                                         rr, xx, qq);
//
//    assert(b_verify_pedersen_batch_with_zero(hh, zz, yy,
//                                             f,
//                                             yy_prime,
//                                             proof));
//
//    RistP3AndBytes h;
//    h.elem = hh.elems[0][0];
//    h.fill_bytes();
//
//    RistP3VecAndBytes bound_keys_1(1), bound_keys_2(1);
//    bound_keys_1.elems[0] = f.elem;
//    bound_keys_1.fill_bytes();
//
//    RistElem bound_key_2_orig;
//    rand_init(bound_key_2_orig);
//    p3_from_bytes(bound_keys_2.elems[0], bound_key_2_orig.element);
//    bound_keys_2.fill_bytes();
//
//    RistScal gamma;
//    rand_init(gamma);
//
//    RistScal v(0);
//    RistP3AndBytes V;
//    pedersen_commit(V.elem, v, h.elem, gamma);
//    V.fill_bytes();
//
//    auto proof_bound = generate_range_proof_power_two(h, V, bound_keys_1, bound_keys_2,
//                                                      1, gamma, v);
//
//    assert(b_verify_range_power_two(h, V, bound_keys_1, bound_keys_2,
//                                    1, proof_bound));
//
//    std::cout << "batch zkp passed!\n";
//
//}
