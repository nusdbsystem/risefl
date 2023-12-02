//
// Created by yizheng on 15/3/23.
//

#include <random>
#include "include/di_zkp_interface_client.h"
#include "include/di_zkp_interface_common.h"
#include "include/shamir.h"
#include "include/rist_fast_computation.h"
#include "include/random_avx/mersenne-twister-avx2.h"


void ClientInterface::generate_batch_shares_and_check_string() {
    batch_shamir_share_with_check = compute_batch_shamir_shares_with_check_string(blinds_to_share, num_clients,
                                                                                  max_malicious_clients + 1);
}

void ClientInterface::encrypt_shares() {
    std::for_each(std::execution::par_unseq, encrypted_shamir_shares.begin() + 1, encrypted_shamir_shares.end(),
                  [this](auto &&encrypted) {
                      auto j = &encrypted - encrypted_shamir_shares.data();
                      if (j != client_id && !server_flags[j]) {
                          for (int i = 0; i < predicate.num_blinds_per_group_element; i++) {
                              encrypted_shamir_shares[j][i] = encrypt(
                                      std::vector<unsigned char>(
                                              batch_shamir_share_with_check.shares[i][j].scalar.begin(),
                                              batch_shamir_share_with_check.shares[i][j].scalar.end()),
                                      dh_public_key_collection[j], dh_private_key);
                          }
                      }
                  });
}

void ClientInterface::decrypt_shamir_shares() {
    std::for_each(std::execution::par_unseq, other_shamir_shares.begin() + 1, other_shamir_shares.end(),
                  [this](auto &&shares) {
                      int j = &shares - other_shamir_shares.data();
                      if (j != client_id && !server_flags[j]) {
                          try {
                              for (int i = 0; i < predicate.num_blinds_per_group_element; i++) {
                                  auto decrpyted = decrypt(other_encrypted_shamir_shares[j][i],
                                                           dh_public_key_collection[j], dh_private_key);
                                  std::copy(decrpyted.begin(), decrpyted.end(),
                                            other_shamir_shares[j][i].scalar.data());
                              }
                          }
                          catch (DecryptionError &) {
                              flags[j] = 1;
                          }
                      }
                  });
}

void ClientInterface::check_shamir_share_integrity() {
    std::for_each(std::execution::par_unseq, other_shamir_shares.begin() + 1, other_shamir_shares.end(),
                  [this](auto &&shares) {
                      int j = &shares - other_shamir_shares.data();
                      if (j != client_id && !server_flags[j]) {
                          if (!b_shamir_verify(other_shamir_shares[j],
                                               shamir_check_string_collection[j].elems,
                                               max_malicious_clients + 1,
                                               client_id)) {
                              flags[j] = 1;
                          }
                      }
                  });
}

void ClientInterface::update_other_shamir_shares_with_dispute() {
    for (int j = 1; j <= num_clients; j++) {
        if (!server_flags[j] && flags[j]) {
            other_shamir_shares[j] = other_dispute_shares[j];
        }
    }
}

void ClientInterface::compute_aggegrated_share() {
    for (int k = 0; k < predicate.num_blinds_per_group_element; k++) {
        aggregates[k] = batch_shamir_share_with_check.shares[k][client_id];
        for (int j = 1; j <= num_clients; j++) {
            if (j != client_id && !server_flags[j]) {
                aggregates[k] = aggregates[k] + other_shamir_shares[j][k];
            }
        }
    }
}


void ClientInterface::generate_dispute_shares() {
    int count = 0;
    for (int j = 1; j <= num_clients; j++) {
        if (dispute_clients[j]) {
            for (int k = 0; k < predicate.num_blinds_per_group_element; k++) {
                dispute_shares[j][k] = batch_shamir_share_with_check.shares[k][j];
            }
            count++;
        }
    }
    if (count > max_malicious_clients)
        throw Abort();
}

void ClientInterface::generate_proof() {
    auto check_type = predicate.check_param.check_type;
    int num_samples_extra = predicate.num_samples + extra_inner_prod(check_type);
    int num_sq = predicate.num_samples + extra_square(check_type);

    RistScalVec inner_prods(num_samples_extra + 1);
//    = mat_mult_vec(predicate.aa_0, predicate.aa_pos, weight_updates);

    // compute inner prods and check correctness of hh_comb

    RistScalMat cc(num_samples_extra + 1, RistScalVec(predicate.num_blinds_per_group_element));
    for (auto &&c: cc)
        rand_init(c);

    RistScalVec aa_cc(predicate.num_weight_keys, c_scal_zero);

    LinearCombCalculator lcc;
    lcc.add(cc, predicate.hh_comb.elems);

    auto seeds_at_a_m = generate_multiple_seeds(predicate.aa_seed, predicate.num_samples + 1);

    RistScalVec top_row(dim);
    for (int l = 0; l < dim; l++) {
        top_row[l] = rist_scalar_from_hash_from_bytes(predicate.aa_seed.hashbytes, bytes_from_int(l));
        aa_cc[l / predicate.num_blinds_per_group_element] +=
                top_row[l] * cc[0][l % predicate.num_blinds_per_group_element];
        if (check_type == CHECK_TYPE::COSINE_SIM) {
            aa_cc[l / predicate.num_blinds_per_group_element] +=
                    predicate.check_param.cosine_param.pivot[l] * cc[cc.size() - 1][l % predicate.num_blinds_per_group_element];
        }
    }

    if (check_type == CHECK_TYPE::L2NORM || check_type == CHECK_TYPE::COSINE_SIM || check_type == CHECK_TYPE::ZENO) {
        inner_prods[0] = inner_prod(weight_updates, top_row);
    }
    RistScalVec weight_updates_rist = rist_scal_vec_from_long_vec(weight_updates);
    if (check_type == CHECK_TYPE::SPHERE) {
        inner_prods[0] = inner_prod(weight_updates_rist - predicate.check_param.sphere_param.center, top_row);
    }
    if (check_type == CHECK_TYPE::COSINE_SIM) {
        inner_prods[inner_prods.size() - 1] = inner_prod(weight_updates, predicate.check_param.cosine_param.pivot);
    }

    std::vector<std::vector<int>> temp_A(NUM_THREADS, std::vector<int>(dim));

    std::vector<IntCombRistScal> aa_cc_ttmath(predicate.num_weight_keys);
    rist_to_ttmath(aa_cc_ttmath, aa_cc);

    for (int m = 1; m <= predicate.num_samples; m += NUM_THREADS) {
        // TODO: non-AVX2
        assert(__AVX2__);

        std::for_each(
                std::execution::par_unseq,
                temp_A.begin(),
                temp_A.end(),
                [this, &inner_prods, &m, &temp_A, &seeds_at_a_m, &check_type, &weight_updates_rist](auto &&temp_A_row) {
                    auto j = &temp_A_row - temp_A.data();
                    if (m + j <= predicate.num_samples) {
                        mt_avx2 rand_gen(seeds_at_a_m[m + j]);
                        rand_gen.rand_normal_avx2_shifted(temp_A[j].data(), dim, predicate.random_normal_bit_shifter);
                        if (check_type == CHECK_TYPE::L2NORM || check_type == CHECK_TYPE::COSINE_SIM || check_type == CHECK_TYPE::ZENO) {
                            inner_prods[m + j] = RistScal(inner_prod(temp_A[j], weight_updates));
                        }
                        if (check_type == CHECK_TYPE::SPHERE) {
                            inner_prods[m + j] = RistScal(inner_prod(temp_A[j], weight_updates_rist - predicate.check_param.sphere_param.center));
                        }
                    }
                }
        );

        std::for_each(std::execution::par_unseq, aa_cc_ttmath.begin(), aa_cc_ttmath.end(),
                      [this, &m, &aa_cc_ttmath, &cc, &temp_A](auto &&aa_cc_ttmath_l) {
                          auto l = &aa_cc_ttmath_l - aa_cc_ttmath.data();
                          for (int j = 0; (j < NUM_THREADS) && (m + j <= predicate.num_samples); j++) {
                              for (int p = 0; (p < predicate.num_blinds_per_group_element) &&
                                              (l * predicate.num_blinds_per_group_element + p < dim); p++) {
                                  if (temp_A[j][l * predicate.num_blinds_per_group_element + p] >= 0)
                                      aa_cc_ttmath_l += multiply_nonnegint(
                                              temp_A[j][l * predicate.num_blinds_per_group_element + p], cc[m + j][p]);
                                  else
                                      aa_cc_ttmath_l -= multiply_nonnegint(
                                              -temp_A[j][l * predicate.num_blinds_per_group_element + p], cc[m + j][p]);
                              }
                          }
                      });
    }

    ttmath_to_rist(aa_cc, aa_cc_ttmath);

    lcc.sub(aa_cc, predicate.hh);
    if (!lcc.is_0()) {
        throw Abort();
    }

    proof.linear_comb_batch_commitments.resize(num_samples_extra + 1);
    for (int i = 0; i < proof.linear_comb_batch_commitments.size(); i++) {
        pedersen_zero_commit(proof.linear_comb_batch_commitments.elems[i], inner_prods[i]);
        proof.linear_comb_batch_commitments.elems[i] += linear_comb(blinds_to_share, predicate.hh_comb.elems[i]);
    }
    proof.linear_comb_batch_commitments.fill_bytes();

    proof.linear_comb_single_commitments.resize(num_samples_extra + 1);
    RistScalVec blinds_on_linear_comb_single(num_samples_extra + 1);
    rand_init(blinds_on_linear_comb_single);
    for (int i = 0; i < proof.linear_comb_batch_commitments.size(); i++) {
        pedersen_commit(proof.linear_comb_single_commitments.elems[i], inner_prods[i], predicate.square_key.elem,
                        blinds_on_linear_comb_single[i]);;
    }
    proof.linear_comb_single_commitments.fill_bytes();

    proof.proof_well_formed = generate_pedersen_batch_with_zero_proof(predicate.hh_comb,
                                                                      batch_shamir_share_with_check.check_strings.slice_col_0(),
                                                                      proof.linear_comb_batch_commitments,
                                                                      predicate.square_key,
                                                                      proof.linear_comb_single_commitments,
                                                                      blinds_to_share, inner_prods,
                                                                      blinds_on_linear_comb_single);

    RistP3VecAndBytes linear_comb_single_commitments_samples;
    RistScalVec blinds_on_linear_comb_single_samples;
    RistScalVec inner_prods_nonzero = no_head(inner_prods);
    RistScalVec inner_prods_samples;
    if (check_type == CHECK_TYPE::L2NORM || check_type == CHECK_TYPE::SPHERE) {
        linear_comb_single_commitments_samples = proof.linear_comb_single_commitments.no_head();
        blinds_on_linear_comb_single_samples = no_head(blinds_on_linear_comb_single);
        inner_prods_samples = inner_prods_nonzero;
    }
    else {
        linear_comb_single_commitments_samples = proof.linear_comb_single_commitments.no_head_no_tail();
        blinds_on_linear_comb_single_samples = no_head_no_tail(blinds_on_linear_comb_single);
        inner_prods_samples = no_tail(inner_prods_nonzero);
    }

    proof.proof_linear_comb_bound = generate_abs_val_range_proof_power_two_agg(predicate.square_key,
                                                                               linear_comb_single_commitments_samples,
                                                                               predicate.bound_elem_keys_1,
                                                                               predicate.bound_elem_keys_2,
                                                                               predicate.inner_prod_bound_bits,
                                                                               predicate.num_samples,
                                                                               blinds_on_linear_comb_single_samples,
                                                                               inner_prods_samples);

    RistScalVec blinds_on_square(num_sq);
    rand_init(blinds_on_square);
    RistP3VecAndBytes linear_comb_single_commitments_to_sq = proof.linear_comb_single_commitments.no_head();
    RistScalVec inner_prods_to_sq = inner_prods_nonzero;
    auto blinds_on_single_to_sq = no_head(blinds_on_linear_comb_single);

    if (check_type == CHECK_TYPE::ZENO) {
        inner_prods_to_sq = no_tail(inner_prods_to_sq);
        linear_comb_single_commitments_to_sq = proof.linear_comb_single_commitments.no_head_no_tail();
        blinds_on_single_to_sq = no_tail(blinds_on_single_to_sq);
    }
    proof.square_commitments.resize(num_sq);
    pedersen_commit(proof.square_commitments.elems, inner_prods_to_sq * inner_prods_to_sq,
                    predicate.square_key.elem,
                    blinds_on_square);
    proof.square_commitments.fill_bytes();

    proof.proof_squares = generate_batch_pedersen_with_square_proof(predicate.square_key,
                                                                    linear_comb_single_commitments_to_sq,
                                                                    proof.square_commitments,
                                                                    inner_prods_to_sq,
                                                                    blinds_on_single_to_sq,
                                                                    blinds_on_square);

    RistP3AndBytes bound_minus_sum_of_square_commitments;
    if (check_type == CHECK_TYPE::L2NORM) {
        auto diff = predicate.check_param.l2_param.bound_sq - inner_prod(inner_prods_nonzero, inner_prods_nonzero);
        bound_minus_sum_of_square_commitments.elem = pedersen_commit(diff,
                                                                     predicate.square_key.elem,
                                                                     -sum(blinds_on_square));
        bound_minus_sum_of_square_commitments.fill_bytes();
        proof.proof_sum_range = generate_range_proof_power_two(predicate.square_key,
                                                               bound_minus_sum_of_square_commitments,
                                                               predicate.bound_elem_keys_1,
                                                               predicate.bound_elem_keys_2,
                                                               predicate.max_bound_sq_bits,
                                                               -sum(blinds_on_square),
                                                               diff);
    }
    if (check_type == CHECK_TYPE::SPHERE) {
        auto diff = predicate.check_param.sphere_param.bound_sq - inner_prod(inner_prods_nonzero, inner_prods_nonzero);
        bound_minus_sum_of_square_commitments.elem = pedersen_commit(diff,
                                                                     predicate.square_key.elem,
                                                                     -sum(blinds_on_square));
        bound_minus_sum_of_square_commitments.fill_bytes();
        proof.proof_sum_range = generate_range_proof_power_two(predicate.square_key,
                                                               bound_minus_sum_of_square_commitments,
                                                               predicate.bound_elem_keys_1,
                                                               predicate.bound_elem_keys_2,
                                                               predicate.max_bound_sq_bits,
                                                               -sum(blinds_on_square),
                                                               diff);
    }
    if (check_type == CHECK_TYPE::COSINE_SIM) {
        auto diff = predicate.check_param.cosine_param.bound_sq - inner_prod(inner_prods_samples, inner_prods_samples);
        bound_minus_sum_of_square_commitments.elem = pedersen_commit(diff,
                                                                     predicate.square_key.elem,
                                                                     -sum(no_tail(blinds_on_square)));
        bound_minus_sum_of_square_commitments.fill_bytes();
        proof.proof_sum_range = generate_range_proof_power_two(predicate.square_key,
                                                               bound_minus_sum_of_square_commitments,
                                                               predicate.bound_elem_keys_1,
                                                               predicate.bound_elem_keys_2,
                                                               predicate.max_bound_sq_bits,
                                                               -sum(no_tail(blinds_on_square)),
                                                               diff);

        auto inner_prod_sq_multiplier = predicate.check_param.cosine_param.inner_prod_sq_multiplier;
        auto l2_sq_multiplier = predicate.check_param.cosine_param.l2_sq_multiplier;
        auto inner_prod_pivot = inner_prods[inner_prods.size() - 1];
        diff =  inner_prod_sq_multiplier * inner_prod_pivot * inner_prod_pivot
                - l2_sq_multiplier *
                    inner_prod(inner_prods_samples, inner_prods_samples) *
                    inner_prod(predicate.check_param.cosine_param.pivot, predicate.check_param.cosine_param.pivot);
        RistP3AndBytes cosine_diff_commitment;
        auto cosine_blind = inner_prod_sq_multiplier * blinds_on_square[blinds_on_square.size() - 1] -
                l2_sq_multiplier * sum(no_tail(blinds_on_square)) * inner_prod(predicate.check_param.cosine_param.pivot, predicate.check_param.cosine_param.pivot);
        cosine_diff_commitment.elem = pedersen_commit(diff, predicate.square_key.elem, cosine_blind);
        cosine_diff_commitment.fill_bytes();
        proof.proof_sum_range_for_cosine = generate_range_proof_power_two(predicate.square_key,
                                                                          cosine_diff_commitment,
                                                                          predicate.bound_elem_keys_1,
                                                                          predicate.bound_elem_keys_2,
                                                                          predicate.max_bound_sq_bits,
                                                                          cosine_blind,
                                                                          diff);
    }
}

//std::vector<unsigned char> ClientInterface::generate_encrypted_shares_and_check_string() {
//    generate_batch_shares_and_check_string();
//    encrypt_shares();
//    std::vector<unsigned char> out;
//    for (auto v : encrypted_shamir_shares) {
//        auto bytes = bytestream_from_cipher_with_nonce(<#initializer#>, 0, v);
//        out.insert(out.end(), bytes.begin(), bytes.end());
//    }
//    auto bytes = bytestream_from_p3_vec(shamir_check_string);
//    out.insert(out.end(), bytes.begin(), bytes.end());
//    return out;
//}

std::vector<unsigned char> ClientInterface::send_1_internal(const CheckParam &check_param, const std::vector<long> &u) {
    if (protocol_type == PROTOCOL_TYPE::NON_PRIV_INT) {
        std::vector<unsigned char> ret(sizeof(long) * dim);
        for (int i = 0; i < u.size(); i++) {
            memcpy(ret.data() + sizeof(long) * i, &u[i], sizeof(long));
        }
        return ret;
    }
    import_weight_updates(u);
    reset_server_flags();
    reset_flags();
    set_param(check_param);
//    set_l2_bound_sq(B);
    generate_dh_key_pair();

    auto signed_public_key = sign(std::vector<unsigned char>(dh_public_key.pub.begin(), dh_public_key.pub.end()),
                                  bul_prv_key);
    assert(signed_public_key.size() == SIGNBYTES + crypto_box_PUBLICKEYBYTES);

    std::vector<unsigned char> ret(signed_public_key.size());
    std::copy(signed_public_key.begin(), signed_public_key.end(), ret.begin());
    return ret;
}

std::vector<unsigned char> ClientInterface::send_1_bytes(const CheckParamFloat &check_param_float,
                                                         const std::vector<float> &u_float) {
    if (protocol_type == PROTOCOL_TYPE::NON_PRIV_FLOAT) {
        std::vector<unsigned char> ret(sizeof(float) * dim);
        for (int i = 0; i < u_float.size(); i++) {
            float temp = u_float[i];
            memcpy(ret.data() + sizeof(float) * i, &temp, sizeof(float));
        }
        return ret;
    }
    auto check_param = set_normalizing_factor_and_compute_check_param_from_gamma(check_param_float);
    assert(u_float.size() == dim);
    std::vector<long> u_int(dim);
    for (int i = 0; i < dim; i++) {
        u_int[i] = static_cast<long>(static_cast<double>(u_float[i] / normalizing_factor) * (1L << (weight_bits - 1)));
    }
    return send_1_internal(check_param, u_int);
}


std::string ClientInterface::send_1(const CheckParamFloat &check_param_float,
                                    const std::vector<float> &u_float) {
    return base64_encode(send_1_bytes(check_param_float, u_float));
}

std::vector<unsigned char> ClientInterface::receive_and_send_2_bytes(const std::vector<unsigned char> &bytes) {
    if (bytes.size() != num_clients * (SIGNBYTES + crypto_box_PUBLICKEYBYTES))
        throw DeserializationError();

    auto it_bytes = bytes.cbegin();
    for (int i = 1; i <= num_clients; i++) {
        std::copy_n(it_bytes, SIGNBYTES + crypto_box_PUBLICKEYBYTES, bul_signed_pub_key_collection[i].begin());
        it_bytes += SIGNBYTES + crypto_box_PUBLICKEYBYTES;
//        dh_public_key_collection[i].import_from_bytestream(it_bytes);
    }
    assert(it_bytes == bytes.end());

    for (int i = 1; i <= num_clients; i++) {
        try {
            auto pub_key = sign_open(bul_signed_pub_key_collection[i], bul_pub_keys[i]);
            std::copy(pub_key.begin(), pub_key.end(), dh_public_key_collection[i].pub.begin());
        }
        catch (InvalidSign &) {
            Abort();
        }
    }

    int size_committed_updates = dim * RISTBYTES;
    int size_encrypted_shares = num_clients * predicate.num_blinds_per_group_element * ENCRISTSCALBYTES;
    int size_check_strings = predicate.num_blinds_per_group_element * (max_malicious_clients + 1) * RISTBYTES;

    std::vector<unsigned char> ret(size_committed_updates + size_encrypted_shares + size_check_strings);

    auto it_ret = ret.begin();

    rand_init(blinds_to_share);
    if (predicate.b_precomp)
        pedersen_commit_p3_to_bytes_from_precomp(it_ret, weight_updates, weight_bits, predicate.hh_precomp,
                                                 blinds_to_share);
    else {
        assert(blinds_to_share.size() == 1);
        pedersen_commit_p3_to_bytes(it_ret,
                                    weight_updates,
                                    weight_bits,
                                    predicate.hh, blinds_to_share[0]);
    }
    generate_batch_shares_and_check_string();

    encrypt_shares();
    for (int j = 1; j <= num_clients; j++) {
        for (int k = 0; k < predicate.num_blinds_per_group_element; k++)
            encrypted_shamir_shares[j][k].export_to_bytestream(it_ret);
    }

    batch_shamir_share_with_check.check_strings.export_to_bytestream(it_ret);
    assert(it_ret == ret.end());
    return ret;
}

std::string ClientInterface::receive_and_send_2(const std::string &bytes_str) {
    return base64_encode(receive_and_send_2_bytes(base64_decode(bytes_str)));
}

std::vector<unsigned char> ClientInterface::receive_and_send_3_bytes(const std::vector<unsigned char> &bytes) {
    int size_encrypted_shares = num_clients * predicate.num_blinds_per_group_element * ENCRISTSCALBYTES;
    int size_check_strings =
            num_clients * predicate.num_blinds_per_group_element * (max_malicious_clients + 1) * RISTBYTES;
    int size_aa_seed = crypto_core_ed25519_HASHBYTES;
    int size_hh_comb = (predicate.num_samples + 1) * predicate.num_blinds_per_group_element * RISTBYTES;

    int size = size_encrypted_shares + size_check_strings + size_aa_seed + size_hh_comb;
    if (bytes.size() != size)
        throw DeserializationError();

    auto it_bytes = bytes.cbegin();
    for (int j = 1; j <= num_clients; j++) {
        for (int k = 0; k < predicate.num_blinds_per_group_element; k++)
            other_encrypted_shamir_shares[j][k].import_from_bytestream(it_bytes);
    }

    for (int j = 1; j <= num_clients; j++) {
        shamir_check_string_collection[j].import_from_bytestream(it_bytes);
    };

    predicate.aa_seed.import_from_bytestream(it_bytes);
    predicate.hh_comb.import_from_bytestream(it_bytes);

    assert(it_bytes == bytes.end());

    generate_proof();

    decrypt_shamir_shares();
    check_shamir_share_integrity();

    std::vector<unsigned char> ret(Proof::size(predicate.check_param.check_type,
                                               predicate.num_blinds_per_group_element,
                                               predicate.num_samples,
                                               predicate.inner_prod_bound_bits,
                                               predicate.max_bound_sq_bits)
                                   + num_clients);

    auto it_ret = ret.begin();
    proof.export_to_bytestream(it_ret);

    export_flags_start_from_1_to_bytestream(flags, it_ret);

    assert(it_ret == ret.end());
    return ret;
}

std::string ClientInterface::receive_and_send_3(const std::string &bytes_str) {
    return base64_encode(receive_and_send_3_bytes(base64_decode(bytes_str)));
}

std::vector<unsigned char> ClientInterface::receive_and_send_4_bytes(const std::vector<unsigned char> &bytes) {
    if (bytes.size() != num_clients)
        throw DeserializationError();

    auto it_bytes = bytes.cbegin();

    import_flags_start_from_1_from_bytestream(dispute_clients, it_bytes);
    assert(it_bytes == bytes.end());

    generate_dispute_shares();

    std::vector<unsigned char> ret(num_clients *predicate.num_blinds_per_group_element * RISTSCALBYTES);

    auto it_ret = ret.begin();
    for (int j = 1; j <= num_clients; j++) {
        if (dispute_clients[j]) {
            export_to_bytestream(dispute_shares[j], it_ret);
        } else {
            export_to_bytestream(RistScalVec(predicate.num_blinds_per_group_element, c_scal_zero), it_ret);
        }
    }
    assert(it_ret == ret.end());
    return ret;
}

std::string ClientInterface::receive_and_send_4(const std::string &bytes_str) {
    return base64_encode(receive_and_send_4_bytes(base64_decode(bytes_str)));
}

std::vector<unsigned char> ClientInterface::receive_and_send_5_bytes(const std::vector<unsigned char> &bytes) {
    if (bytes.size() != num_clients * predicate.num_blinds_per_group_element * RISTSCALBYTES + num_clients)
        throw DeserializationError();
    auto it_bytes = bytes.cbegin();
    for (int j = 1; j <= num_clients; j++) {
        import_from_bytestream(other_dispute_shares[j], it_bytes);
    }
    import_flags_start_from_1_from_bytestream(server_flags, it_bytes);
    assert(it_bytes == bytes.end());

    update_other_shamir_shares_with_dispute();
    compute_aggegrated_share();

    std::vector<unsigned char> ret(predicate.num_blinds_per_group_element * RISTSCALBYTES);
    auto it_ret = ret.begin();
    export_to_bytestream(aggregates, it_ret);
    assert(it_ret == ret.end());
    return ret;
}

std::string ClientInterface::receive_and_send_5(const std::string &bytes_str) {
    return base64_encode(receive_and_send_5_bytes(base64_decode(bytes_str)));
}