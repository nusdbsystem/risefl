//
// Created by yizheng on 20/3/23.
//
#include <sodium.h>
#include <iostream>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <chrono>
#include <cassert>
#include <vector>
#include <string>
#include <algorithm>
#include <execution>
#include <thread>
#include "../include/ristretto.h"
#include "../include/ristretto_vector.h"
#include "../include/risefl_interface_common.h"
#include "../include/risefl_interface_client.h"
#include "../include/risefl_interface_server.h"
#include <random>

void test_discrete_log() {
    // server generates small mult base table and big mult base vals
    int small_mult_base_table_bit_size = 16;
    std::unordered_map<RistElem, int, RistElemHasher> small_table;
    auto curr = scalar_mult_base(RistScal(-(1 << (small_mult_base_table_bit_size - 1))));
    auto g = scalar_mult_base(c_scal_one);
    for (int i = -(1 << (small_mult_base_table_bit_size - 1)); i < (1 << (small_mult_base_table_bit_size - 1)); i++) {
        small_table.emplace(curr, i);
        add(curr, curr, g);
    }


    assert(discrete_log(scalar_mult_base(RistScal(1 << 10)), small_table, 1 << small_mult_base_table_bit_size) == 1 << 10);
    assert(discrete_log(scalar_mult_base(-RistScal(1 << 10)), small_table, 1 << small_mult_base_table_bit_size) == - (1 << 10));
    assert(discrete_log(scalar_mult_base(RistScal(1 << 20)), small_table, 1 << small_mult_base_table_bit_size) == 1 << 20);
    assert(discrete_log(scalar_mult_base(-RistScal(1 << 20)), small_table, 1 << small_mult_base_table_bit_size) == - (1 << 20));
    assert(discrete_log(scalar_mult_base(RistScal((1 << 24) + 33)), small_table, 1 << small_mult_base_table_bit_size) == (1 << 24) + 33);
    assert(discrete_log(scalar_mult_base(-RistScal((1 << 24) + 33)), small_table, 1 << small_mult_base_table_bit_size) == - (1 << 24) - 33);

    std::cout << "discrete log test passed!\n";
}

// no bad clients
void test_interface_1(int num_clients,
                      int max_malicious_clients,
                      int dim,
                      int num_norm_bound_samples,
                      const std::vector<std::vector<int>> &weight_updates_collection,
                      const RistScal &B) {

    assert(weight_updates_collection.size() == num_clients + 1);
    for (int i = 1; i <= num_clients; i++)
        assert(weight_updates_collection[i].size() == dim);

    int weight_bits = 16;
    int random_normal_bit_shifter = 16;
    int max_bound_sq_bits = 128;

//    RistElemVec hh_bytes(dim), bound_keys_1_bytes(max_bound_sq_bits), bound_keys_2_bytes(max_bound_sq_bits);
    RistElemP3Vec hh(dim);
    RistElemP3Vec bound_keys_1(max_bound_sq_bits);
    RistElemP3Vec bound_keys_2(max_bound_sq_bits);
    std::array<unsigned, 64> random_bytes;
    randombytes_buf(random_bytes.data(), random_bytes.size());

    // random_bytes are agreed between server and clients.
    // Before the entire training process,
    //      use the same algorithm to compute public keys hh_bytes from random_bytes on server and clients,
    //      and reuse hh_bytes throughout all the iterations.
    for (int i = 0; i < hh.size(); i++) {
        hh[i] = p3_from_hash_from_bytes(random_bytes, bytes_from_int(i));
    }
    for (int i = 0; i < max_bound_sq_bits; i++) {
        bound_keys_1[i] = p3_from_hash_from_bytes(random_bytes, bytes_from_int(2 * i + dim));
        bound_keys_2[i] = p3_from_hash_from_bytes(random_bytes, bytes_from_int(2 * i + 1 + dim));
    }

    // server generates small mult base table
    int small_mult_base_table_bit_size = 24;
    std::unordered_map<RistElem, int, RistElemHasher> small_table;
    auto curr = scalar_mult_base(RistScal(-(1 << (small_mult_base_table_bit_size - 1))));
    auto g = scalar_mult_base(c_scal_one);
    for (int i = -(1 << (small_mult_base_table_bit_size - 1)); i < (1 << (small_mult_base_table_bit_size - 1)); i++) {
        small_table.emplace(curr, i);
        add(curr, curr, g);
    }

    // Initialize server
    ServerInterface server(num_clients, max_malicious_clients, dim,
                           weight_bits, random_normal_bit_shifter,
                           num_norm_bound_samples, max_bound_sq_bits);

    // import public keys hh_bytes and bound keys. generate lookup table of multiplications of small steps
    // Needed at 1st iteration, not needed from 2nd iteration onwards.
    server.import_group_element_keys(hh, bound_keys_1, bound_keys_2);
    server.import_small_mult_base_table(small_table);
    assert(server.predicate.hh == hh);
    assert(server.predicate.bound_elem_keys_1 == bound_keys_1);
    assert(server.predicate.bound_elem_keys_2 == bound_keys_2);

    for (int i = -(1 << (small_mult_base_table_bit_size - 1)); i < (1 << (small_mult_base_table_bit_size - 1)); i++) {
        auto val = server.small_mult_base_table.find(scalar_mult_base(RistScal(i)));
        assert(val->first == scalar_mult_base(RistScal(i)));
        assert(val->second == i);
    }



    // server resets malicious clients flags and dispute table
    //  (not needed at 1st iteration, but needed from 2nd iteration onwards)
    server.reset_server_flags();
    assert(server.server_flags == std::vector<int>(num_clients + 1, 0));
    server.reset_dispute_table();
    assert(server.dispute_table ==
           std::vector<std::vector<int>>(num_clients + 1, std::vector<int>(num_clients + 1, 0)));

    // server sets a bound (needed at every iteration)
    server.set_l2_bound_sq(B);
    assert(server.predicate.check_param.l2_param.bound_sq == B);

    // Initialize clients. Client 0 is never used. Only clients 1, 2, ..., num_clients are used.
    std::vector<ClientInterface> client;
    for (int i = 0; i <= num_clients; i++) {
        client.emplace_back(ClientInterface(num_clients, max_malicious_clients, dim,
                                            weight_bits, random_normal_bit_shifter,
                                            num_norm_bound_samples, max_bound_sq_bits,
                                            CHECK_TYPE::L2NORM,
                                            i));
    }

    // client i does:
    for (int i = 1; i <= num_clients; i++) {
        // import updates
        client[i].import_weight_updates(weight_updates_collection[i]);
        assert(client[i].weight_updates.size() == dim);
        for (int l = 0; l < dim; l++)
            assert(client[i].weight_updates[l] == RistScal(weight_updates_collection[i][l]));

        // import public keys hh_bytes and bound keys.
        // Needed at 1st iteration, not needed from 2nd iteration onwards.
        client[i].import_group_element_keys(hh, bound_keys_1, bound_keys_2);
        assert(client[i].predicate.hh == hh);
        assert(client[i].predicate.bound_elem_keys_1 == bound_keys_1);
        assert(client[i].predicate.bound_elem_keys_2 == bound_keys_2);

        // reset malicious clients flags
        // not needed at 1st iteration, but needed from 2nd iteration onwards
        client[i].reset_server_flags();
        assert(client[i].server_flags == std::vector<int>(num_clients + 1, 0));
        client[i].reset_flags();
        assert(client[i].flags == std::vector<int>(num_clients + 1, 0));

        // import the bound (needed at every iteration)
        client[i].set_l2_bound_sq(B);
        assert(client[i].predicate.check_param.l2_param.bound_sq == B);


        // generate Diffe-Hellman key pair
        client[i].generate_dh_key_pair();
    }

    // clients send their Diffe-Hellman public keys to the server
    for (int i = 1; i <= num_clients; i++) {
        server.dh_public_key_collection[i] = client[i].dh_public_key;
    }

    // clients download all the Diffe-Hellman public keys from the server
    for (int i = 1; i <= num_clients; i++) {
        client[i].dh_public_key_collection = server.dh_public_key_collection;
        assert(client[i].dh_public_key == client[i].dh_public_key_collection[i]);
    }

    std::vector<std::vector<unsigned char>> committed_updates_collection(num_clients + 1);
    // client i commit weight updates, generate shares of its blind_to_share, encrypt these shares
    for (int i = 1; i <= num_clients; i++) {
        committed_updates_collection[i] = client[i].commit_updates();
        assert(committed_updates_collection[i].size() == dim * RISTBYTES);
        client[i].generate_shares_and_check_string();
        assert(client[i].shamir_shares.size() == num_clients + 1);
        assert(client[i].shamir_check_string.size() == max_malicious_clients + 1);
        client[i].encrypt_shares();
        assert(client[i].encrypted_shamir_shares.size() == num_clients + 1);
    }

    // client i sends committed weight updates, encrypted shares of its blind_to_share and check string to the server
    for (int i = 1; i <= num_clients; i++) {
        server.import_committed_updates(committed_updates_collection[i], i);
        assert(server.committed_updates_collection.size() == num_clients + 1);
        server.encrypted_shamir_shares_collection[i] = client[i].encrypted_shamir_shares;
        assert(server.encrypted_shamir_shares_collection.size() == num_clients + 1);
        server.shamir_check_string_collection[i] = client[i].shamir_check_string;
        assert(server.shamir_check_string_collection.size() == num_clients + 1);
    }

    // server sends encrypted shamir shares and check strings back to clients
    // server also sends the seed generating aa and linear comb of (aa, hh) to the clients
    server.transpose_encrypted_shamir_shares();
    assert(server.encrypted_shamir_shares_collection_transposed.size() == num_clients + 1);
    server.generate_aa_seed();
    server.generate_aa();
    server.compute_hh_comb();
    server.compute_hh_comb_nonzero_sum();
    assert(sum_single_thread(server.predicate.hh_comb) == server.predicate.hh_comb[0] + server.predicate.hh_comb_nonzero_sum);

    // client i downloads other clients' encrypted shamir shares and check strings
    //  downloads the seed that generates linear combination coefficients,
    //  generates the linear combination coefficient using the seed,
    //  downloads server's computation of hh_comb,
    //  checks the correctness of hh_comb,
    //  and compute the sum_single_thread of hh_comb
    for (int i = 1; i <= num_clients; i++) {
        client[i].other_encrypted_shamir_shares = server.encrypted_shamir_shares_collection_transposed[i];
        for (int j = 1; j <= num_clients; j++) {
            assert(client[i].other_encrypted_shamir_shares[j] == client[j].encrypted_shamir_shares[i]);
        }
        client[i].shamir_check_string_collection = server.shamir_check_string_collection;
        for (int j = 1; j <= num_clients; j++) {
            assert(client[i].shamir_check_string_collection[j] == client[j].shamir_check_string);
        }
        client[i].set_aa_seed(server.aa_seed);
        assert(client[i].aa_seed == server.aa_seed);
        client[i].generate_aa();
        assert(client[i].predicate.aa == server.predicate.aa);
        client[i].set_hh_comb(server.predicate.hh_comb);
        assert(client[i].predicate.hh_comb == server.predicate.hh_comb);
        client[i].check_correctness_hh_comb();
        client[i].compute_hh_comb_nonzero_sum();
        assert(client[i].predicate.hh_comb_nonzero_sum == server.predicate.hh_comb_nonzero_sum);
    }

    // client i generates proof and check integrity of other clients' shamir shares,
    //      flag clients whose shamir shares didn't pass the integrity check
    for (int i = 1; i <= num_clients; i++) {
        client[i].generate_proof();
        client[i].decrypt_shamir_shares();
        for (int j = 1; j <= num_clients; j++) {
            if (j != i)
                assert(client[i].other_shamir_shares[j] == client[j].shamir_shares[i]);
        }
        client[i].check_shamir_share_integrity();
        for (int j = 1; j <= num_clients; j++)
            assert(!client[i].flags[j]);
    }

    // client i sends proof and flags to the server
    for (int i = 1; i <= num_clients; i++) {
        server.proof_collection[i] = client[i].proof;
        server.flags_collection[i] = client[i].flags;
    }

    // server generates dispute table from client flags collection, check clients' proofs
    server.generate_dispute_table();
    for (int i = 1; i <= num_clients; i++)
        for (int j = 1; j <= num_clients; j++)
            assert(!server.dispute_table[i][j]);
    server.check_proofs();
    for (int i = 1; i <= num_clients; i++)
        assert(!server.server_flags[i]);

    // server sends dispute table to clients
    for (int i = 1; i <= num_clients; i++) {
        client[i].dispute_clients = server.dispute_table[i];
        client[i].generate_dispute_shares();
    }

    // client i sends dispute shares to server
    for (int i = 1; i <= num_clients; i++) {
        server.dispute_shares_collection[i] = client[i].dispute_shares;
    }

    // server checks disputes
    server.check_disputes();
    server.transpose_dispute_shares();

    // server sends disputes and flags to clients
    for (int i = 1; i <= num_clients; i++) {
        client[i].other_dispute_shares = server.dispute_shares_collection_transposed[i];
        client[i].server_flags = server.server_flags;
        for (int j = 1; j <= num_clients; j++) {
            assert(!client[i].server_flags[j]);
        }
    }

    // client i uses disputes (if any) to compute non-server-flagged aggregation
    for (int i = 1; i <= num_clients; i++) {
        client[i].update_other_shamir_shares_with_dispute();
        client[i].compute_aggegrated_share();
    }

    // client i sends aggregated share to server
    for (int i = 1; i <= num_clients; i++) {
        server.aggregates_collection[i] = client[i].aggregate;
    }

    // server computes final update
    server.compute_final_update();
    for (int l = 0; l < dim; l++) {
        auto sum = 0;
        for (int i = 1; i <= num_clients; i++)
            sum += weight_updates_collection[i][l];
        assert(server.final_update[l] == sum);
    }

    std::cout << "Interface test 1 passed!\n";
}

// has out of bound clients
// dishonest flagging clients (flags max_malicious_clients+1 clients)
// partial dishonest flagging clients  (flags max_malicious_clients clients)
void test_interface_2(int num_clients,
                      int max_malicious_clients,
                      int dim,
                      int num_norm_bound_samples,
                      const std::vector<std::vector<int>> &weight_updates_collection,
                      const RistScal &B,
                      const std::vector<int> &out_of_bound_flags,
                      const std::vector<int> &dishonest_flagging_flags,
                      const std::vector<int> &partial_dishonest_flagging_flags,
                      const std::vector<int> &bad_share_flags) {

    assert(weight_updates_collection.size() == num_clients + 1);
    for (int i = 1; i <= num_clients; i++)
        assert(weight_updates_collection[i].size() == dim);

    int weight_bits = 16;
    int random_normal_bit_shifter = 16;
    int max_bound_sq_bits = 128;

    RistElemP3Vec hh(dim);
    RistElemP3Vec bound_keys_1(max_bound_sq_bits);
    RistElemP3Vec bound_keys_2(max_bound_sq_bits);
    std::array<unsigned, 64> random_bytes;
    randombytes_buf(random_bytes.data(), random_bytes.size());

    // random_bytes are agreed between server and clients.
    // Before the entire training process,
    //      use the same algorithm to compute public keys hh from random_bytes on server and clients,
    //      and reuse hh throughout all the iterations.
    for (int i = 0; i < hh.size(); i++) {
        hh[i] = p3_from_hash_from_bytes(random_bytes, bytes_from_int(i));
    }
    for (int i = 0; i < max_bound_sq_bits; i++) {
        bound_keys_1[i] = p3_from_hash_from_bytes(random_bytes, bytes_from_int(2 * i + dim));
        bound_keys_2[i] = p3_from_hash_from_bytes(random_bytes, bytes_from_int(2 * i + 1 + dim));
    }

    // server generates small mult base table and big mult base vals
    int small_mult_base_table_bit_size = 24;
    std::unordered_map<RistElem, int, RistElemHasher> small_table;
    auto curr = scalar_mult_base(RistScal(-(1 << (small_mult_base_table_bit_size - 1))));
    auto g = scalar_mult_base(c_scal_one);
    for (int i = -(1 << (small_mult_base_table_bit_size - 1)); i < (1 << (small_mult_base_table_bit_size - 1)); i++) {
        small_table.emplace(curr, i);
        add(curr, curr, g);
    }

    // Initialize server
    ServerInterface server(num_clients, max_malicious_clients, dim,
                           weight_bits, random_normal_bit_shifter,
                           num_norm_bound_samples, max_bound_sq_bits);

    // import public keys hh and bound keys. generate lookup table of multiplications of small steps
    // Needed at 1st iteration, not needed from 2nd iteration onwards.
    server.import_group_element_keys(hh, bound_keys_1, bound_keys_2);
    server.import_small_mult_base_table(small_table);
    assert(server.predicate.hh == hh);
    assert(server.predicate.bound_elem_keys_1 == bound_keys_1);
    assert(server.predicate.bound_elem_keys_2 == bound_keys_2);

    for (int i = -(1 << (small_mult_base_table_bit_size - 1)); i < (1 << (small_mult_base_table_bit_size - 1)); i++) {
        auto val = server.small_mult_base_table.find(scalar_mult_base(RistScal(i)));
        assert(val->first == scalar_mult_base(RistScal(i)));
        assert(val->second == i);
    }



    // server resets malicious clients flags and dispute table
    //  (not needed at 1st iteration, but needed from 2nd iteration onwards)
    server.reset_server_flags();
    assert(server.server_flags == std::vector<int>(num_clients + 1, 0));
    server.reset_dispute_table();
    assert(server.dispute_table ==
           std::vector<std::vector<int>>(num_clients + 1, std::vector<int>(num_clients + 1, 0)));

    // server sets a bound (needed at every iteration)
    server.set_l2_bound_sq(B);
    assert(server.predicate.check_param.l2_param.bound_sq == B);

    // Initialize clients. Client 0 is never used. Only clients 1, 2, ..., num_clients are used.
    std::vector<ClientInterface> client;
    for (int i = 0; i <= num_clients; i++) {
        client.emplace_back(ClientInterface(num_clients, max_malicious_clients, dim,
                                            weight_bits, random_normal_bit_shifter,
                                            num_norm_bound_samples, max_bound_sq_bits,
                                            CHECK_TYPE::L2NORM,
                                            i));
    }

    // client i does:
    for (int i = 1; i <= num_clients; i++) {
        // import updates
        client[i].import_weight_updates(weight_updates_collection[i]);
        assert(client[i].weight_updates.size() == dim);
        for (int l = 0; l < dim; l++)
            assert(client[i].weight_updates[l] == RistScal(weight_updates_collection[i][l]));

        // import public keys hh and bound keys.
        // Needed at 1st iteration, not needed from 2nd iteration onwards.
        client[i].import_group_element_keys(hh, bound_keys_1, bound_keys_2);
        assert(client[i].predicate.hh == hh);
        assert(client[i].predicate.bound_elem_keys_1 == bound_keys_1);
        assert(client[i].predicate.bound_elem_keys_2 == bound_keys_2);

        // reset malicious clients flags
        // not needed at 1st iteration, but needed from 2nd iteration onwards
        client[i].reset_server_flags();
        assert(client[i].server_flags == std::vector<int>(num_clients + 1, 0));
        client[i].reset_flags();
        assert(client[i].flags == std::vector<int>(num_clients + 1, 0));

        // import the bound (needed at every iteration)
        client[i].set_l2_bound_sq(B);
        assert(client[i].predicate.check_param.l2_param.bound_sq == B);


        // generate Diffe-Hellman key pair
        client[i].generate_dh_key_pair();
    }

    // clients send their Diffe-Hellman public keys to the server
    for (int i = 1; i <= num_clients; i++) {
        server.dh_public_key_collection[i] = client[i].dh_public_key;
    }

    // clients download all the Diffe-Hellman public keys from the server
    for (int i = 1; i <= num_clients; i++) {
        client[i].dh_public_key_collection = server.dh_public_key_collection;
        assert(client[i].dh_public_key == client[i].dh_public_key_collection[i]);
    }

    std::vector<std::vector<unsigned char>> committed_updates_collection(num_clients + 1);
    // client i commit weight updates, generate shares of its blind_to_share, encrypt these shares
    for (int i = 1; i <= num_clients; i++) {
        committed_updates_collection[i] = client[i].commit_updates();
        assert(committed_updates_collection[i].size() == dim * RISTBYTES);
        client[i].generate_shares_and_check_string();
        assert(client[i].shamir_shares.size() == num_clients + 1);
        assert(client[i].shamir_check_string.size() == max_malicious_clients + 1);
        if (bad_share_flags[i]) {
            client[i].shamir_shares[1].scalar[0] += 1000;
        }
        client[i].encrypt_shares();
        assert(client[i].encrypted_shamir_shares.size() == num_clients + 1);
    }

    // client i sends committed weight updates, encrypted shares of its blind_to_share and check string to the server
    for (int i = 1; i <= num_clients; i++) {
        server.import_committed_updates(committed_updates_collection[i], i);
        assert(server.committed_updates_collection.size() == num_clients + 1);
        server.encrypted_shamir_shares_collection[i] = client[i].encrypted_shamir_shares;
        assert(server.encrypted_shamir_shares_collection.size() == num_clients + 1);
        server.shamir_check_string_collection[i] = client[i].shamir_check_string;
        assert(server.shamir_check_string_collection.size() == num_clients + 1);
    }

    // server sends encrypted shamir shares and check strings back to clients
    // server also sends the seed generating aa and linear comb of (aa, hh) to the clients
    server.transpose_encrypted_shamir_shares();
    assert(server.encrypted_shamir_shares_collection_transposed.size() == num_clients + 1);
    server.generate_aa_seed();
    server.generate_aa();
    server.compute_hh_comb();
    server.compute_hh_comb_nonzero_sum();
    assert(sum_single_thread(server.predicate.hh_comb) == server.predicate.hh_comb[0] + server.predicate.hh_comb_nonzero_sum);

    // client i downloads other clients' encrypted shamir shares and check strings
    //  downloads the seed that generates linear combination coefficients,
    //  generates the linear combination coefficient using the seed,
    //  downloads server's computation of hh_comb,
    //  checks the correctness of hh_comb,
    //  and compute the sum_single_thread of hh_comb
    for (int i = 1; i <= num_clients; i++) {
        client[i].other_encrypted_shamir_shares = server.encrypted_shamir_shares_collection_transposed[i];
        for (int j = 1; j <= num_clients; j++) {
            assert(client[i].other_encrypted_shamir_shares[j] == client[j].encrypted_shamir_shares[i]);
        }
        client[i].shamir_check_string_collection = server.shamir_check_string_collection;
        for (int j = 1; j <= num_clients; j++) {
            assert(client[i].shamir_check_string_collection[j] == client[j].shamir_check_string);
        }
        client[i].set_aa_seed(server.aa_seed);
        assert(client[i].aa_seed == server.aa_seed);
        client[i].generate_aa();
        assert(client[i].predicate.aa == server.predicate.aa);
        client[i].set_hh_comb(server.predicate.hh_comb);
        assert(client[i].predicate.hh_comb == server.predicate.hh_comb);
        client[i].check_correctness_hh_comb();
        client[i].compute_hh_comb_nonzero_sum();
        assert(client[i].predicate.hh_comb_nonzero_sum == server.predicate.hh_comb_nonzero_sum);
    }

    // client i generates proof and check integrity of other clients' shamir shares,
    //      flag clients whose shamir shares didn't pass the integrity check
    for (int i = 1; i <= num_clients; i++) {
        client[i].generate_proof();
        client[i].decrypt_shamir_shares();
        for (int j = 1; j <= num_clients; j++) {
            if (j != i)
                assert(client[i].other_shamir_shares[j] == client[j].shamir_shares[i]);
        }
        client[i].check_shamir_share_integrity();
        for (int j = 1; j <= num_clients; j++)
            assert(client[i].flags[j] == ((i == 1) && bad_share_flags[j]));

        // dishonest flags
        if (dishonest_flagging_flags[i]) {
            for (int j = 1; j <= max_malicious_clients + 1; j++)
                client[i].flags[j] = 1;
        }
        if (partial_dishonest_flagging_flags[i]) {
            for (int j = 1; j <= max_malicious_clients; j++)
                client[i].flags[j] = 1;
        }
    }

    // client i sends proof and flags to the server
    for (int i = 1; i <= num_clients; i++) {
        server.proof_collection[i] = client[i].proof;
        server.flags_collection[i] = client[i].flags;
    }

    // server generates dispute table from client flags collection, check clients' proofs
    server.generate_dispute_table();
    for (int i = 1; i <= num_clients; i++)
        for (int j = 1; j <= num_clients; j++) {
            if (dishonest_flagging_flags[i] || dishonest_flagging_flags[j] || out_of_bound_flags[i] ||
                out_of_bound_flags[j])
                assert(!server.dispute_table[j][i]);
            else if (partial_dishonest_flagging_flags[i]) {
                assert(server.dispute_table[j][i] == (j <= max_malicious_clients));
            } else if (bad_share_flags[j] && i == 1) {
                assert(server.dispute_table[j][i]);
            } else
                assert(!server.dispute_table[j][i]);
        }
    server.check_proofs();
    for (int i = 1; i <= num_clients; i++)
        assert(server.server_flags[i] == out_of_bound_flags[i] || dishonest_flagging_flags[i]);
//    for (int i = 1; i <= num_clients; i++)
//        assert(!server.server_flags[i]);

    // server sends dispute table to clients
    for (int i = 1; i <= num_clients; i++) {
        client[i].dispute_clients = server.dispute_table[i];

        if (i <= max_malicious_clients && !bad_share_flags[i]) {
            for (int j = 1; j <= num_clients; j++) {
                assert(client[i].dispute_clients[j] == partial_dishonest_flagging_flags[j] || !out_of_bound_flags[j]);
            }
        } else if (bad_share_flags[i]) {
            for (int j = 1; j <= num_clients; j++) {
                assert(client[i].dispute_clients[j] == (j == 1));
            }
        } else {
            for (int j = 1; j <= num_clients; j++) {
                assert(client[i].dispute_clients[j] == 0);
            }
        }
        client[i].generate_dispute_shares();
        if (i <= max_malicious_clients) {
            for (int j = 1; j <= num_clients; j++) {
                if (client[i].dispute_clients[j])
                    assert(client[i].dispute_shares[j] == client[i].shamir_shares[j]);
            }
        }
    }

    // client i sends dispute shares to server
    for (int i = 1; i <= num_clients; i++) {
        server.dispute_shares_collection[i] = client[i].dispute_shares;
    }

    // server checks disputes
    server.check_disputes();
    for (int i = 1; i <= num_clients; i++) {
        assert(server.server_flags[i] == out_of_bound_flags[i] || dishonest_flagging_flags[i] || bad_share_flags[i]);
    }
    server.transpose_dispute_shares();

    // server sends disputes and flags to clients
    for (int i = 1; i <= num_clients; i++) {
        client[i].other_dispute_shares = server.dispute_shares_collection_transposed[i];
        if (partial_dishonest_flagging_flags[i]) {
            for (int j = 1; j <= max_malicious_clients; j++) {
                if (!server.server_flags[i])
                    assert(client[i].other_dispute_shares[j] == client[j].shamir_shares[i]);
            }
        }
        client[i].server_flags = server.server_flags;
//        for (int j = 1; j <= num_clients; j++) {
//            assert(!client[i].server_flags[j]);
//        }
    }

    // client i uses disputes (if any) to compute non-server-flagged aggregation
    for (int i = 1; i <= num_clients; i++) {
        client[i].update_other_shamir_shares_with_dispute();
        client[i].compute_aggegrated_share();
    }

    // client i sends aggregated share to server
    for (int i = 1; i <= num_clients; i++) {
        server.aggregates_collection[i] = client[i].aggregate;
    }

    // server computes final update
    server.compute_final_update();
    for (int l = 0; l < dim; l++) {
        auto sum = 0;
        for (int i = 1; i <= num_clients; i++)
            if (!server.server_flags[i])
                sum += weight_updates_collection[i][l];
        assert(server.final_update[l] == sum);
    }

    std::cout << "Interface test 2 passed!\n";
}


void test_interface() {

    std::cout << "testing interface!\n";

    test_discrete_log();

    int num_clients = 3;
    int max_malicious_clients = 1;
    int dim = 2;
    int num_norm_bound_samples = 25;
    std::vector<std::vector<int>> weight_updates_collection{{0, 0}, // weight_updates_collection[0] is never used
                                                            {0, 1},
                                                            {3, 4},
                                                            {6, -9}};
    RistScal B(1L << 50);

    test_interface_1(num_clients, max_malicious_clients, dim, num_norm_bound_samples,
                     weight_updates_collection, B);

    num_clients = 10;
    max_malicious_clients = 3;
    weight_updates_collection.resize(num_clients + 1);

    std::mt19937 gen;
    std::uniform_int_distribution<> distrib(-10, 10);
    for (int i = 1; i <= num_clients; i++) {
        weight_updates_collection[i].resize(dim);
        for (int j = 0; j <= dim; j++)
            weight_updates_collection[i][j] = distrib(gen);
    }

    test_interface_1(num_clients, max_malicious_clients, dim, num_norm_bound_samples,
                     weight_updates_collection, B);

    std::vector<int> out_of_bound_flags(num_clients + 1, 0);
    out_of_bound_flags[10] = 1;
    weight_updates_collection[10][0] = 1L << 30;

    std::vector<int> dishonest_flagging_flags(num_clients + 1, 0);
    dishonest_flagging_flags[5] = 1; // 5 flags 1,2,3,4

    std::vector<int> partial_dishonest_flagging_flags(num_clients + 1, 0);
    partial_dishonest_flagging_flags[4] = 1; // 4 flags 1,2,3
    partial_dishonest_flagging_flags[3] = 1; // 3 flags 1,2,3

    std::vector<int> bad_share_flags(num_clients + 1, 0);
    bad_share_flags[8] = 1; // 8 sends 1 a bad share

    test_interface_2(num_clients, max_malicious_clients, dim, num_norm_bound_samples,
                     weight_updates_collection, B,
                     out_of_bound_flags,
                     dishonest_flagging_flags,
                     partial_dishonest_flagging_flags,
                     bad_share_flags);

    std::cout << "Interface test passed!\n";
}