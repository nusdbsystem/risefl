//
// Created by yizheng on 31/3/23.
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
#include "../include/di_zkp_interface_common.h"
#include "../include/di_zkp_interface_client.h"
#include "../include/di_zkp_interface_server.h"
#include "../include/bulletin.h"
#include <random>

void test_interface_abstract_1(int num_clients,
                               int max_malicious_clients,
                               int dim,
                               int num_blinds_per_group_element,
                               int num_norm_bound_samples,
                               const std::vector<std::vector<float>> &weight_updates_collection,
                               float norm_bound) {
    int weight_bits = 16;
    int random_normal_bit_shifter = 16;
    int inner_prod_bound_bits = weight_bits + random_normal_bit_shifter + 4;
    int max_bound_sq_bits = 2 * inner_prod_bound_bits + 64;

    CheckParamFloat check_param_float(CHECK_TYPE::L2NORM);
    check_param_float.l2_param.bound = norm_bound;

    // create a seed generates public group element keys
    std::vector<unsigned char> random_bytes(64);
    randombytes_buf(random_bytes.data(), random_bytes.size());
    std::string random_bytes_str = base64_encode(random_bytes.data(), random_bytes.size());

    // Initialize server
    ServerInterface server(num_clients, max_malicious_clients, dim, num_blinds_per_group_element,
                           weight_bits, random_normal_bit_shifter,
                           num_norm_bound_samples, inner_prod_bound_bits, max_bound_sq_bits,
                           CHECK_TYPE::L2NORM);
    server.initialize_from_seed(random_bytes_str);
    server.initialize_new_iteration(check_param_float);

    // Initialize clients. Client 0 is never used. Only clients 1, 2, ..., num_clients are used.
    std::vector<ClientInterface> client;
    std::vector<SignPubKey> bul_pub_keys(num_clients + 1);
    std::vector<SignPrvKey> bul_prv_keys(num_clients + 1);
    for (int i = 0; i <= num_clients; i++) {
        auto temp = gen_sign_key_pair();
        bul_pub_keys[i] = temp.first;
        bul_prv_keys[i] = temp.second;
    }
    for (int i = 0; i <= num_clients; i++) {
        client.emplace_back(ClientInterface(num_clients, max_malicious_clients, dim, num_blinds_per_group_element,
                                            weight_bits, random_normal_bit_shifter,
                                            num_norm_bound_samples, inner_prod_bound_bits, max_bound_sq_bits,
                                            CHECK_TYPE::L2NORM,
                                            i, bul_pub_keys, bul_prv_keys[i]));
    }
    for (int i = 1; i <= num_clients; i++) {
        client[i].initialize_from_seed(random_bytes_str);
    }

    // step 1: clients send messages to server
    for (int i = 1; i <= num_clients; i++) {
        server.receive_1(client[i].send_1(check_param_float, weight_updates_collection[i]), i);
    }

    // step 2: server sends (the same) messages to clients, and gets messages back
    auto bytes_str_sent_2 = server.send_2();
    for (int i = 1; i <= num_clients; i++) {
        server.receive_2(client[i].receive_and_send_2(bytes_str_sent_2), i);
    }

    // server can concurrently do computations while receiving messages in steps 1 & 2, before step 3
    server.concurrent_process_before_send_3();

    // step 3: server sends messages to clients, and gets messages back
    for (int i = 1; i <= num_clients; i++) {
        server.receive_3(client[i].receive_and_send_3(server.send_3(i)), i);
    }

    // step 4: server sends messages to clients, and gets messages back
    server.process_before_send_4();
    for (int i = 1; i <= num_clients; i++) {
        server.receive_4(client[i].receive_and_send_4(server.send_4(i)), i);
    }

    // step 5: server sends messages to clients, and gets messages back
    server.process_before_send_5();
    for (int i = 1; i <= num_clients; i++) {
        server.receive_5(client[i].receive_and_send_5(server.send_5(i)), i);
    }

    // finish one iteration. aggregated sum is in server.final_update
    server.finish_iteration();

    for (int l = 0; l < dim; l++) {
        float sum = 0;
        for (int i = 1; i <= num_clients; i++)
            sum += weight_updates_collection[i][l];
        assert(abs(server.final_update_float[l] - sum) < 1e-4);
    }

    std::cout << "test interface abstract 1 success!\n";
}

void test_interface_abstract() {

    int num_clients = 3;
    int max_malicious_clients = 1;
    int dim = 2;
    int num_blinds_per_group_element = 1;
    int num_norm_bound_samples = 3;
    std::vector<std::vector<float>> weight_updates_collection{{0, 0}, // weight_updates_collection[0] is never used
                                                            {0, 1},
                                                            {0.3, 0.4},
                                                            {0.6, -0.9}};
//    RistScal B(1L << 50);
    float norm_bound = 2;
//    float standard_deviation_factor = 6;

    test_interface_abstract_1(num_clients, max_malicious_clients, dim, num_blinds_per_group_element, num_norm_bound_samples,
                              weight_updates_collection, norm_bound);

    std::cout << "test interface abstract success!\n";
}

