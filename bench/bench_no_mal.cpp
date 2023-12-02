//
// Created by yizheng on 4/4/23.
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
#include "../include/di_zkp_interface_client.h"
#include "../include/di_zkp_interface_common.h"
#include "../include/di_zkp_interface_server.h"
#include "simulate.h"
#include "../include/bulletin.h"

double average_from_1(std::vector<double> vv) {
    return std::reduce(vv.begin() + 1, vv.end()) / (vv.size() - 1);
}

struct BenchResult {
    int num_clients;

    std::vector<double> c1, c2, c3, c4, c5;
    double s2s, s3p, s4p, s5p, s6;
    std::vector<double> s1r, s2r, s3r, s3s, s4r, s4s, s5r, s5s;

    double c1a, c2a, c3a, c4a, c5a;
    double s1ra, s2ra, s3ra, s3sa, s4ra, s4sa, s5ra, s5sa;

    std::shared_ptr<double> cost_commitment;
    std::shared_ptr<double> cost_encrypt_vss;
    std::shared_ptr<double> cost_generate_proof;
    std::shared_ptr<double> cost_verify_shamir;

    BenchResult(int num_clients) :
            num_clients(num_clients),
            c1(num_clients + 1, 0),
            c2(num_clients + 1, 0),
            c3(num_clients + 1, 0),
            c4(num_clients + 1, 0),
            c5(num_clients + 1, 0),
            s1r(num_clients + 1, 0),
            s2r(num_clients + 1, 0),
            s3s(num_clients + 1, 0),
            s3r(num_clients + 1, 0),
            s4s(num_clients + 1, 0),
            s4r(num_clients + 1, 0),
            s5s(num_clients + 1, 0),
            s5r(num_clients + 1, 0),
            s6(0),
            cost_commitment(std::make_shared<double>(0)),
            cost_encrypt_vss(std::make_shared<double>(0)),
            cost_generate_proof(std::make_shared<double>(0)),
            cost_verify_shamir(std::make_shared<double>(0))
    {}

    void compute_average() {
        c1a = average_from_1(c1);
        c2a = average_from_1(c2);
        c3a = average_from_1(c3);
        c4a = average_from_1(c4);
        c5a = average_from_1(c5);
        s1ra = average_from_1(s1r);
        s2ra = average_from_1(s2r);
        s3ra = average_from_1(s3r);
        s3sa = average_from_1(s3s);
        s4ra = average_from_1(s4r);
        s4sa = average_from_1(s4s);
        s5ra = average_from_1(s5r);
        s5sa = average_from_1(s5s);
    }

    void print() {
        std::cout << "client average: \n"
                  << "\tstep 1: " << c1a << " seconds \n"
                  << "\tstep 2 (Commitment + VSS): " << c2a << " seconds \n"
                  << "\t\tCommitment + VSS: " << *cost_commitment / num_clients << " seconds \n"
                  << "\t\tEncrypt shares per client: " << *cost_encrypt_vss / num_clients / (num_clients - 1) << " seconds \n"
                  << "\tstep 3 (Generate Proof + VSS check): " << c3a << " seconds \n"
                  << "\t\tGenerate proof: " << *cost_generate_proof / num_clients << " seconds \n"
                  << "\t\tVerify VSS per client: " << *cost_verify_shamir / num_clients / (num_clients - 1) << " seconds \n"
                  << "\tstep 4: " << c4a << " seconds \n"
                  << "\tstep 5: " << c5a << " seconds \n"
                  << "server average: \n"
                  << "\t step 1 receive: " << s1ra << " seconds \n"
                  << "\t step 2 receive (Store commitments): " << s2ra << " seconds \n"
                  << "\t step 2 send (total): " << s2s << " seconds \n"
                  << "\t step 3 prepare (Compute combinations of group element keys: concurrent, total): " << s3p << " seconds \n"
                  << "\t step 3 receive (Verify proofs): " << s3ra << " seconds \n"
                  << "\t step 3 send: " << s3sa << " seconds \n"
                  << "\t step 4 prepare (total): " << s4p << " seconds \n"
                  << "\t step 4 receive " << s4ra << " seconds \n"
                  << "\t step 4 send " << s4sa << " seconds \n"
                  << "\t step 5 prepare (total): " << s5p << " seconds \n"
                  << "\t step 5 receive " << s5ra << " seconds \n"
                  << "\t step 5 send " << s5sa << " seconds \n"
                  << "\t step 6 (Aggregation: total) " << s6 << " seconds \n";
    }
    void print2() {
        std::cout << "client average: \n"
                  << "\t\tCommitment + VSS: " << *cost_commitment << " seconds \n"
                  << "\t\tEncrypt shares total: " << *cost_encrypt_vss * (num_clients - 1) << " seconds \n"
                  << "\t\tGenerate proof: " << *cost_generate_proof << " seconds \n"
                  << "\t\tVerify VSS total: " << *cost_verify_shamir   * (num_clients - 1) << " seconds \n"
                  << "\t\t\tTotal: " << (*cost_commitment) + (*cost_encrypt_vss * (num_clients - 1)) + (*cost_generate_proof) + (*cost_verify_shamir   * (num_clients - 1)) << " seconds \n"
                  << "server average: \n"
                  << "\t step 3 prepare (Compute combinations of group element keys: concurrent, total): " << s3p << " seconds \n"
                  << "\t step 3 receive (Verify proofs total): " << s3r[1] * num_clients << " seconds \n"
                  << "\t step 6 (Aggregation: total) " << s6 << " seconds \n"
                  << "\t\tTotal: " <<  s3p + s3r[1] * num_clients + s6 <<" seconds \n";
    }
};

void clear_aa_pos_and_hh_and_precomp(ClientInterface& client) {
//    client.predicate.aa_0.clear();
//    client.predicate.aa_pos.clear();
    client.predicate.hh_precomp.clear();
    client.predicate.hh.clear();
//    client.predicate.coeffs.clear();
//    client.predicate.coeffs_bound.clear();
}

void recover_aa_pos(ClientInterface& client, const CommonInterface& r) {
//    client.predicate.aa_0 = r.predicate.aa_0;
//    client.predicate.aa_pos = r.predicate.aa_pos;
}

void recover_hh(ClientInterface& client, const CommonInterface& r) {
    client.predicate.hh = r.predicate.hh;
}

void recover_precomp(ClientInterface& client, const CommonInterface& r) {
    client.predicate.hh_precomp = r.predicate.hh_precomp;
}

//BenchResult bench_no_mal_1(int num_clients,
//                           int max_malicious_clients,
//                           int dim,
//                           int num_blinds_per_group_element,
//                           int num_norm_bound_samples,
//                           const RistScal &B,
//                           int weight_bits,
//                           int random_normal_bit_shifter,
//                           int inner_prod_bound_bits,
//                           int max_bound_sq_bits,
//                           bool b_precomp) {
//    if (b_precomp)
//        std::cout << "Using precomputation!\n";
//    else
//        std::cout << "Not using precomputation!\n";
//
//    // create a seed generates public group element keys
//    std::vector<unsigned char> random_bytes(64);
//    randombytes_buf(random_bytes.data(), random_bytes.size());
//    std::string random_bytes_str = base64_encode(random_bytes.data(), random_bytes.size());
//
//    // Initialize server
//    ServerInterface server(num_clients, max_malicious_clients, dim, num_blinds_per_group_element,
//                           weight_bits, random_normal_bit_shifter,
//                           num_norm_bound_samples, inner_prod_bound_bits, max_bound_sq_bits);
//    server.initialize_from_seed(random_bytes_str);
//    server.initialize_new_iteration(B);
//
//    std::cout << "server initialization finish!\n";
//
//    // Initialize clients. Client 0 is never used. Only clients 1, 2, ..., num_clients are used.
//    std::vector<ClientInterface> client;
//    for (int i = 0; i <= num_clients; i++) {
//        client.emplace_back(ClientInterface(num_clients, max_malicious_clients, dim, num_blinds_per_group_element,
//                                            weight_bits, random_normal_bit_shifter,
//                                            num_norm_bound_samples, inner_prod_bound_bits, max_bound_sq_bits,
//                                            i,
//                                            b_precomp));
//        clear_aa_pos_and_hh_and_precomp(client.back()); // clear aa_pos to save memory
////        std::cout << "clients initialization " << i << " finish!\n";
//
//    }
//    for (int i = 1; i <= num_clients; i++) {
//        client[i].predicate.bound_elem_keys_1 = server.predicate.bound_elem_keys_1;
//        client[i].predicate.bound_elem_keys_2 = server.predicate.bound_elem_keys_2;
//        client[i].predicate.square_key = server.predicate.square_key;
////        client[i].initialize_from_seed(random_bytes);
//        clear_aa_pos_and_hh_and_precomp(client[i]);
//    }
//
//    std::cout << "clients initialization finish!\n";
//
//    std::vector<long> weight_update(dim);
//    std::string message;
//    BenchResult res(num_clients);
//
//    // step 1: clients send messages to server
//    for (int i = 1; i <= num_clients; i++) {
////        rand_init(weight_update, - (1 << (weight_bits - 1)), 1 << (weight_bits - 1));
//
//        std::random_device rd{};
//        std::mt19937 gen{rd()};
//        std::normal_distribution<> d{0, (1 << (weight_bits - 1)) / sqrt(dim)};
//        for (auto && w : weight_update) {
//            w = static_cast<long>(d(gen));
//        }
//
//        recover_aa_pos(client[i], server);
//        recover_hh(client[i], server);
//        res.c1[i] = measure_time([&](){ message = client[i].send_1_internal(B, weight_update); });
//        clear_aa_pos_and_hh_and_precomp(client[i]);
//        res.s1r[i] = measure_time([&]() { server.receive_1(message, i); });
//    }
//
//    std::cout << "step 1 finish!\n";
//
//    // step 2: server sends (the same) messages to clients, and gets messages back
//    res.s2s = measure_time([&]() { message = server.send_2(); });
//    auto bytes_sent_2 = message;
//    for (int i = 1; i <= num_clients; i++) {
//        recover_aa_pos(client[i], server);
//        recover_precomp(client[i], server);
//        recover_hh(client[i], server);
//        res.c2[i] = measure_time([&]() { message = client[i].receive_and_send_2(bytes_sent_2); });
//        clear_aa_pos_and_hh_and_precomp(client[i]);
//        res.s2r[i] = measure_time([&]() { server.receive_2(message, i); });
//    }
//
//    std::cout << "step 2 finish!\n";
//
//    // server can concurrently do computations while receiving messages in steps 1 & 2, before step 3
//    res.s3p = measure_time([&]() { server.concurrent_process_before_send_3(); });
//
//    // step 3: server sends messages to clients, and gets messages back
//    for (int i = 1; i <= num_clients; i++) {
//        res.s3s[i] = measure_time([&]() { message = server.send_3(i); });
//        recover_aa_pos(client[i], server);
//        recover_hh(client[i], server);
//        res.c3[i] = measure_time([&]() { message = client[i].receive_and_send_3(message); });
//        clear_aa_pos_and_hh_and_precomp(client[i]);
//        res.s3r[i] = measure_time([&]() { server.receive_3(message, i); });
//    }
//
//    std::cout << "step 3 finish!\n";
//
//    // step 4: server sends messages to clients, and gets messages back
//    res.s4p = measure_time([&]() { server.process_before_send_4(); });
//    for (int i = 1; i <= num_clients; i++) {
//        res.s4s[i] = measure_time([&]() { message = server.send_4(i); });
//        recover_aa_pos(client[i], server);
//        res.c4[i] = measure_time([&]() { message = client[i].receive_and_send_4(message); });
//        clear_aa_pos_and_hh_and_precomp(client[i]);
//        res.s4r[i] = measure_time([&]() { server.receive_4(message, i); });
//    }
//
//    std::cout << "step 4 finish!\n";
//
//    // step 5: server sends messages to clients, and gets messages back
//    res.s5p = measure_time([&]() { server.process_before_send_5(); });
//    for (int i = 1; i <= num_clients; i++) {
//        res.s5s[i] = measure_time([&]() { message = server.send_5(i); });
//        recover_aa_pos(client[i], server);
//        res.c5[i] = measure_time([&]() { message = client[i].receive_and_send_5(message); });
//        clear_aa_pos_and_hh_and_precomp(client[i]);
//        res.s5r[i] = measure_time([&]() { server.receive_5(message, i); });
//    }
//
//    std::cout << "step 5 finish!\n";
//
//    // finish one iteration. aggregated sum is in server.final_update
//    res.s6 = measure_time([&]() { server.finish_iteration(); });
//
////    for (int l = 0; l < dim; l++) {
////        auto sum = 0;
////        for (int i = 1; i <= num_clients; i++)
////            sum += client[i].weight_updates[l];
////        assert(server.final_update[l] == sum);
////    }
//
//    return res;
//}


BenchResult bench_no_mal_2(int num_clients,
                           int max_malicious_clients,
                           int dim,
                           int num_blinds_per_group_element,
                           int num_norm_bound_samples,
                           const RistScal &B,
                           int weight_bits,
                           int random_normal_bit_shifter,
                           int inner_prod_bound_bits,
                           int max_bound_sq_bits,
                           bool b_precomp) {
    if (b_precomp)
        std::cout << "Using precomputation!\n";
    else
        std::cout << "Not using precomputation!\n";

    // create a seed generates public group element keys
    std::vector<unsigned char> random_bytes(64);
    randombytes_buf(random_bytes.data(), random_bytes.size());
    std::string random_bytes_str = base64_encode(random_bytes.data(), random_bytes.size());

    CheckParam check_param(CHECK_TYPE::L2NORM);
    check_param.l2_param.bound_sq = B;
    CheckParamFloat check_param_float(CHECK_TYPE::L2NORM);
    check_param_float.l2_param.bound = 1.0;

    // Initialize server
    ServerInterface server(1, max_malicious_clients, dim, num_blinds_per_group_element,
                           weight_bits, random_normal_bit_shifter,
                           num_norm_bound_samples, inner_prod_bound_bits, max_bound_sq_bits,
                           CHECK_TYPE::L2NORM,
                           b_precomp);
    server.num_clients = num_clients;
    server.initialize_from_seed(random_bytes_str);
    server.initialize_new_iteration(check_param_float);

    std::cout << "server initialization finish!\n";

    std::vector<SignPubKey> bul_pub_keys(num_clients + 1);
    std::vector<SignPrvKey> bul_prv_keys(num_clients + 1);
    for (int i = 0; i <= num_clients; i++) {
        auto temp = gen_sign_key_pair();
        bul_pub_keys[i] = temp.first;
        bul_prv_keys[i] = temp.second;
    }

    // Initialize clients. Client 0 is never used. Only clients 1, 2, ..., num_clients are used.
    ClientInterface client(num_clients, max_malicious_clients, dim, num_blinds_per_group_element,
                           weight_bits, random_normal_bit_shifter,
                           num_norm_bound_samples, inner_prod_bound_bits, max_bound_sq_bits,
                           CHECK_TYPE::L2NORM,
                           1, bul_pub_keys, bul_prv_keys[1],
                           b_precomp);
    client.initialize_from_seed(random_bytes_str);


    std::cout << "clients initialization finish!\n";

    std::vector<long> weight_update(dim);
//    std::vector<unsigned char> message;
    BenchResult res(num_clients);

    // step 1: clients send messages to server
//        rand_init(weight_update, - (1 << (weight_bits - 1)), 1 << (weight_bits - 1));

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d{0, (1 << (weight_bits - 1)) / sqrt(dim)};
    for (auto && w : weight_update) {
        w = static_cast<long>(d(gen));
    }

    client.send_1_internal(check_param, weight_update);

    std::cout << "step 1 finish!\n";

    RistP3VecAndBytes commitment(dim);
    auto it_commit = commitment.bytes.begin();

    // step 2: server sends (the same) messages to clients, and gets messages back
//    client.dh_public_key_collection[1] = client.dh_public_key;
    measure_time_and_add_to_bench([&]() {
                                      rand_init(client.blinds_to_share);
                                      if (client.predicate.b_precomp)
                                          pedersen_commit_p3_to_bytes_from_precomp(it_commit, client.weight_updates, weight_bits, client.predicate.hh_precomp,
                                                                                   client.blinds_to_share);
                                      else {
                                          assert(client.blinds_to_share.size() == 1);
                                          pedersen_commit_p3_to_bytes(it_commit,
                                                                      client.weight_updates,
                                                                      weight_bits,
                                                                      client.predicate.hh, client.blinds_to_share[0]);
                                      }
                                      client.generate_batch_shares_and_check_string();
                                  },
                                  res.cost_commitment);

    measure_time_and_add_to_bench([&]() {
                                      for (int i = 0; i < client.predicate.num_blinds_per_group_element; i++) {
                                          client.encrypted_shamir_shares[1][i] = encrypt(
                                                  std::vector<unsigned char>(
                                                          client.batch_shamir_share_with_check.shares[i][1].scalar.begin(),
                                                          client.batch_shamir_share_with_check.shares[i][1].scalar.end()),
                                                  client.dh_public_key, client.dh_private_key);
                                      }
                                  },
                                  res.cost_encrypt_vss);

    auto const_it_commit = commitment.bytes.cbegin();
    p3_vec_from_bytestream(server.committed_updates_collection[1], const_it_commit);

    std::cout << "step 2 finish!\n";

    // server can concurrently do computations while receiving messages in steps 1 & 2, before step 3
    res.s3p = measure_time([&]() { server.concurrent_process_before_send_3(); });

    // step 3: server sends messages to clients, and gets messages back
    client.other_encrypted_shamir_shares[1] = client.encrypted_shamir_shares[1];
    client.shamir_check_string_collection[1] = client.batch_shamir_share_with_check.check_strings;

    client.predicate.aa_seed = server.predicate.aa_seed;
    client.predicate.hh_comb = server.predicate.hh_comb;
    client.predicate.hh = server.predicate.hh;
    client.dh_public_key_collection = server.dh_public_key_collection;
    measure_time_and_add_to_bench([&]() {
//                                      client.generate_aa();
//                                      client.check_correctness_hh_comb();
                                      client.generate_proof();
                                  },
                                  res.cost_generate_proof);

    measure_time_and_add_to_bench([&]() {
                                      for (int i = 0; i < client.predicate.num_blinds_per_group_element; i++) {
                                          auto decrpyted = decrypt(client.other_encrypted_shamir_shares[1][i],
                                                                   client.dh_public_key, client.dh_private_key);
                                          std::copy(decrpyted.begin(), decrpyted.end(),
                                                    client.other_shamir_shares[1][i].scalar.data());
                                      }
                                      b_shamir_verify(client.other_shamir_shares[1],
                                                      client.shamir_check_string_collection[1].elems,
                                                      max_malicious_clients + 1,
                                                      1);
                                  },
                                  res.cost_verify_shamir);
    res.s3r[1] = measure_time([&]() {
        server.proof_collection[1] = client.proof;
        server.shamir_check_string_collection[1] = client.batch_shamir_share_with_check.check_strings;
        server.check_proof(1);
        assert(!server.server_flags[1]);
    });

    std::cout << "step 3 finish!\n";

    // step 4: server sends messages to clients, and gets messages back

    std::cout << "step 4 finish!\n";

    // step 5: server sends messages to clients, and gets messages back

    std::cout << "step 5 finish!\n";

    // finish one iteration. aggregated sum is in server.final_update
    int pieces = 1000;
    for (int i = 0; i < dim / pieces; i++) {
        RistScalVec agg_blinds(num_clients);
        rand_init(agg_blinds);
        std::vector<int> uu_int_agg(num_clients);
        for (auto && w : uu_int_agg) {
            w = static_cast<int>(d(gen));
        }

        RistP3VecAndBytes ped_agg(num_clients);

        auto ped_it_agg = ped_agg.bytes.begin();

        if (b_precomp) {
            pedersen_commit_p3_to_bytes_from_precomp(ped_it_agg, uu_int_agg, 16, server.predicate.hh_precomp[0], agg_blinds);
            ped_agg.fill_elems();
        }
        else {
            for (int j = 0; j < num_clients; j++) {
                pedersen_commit(ped_agg.elems[j], uu_int_agg[j], 16, server.predicate.hh[0], agg_blinds[j]);
            }
        }

        RistP3AndBytes sum_agg;
        RistScal sum_agg_blinds = sum(agg_blinds);
        res.s6 += measure_time([&](){
            sum_single_thread(sum_agg.elem, ped_agg.elems);
//            sum_agg.fill_bytes();
            sum_agg.elem = sum_agg.elem - sum_agg_blinds * server.predicate.hh[0];
        });

        RistElem sum_elem;
        res.s6 += measure_time([&](){
            bytes_from_p3(sum_elem.element, sum_agg.elem);
        });

        int agg_int;
        res.s6 += measure_time([&](){
            agg_int = discrete_log(sum_elem, server.small_mult_base_table, num_clients);
        });

    }
    res.s6 *= static_cast<double>(pieces);

//    for (int l = 0; l < dim; l++) {
//        auto sum = 0;
//        for (int i = 1; i <= num_clients; i++)
//            sum += client[i].weight_updates[l];
//        assert(server.final_update[l] == sum);
//    }

    return res;
}

void bench_no_mal() {
    int num_clients;
//    std::cout << "number of clients: ";
//    std::cin >> num_clients;

    int max_malicious_clients;
//    std::cout << "max number of malicious clients: ";
//    std::cin >> max_malicious_clients;

    int dim;
//    std::cout << "dimension: ";
//    std::cin >> dim;
//    dim = 1000;

    int num_blinds_per_group_element;


    int num_norm_bound_samples;
//    std::cout << "number of samples in chi-square test: ";
//    std::cin >> num_norm_bound_samples;

    int precomp = 0;
//    std::cout << "precomputation? 1 for yes, 0 for no: ";
//    std::cin >> precomp;
//    precomp = 0;
    if (!precomp) {
        num_blinds_per_group_element = 1;
    }
    else {
        std::cout << "number of private keys: ";
        std::cin >> num_blinds_per_group_element;
//    num_blinds_per_group_element = 10;
    }

    int weight_bits = 24;
    int random_normal_bit_shifter = 24;
    int inner_prod_bound_bits = weight_bits + random_normal_bit_shifter + 4;
    int max_bound_sq_bits = 2 * (weight_bits + random_normal_bit_shifter) + 20;
//    std::vector<std::vector<int>> weight_updates_collection{{0, 0}, // weight_updates_collection[0] is never used
//                                                            {0, 1},
//                                                            {3, 4},
//                                                            {6, -9}};
    RistScal B = power_of_two(max_bound_sq_bits);
    float standard_deviation_factor = 6;

//    std::cout << "number of clients: " << num_clients << std::endl;
//    std::cout << "max maclicious clients: " << max_malicious_clients << std::endl;
//    std::cout << "dimension of weights: " << dim << std::endl;
//    std::cout << "number of private keys to share: " << num_blinds_per_group_element << std::endl;
    int num_pub_keys = (dim % num_blinds_per_group_element == 0) ? dim / num_blinds_per_group_element : dim / num_blinds_per_group_element + 1;
    std::cout << "number of public keys: " << num_pub_keys << std::endl;
//    std::cout << "number of samples in chi-square test: " << num_norm_bound_samples << std::endl;
    std::cout << "gradient bits: " << weight_bits << std::endl;
    std::cout << "normal distribution N(0, 1) samples are multiplied by 2 ^ " << weight_bits << " and rounded to integer" <<std::endl;
    std::cout << "inner prods between grads and samples are bounded by 2 ^ " << inner_prod_bound_bits << std::endl;
    std::cout << "square of l2 bound is at most 2 ^ " << max_bound_sq_bits << std::endl;
    std::cout << "benchmarking ... " << std::endl;


//    auto bench_result = bench_no_mal_1(num_clients, max_malicious_clients, dim, num_blinds_per_group_element,
//                                       num_norm_bound_samples, B,
//                                       weight_bits, random_normal_bit_shifter, inner_prod_bound_bits,
//                                       max_bound_sq_bits,
//                                       static_cast<bool>(precomp));

//    std::cout << "num of norm bound samples:";
//    std::cin >> num_norm_bound_samples;
//    std::cout << num_norm_bound_samples << " samples.....\n";

    dim = 100;
    num_clients = 100;
    max_malicious_clients = 10;
    for (int i = 0; i <= 16; i += 8) {
        for (int sample_mult = 1; sample_mult <= 9; sample_mult *= 3) {
            num_norm_bound_samples = 1000 * sample_mult;
            weight_bits = 16 + i;
            inner_prod_bound_bits = weight_bits + random_normal_bit_shifter + 4;
            max_bound_sq_bits = 2 * (weight_bits + random_normal_bit_shifter) + 20;
            B = power_of_two(max_bound_sq_bits);
            std::cout << "dim: " << dim << std::endl;
            std::cout << "weight bits: " << weight_bits << std::endl;
            std::cout << "inner_prod_bound_bits: " << inner_prod_bound_bits << std::endl;
            std::cout << "max_bound_sq_bits: " << max_bound_sq_bits << std::endl;
            std::cout << "num samples: " << num_norm_bound_samples << std::endl;
            std::cout << "num clients: " << num_clients << std::endl;
            std::cout << "max mal: " << max_malicious_clients << std::endl;

            auto bench_result = bench_no_mal_2(num_clients, max_malicious_clients, dim, num_blinds_per_group_element,
                                               num_norm_bound_samples, B,
                                               weight_bits, random_normal_bit_shifter, inner_prod_bound_bits,
                                               max_bound_sq_bits,
                                               static_cast<bool>(precomp));
            bench_result.print2();
        }
    }

//    for (int dim_mult = 1; dim_mult <= 1000; dim_mult *= 10) {
//        dim = 1000 * dim_mult;
//        num_clients = 100;
//        max_malicious_clients = 10;
//        for (int sample_mult = 1; sample_mult <= 9; sample_mult *= 3) {
//            num_norm_bound_samples = 1000 * sample_mult;
//            std::cout << "dim: " << dim << std::endl;
//            std::cout << "num samples: " << num_norm_bound_samples << std::endl;
//            std::cout << "num clients: " << num_clients << std::endl;
//            std::cout << "max mal: " << max_malicious_clients << std::endl;
//
//            auto bench_result = bench_no_mal_2(num_clients, max_malicious_clients, dim, num_blinds_per_group_element,
//                                               num_norm_bound_samples, B,
//                                               weight_bits, random_normal_bit_shifter, inner_prod_bound_bits,
//                                               max_bound_sq_bits,
//                                               static_cast<bool>(precomp));
//            bench_result.print2();
//        }
//    }

    std::cout << "benchmark complete!\n";//        dim = 100000;
//        for (int num_clients_mult = 1; num_clients_mult <= 5; num_clients_mult++) {
//            num_clients = num_clients_base * num_clients_mult;
//            std::cout << "number of clients:" << num_clients << std::endl;
//            std::cout << "dim:" << dim << std::endl;
//            int count = 1;
//            MicroBench bench;
//            for (int i = 0; i < count; i++) {
//                bench_once(bench, num_clients, dim);
//                std::cout << "bench no. " << i + 1 << " complete!\n";
//            }
//            bench.print();
//        }
}