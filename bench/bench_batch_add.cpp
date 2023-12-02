//
// Created by yizheng on 17/4/23.
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

void bench_batch_add() {
    int length = 1000;
    int num_clients = 100;
    RistElemP3Mat aa(num_clients, RistElemP3Vec(length));
    RistElemP3Vec sum(length);
    RistHashbytes hash;
    for (int j = 0; j < num_clients; j++)
        for (int i = 0; i < length; i++) {
            rand_init(hash);
            p3_from_hash(aa[j][i], hash);
    }
    for (int i = 0; i < length; i++) {
        p3_0(sum[i]);
    }
    auto seq_time = measure_time([&](){
        for (int j = 0; j < num_clients; j++)
            for (int i = 0; i < length; i++) {
                sum[i] += aa[j][i];
            }
    });
    std::cout << num_clients << " clients, " << length << " dim, sequential time: " << seq_time << " seconds.\n";

    RistElemP3Vec sum_par(length);
    for (int i = 0; i < length; i++) {
        p3_0(sum_par[i]);
    }

    auto par_time = measure_time([&](){
            sum_par = std::reduce(std::execution::par_unseq, aa.begin(), aa.end(), sum_par, [&](const RistElemP3Vec &vv1, const RistElemP3Vec &vv2){
            RistElemP3Vec res(vv1.size());
            for (auto i = 0; i < vv1.size(); i++){
//                vv1[i] += vv2[i];
                res[i] = vv1[i] + vv2[i];
            }
            return res;
        });
    });
    std::cout << "par time:" << par_time << " seconds.\n";

    for (int i = 0; i < length; i++) {
        assert(sum[i] == sum_par[i]);
    }

}