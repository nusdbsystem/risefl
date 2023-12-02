//
// Created by yizheng on 4/4/23.
//

#include "include/mpz_conversion.h"

void MyMpz::reduce_to_rist() {
    mpz_mod(pointer, pointer, c_mpz_rist_order.pointer);
}

void mpz_from_rist(MyMpz &x, const RistScal &s) {
    mpz_import(x.pointer, 32, -1, 1, 1, 0, s.scalar.data());
}

void mpz_from_int(MyMpz &x, const int n) {
    mpz_set_si(x.pointer, n);
}

void rist_from_mpz(RistScal &s, const MyMpz &x) {
    memset(s.scalar.data(), 0, s.scalar.size());
    MyMpz temp(x);
    temp.reduce_to_rist();
    mpz_export(s.scalar.data(), NULL, -1, 1, 1, 0, temp.pointer);
//    if (mpz_sgn(x.pointer) >= 0) {
//        mpz_export(s.scalar.data(), NULL, -1, 1, 1, 0, x.pointer);
//    }
//    else {
//        mpz_export(s.scalar.data(), NULL, -1, 1, 1, 0, x.pointer);
//        s = -s;
//    }
}


void mpz_vec_from_int(MpzVec &xx, const std::vector<int> &nn) {
    for (int i = 0; i < xx.size(); i++) {
        mpz_from_int(xx[i], nn[i]);
    }
}

//void mpz_mat_from_int(MpzMat& xx, const std::vector<std::vector<int>> &nn) {
//    for (int i = 0; i < xx.size(); i++) {
//        mpz_vec_from_int(xx[i], nn[i]);
//    }
//}


void mpz_vec_from_rist(MpzVec &xx, const RistScalVec &ss) {
    std::for_each(
            std::execution::par_unseq,
            xx.begin(),
            xx.end(),
            [&ss, &xx](auto &&x) {
                auto i = &x - xx.data();
                mpz_from_rist(xx[i], ss[i]);
            }
    );
}

void mpz_mat_from_rist(MpzMat &xx, const RistScalMat &ss) {
    std::for_each(
            std::execution::par_unseq,
            xx.begin(),
            xx.end(),
            [&ss, &xx](auto &&x) {
                auto i = &x - xx.data();
                mpz_vec_from_rist(xx[i], ss[i]);
            }
    );
}

void rist_from_mpz_vec(RistScalVec &ss, const MpzVec &xx) {
    std::for_each(
            std::execution::par_unseq,
            xx.begin(),
            xx.end(),
            [&ss, &xx](auto &&x) {
                auto i = &x - xx.data();
                rist_from_mpz(ss[i], xx[i]);
            }
    );
}

void rist_from_mpz_mat(RistScalMat &ss, const MpzMat &xx) {
    std::for_each(
            std::execution::par_unseq,
            xx.begin(),
            xx.end(),
            [&ss, &xx](auto &&x) {
                auto i = &x - xx.data();
                rist_from_mpz_vec(ss[i], xx[i]);
            }
    );
}

