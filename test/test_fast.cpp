//
// Created by yizheng on 27/3/23.
//
#include "../include/rist_fast_computation.h"
#include <chrono>
#include <cassert>

void test_linear_comb();

void test_sum() {
    int length = 1000;
    RistElem a1, a2;
    RistElemP3 b1, b2;
    rand_init(a1);
    p3_from_bytes(b1, a1.element);
    bytes_from_p3(a2.element, b1);
    assert(a1 == a2);
    p3_from_bytes(b2, a2.element);
//    assert(b1.T == b2.T);

    RistElemVec hh_bytes(length);
    rand_init(hh_bytes);
    RistElemP3Vec hh_p3(length);

    typedef std::chrono::high_resolution_clock Clock;

    auto t1 = Clock::now();
    for (int i = 0; i < length; i++) {
        p3_from_bytes(hh_p3[i], hh_bytes[i].element);
    }
    auto t2 = Clock::now();

    long double mytime;
    mytime = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / double(1000000000);
    std::cout << length << " conversions bytes ==> p3 single thread took " << mytime << " seconds" <<std::endl;


    t1 = Clock::now();
    for (int i = 0; i < length; i++) {
        bytes_from_p3(hh_bytes[i].element, hh_p3[i]);
    }
    t2 = Clock::now();

    mytime = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / double(1000000000);
    std::cout << length << " conversions p3 ==> bytes single thread took " << mytime << " seconds" <<std::endl;

    auto a_sum = sum(hh_bytes);

    RistElemP3 p_sum;
    sum_single_thread(p_sum, hh_p3);

    RistElem a_sum_from_p;
    bytes_from_p3(a_sum_from_p.element, p_sum);

    assert (a_sum_from_p == a_sum);


    RistScalVec xx(length), rr(length);
    rand_init(xx);
    rand_init(rr);
    RistElemP3Vec yy_p3(length);

    t1 = Clock::now();
    for (int i = 0; i < length; i++) {
        pedersen_zero_commit(yy_p3[i], xx[i]);
    }
    t2 = Clock::now();

    mytime = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / double(1000000000);
    std::cout << length << " pedersen zero commits in P3 single thread took " << mytime << " seconds" <<std::endl;

    RistElemVec yy_bytes(length);
    t1 = Clock::now();
    for (int i = 0; i < length; i++) {
        pedersen_zero_commit(yy_bytes[i], xx[i]);
    }
    t2 = Clock::now();

    mytime = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / double(1000000000);
    std::cout << length << " pedersen zero commits in bytes single thread took " << mytime << " seconds" <<std::endl;

    RistElemVec yy_converted(length);
    for (int i = 0; i < length; i++) {
        bytes_from_p3(yy_converted[i].element, yy_p3[i]);
    }
    assert(yy_converted == yy_bytes);



    t1 = Clock::now();
    for (int i = 0; i < length; i++) {
        pedersen_commit(yy_p3[i], xx[i], hh_p3[i], rr[0]);
    }
    t2 = Clock::now();

    mytime = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / double(1000000000);
    std::cout << length << " pedersen commits in P3 single thread took " << mytime << " seconds" <<std::endl;

    t1 = Clock::now();
    for (int i = 0; i < length; i++) {
        pedersen_commit(yy_bytes[i], xx[i], hh_bytes[i], rr[0]);
    }
    t2 = Clock::now();

    mytime = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / double(1000000000);
    std::cout << length << " pedersen commits in bytes single thread took " << mytime << " seconds" <<std::endl;

    for (int i = 0; i < length; i++) {
        bytes_from_p3(yy_converted[i].element, yy_p3[i]);
    }
    assert(yy_converted == yy_bytes);


    RistElemP3Vec yy_p3_parallel(length);
    t1 = Clock::now();
    pedersen_commit(yy_p3_parallel, xx, hh_p3, rr[0]);
    t2 = Clock::now();

    mytime = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / double(1000000000);
    std::cout << length << " pedersen commits in P3 parallel took " << mytime << " seconds" <<std::endl;
    assert(yy_p3_parallel == yy_p3);

    t1 = Clock::now();
    pedersen_commit(yy_bytes, xx, hh_bytes, rr[0]);
    t2 = Clock::now();

    mytime = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / double(1000000000);
    std::cout << length << " pedersen commits in bytes parallel took " << mytime << " seconds" <<std::endl;

    for (int i = 0; i < length; i++) {
        bytes_from_p3(yy_converted[i].element, yy_p3_parallel[i]);
    }
    assert(yy_converted == yy_bytes);


    t1 = Clock::now();
    pedersen_zero_commit(yy_p3_parallel, xx);
    t2 = Clock::now();

    mytime = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / double(1000000000);
    std::cout << length << " pedersen zero commits in P3 parallel took " << mytime << " seconds" <<std::endl;
//    assert(yy_p3_parallel == yy_p3);

    t1 = Clock::now();
    pedersen_zero_commit(yy_bytes, xx);
    t2 = Clock::now();

    mytime = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / double(1000000000);
    std::cout << length << " pedersen zero commits in bytes parallel took " << mytime << " seconds" <<std::endl;

    for (int i = 0; i < length; i++) {
        bytes_from_p3(yy_converted[i].element, yy_p3_parallel[i]);
    }
    assert(yy_converted == yy_bytes);


    RistElemP2 y_p2, y_p2_another;
    p3_to_p2(y_p2, yy_p3[0]);
    p2_to_p3(yy_p3[1], y_p2);
    p3_to_p2(y_p2_another, yy_p3[1]);
//    assert(y_p2 == y_p2_another);
    assert(yy_p3[0] == yy_p3[1]);


    int n = 10;
    rand_init(a2);
    p3_from_bytes(b2, a2.element);
    auto a3 = RistScal(1L << n) * a1 + a2;
    power_2_add(b1, n, b2);
    RistElem a4;
    bytes_from_p3(a4.element, b1);
    assert(a3 == a4);

//    {
//        RistElemP3Vec aa(length);
//        std::vector<RistHashbytes> hh(length);
//        for (auto && h : hh)
//            rand_init(h);
//        for (int i = 0; i < aa.size(); i++){
//            aa[i] = p3_from_hash(hh[i]);
//        }
//        RistElemP3 s;
//        sum_single_thread(s, aa);
//        assert(sum_parallel(aa) == s);
//
//        std::cout << "sum parallel test passed!\n";
//    }
}

void test_fast() {

    test_sum();

    test_linear_comb();
    std::cout << "test fast complete!\n";
}

void test_linear_comb() {
    int length = 1000;
    RistScalVec qq(length);
    rand_init(qq);
    RistElemP3 comb;
    RistElemP3Vec yy_p3(length);
    RistElemVec yy_bytes(length);
    rand_init(yy_bytes);
    for (int i = 0; i < length; i++) {
        p3_from_bytes(yy_p3[i], yy_bytes[i].element);
    }

    auto t1 = Clock::now();
    linear_comb(comb, qq, yy_p3);
    auto t2 = Clock::now();
    auto mytime = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / double(1000000000);
    std::cout << length << " linear combs in P3 in parallel took " << mytime << " seconds" <<std::endl;

    RistElem comb_conv;
    bytes_from_p3(comb_conv.element, comb);

    for (int i = 0; i < length; i++) {
        bytes_from_p3(yy_bytes[i].element, yy_p3[i]);
    }

    t1 = Clock::now();
    RistElem comb_bytes = linear_comb(qq, yy_bytes);
    t2 = Clock::now();
    mytime = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / double(1000000000);
    std::cout << length << " linear combs bytes in parallel took " << mytime << " seconds" <<std::endl;

    assert(comb_bytes == comb_conv);

}