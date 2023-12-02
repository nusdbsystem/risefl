////
//// Created by yizheng on 14/3/23.
////
#include <cassert>
#include "include/ristretto_vector.h"

std::vector<unsigned char> bytes_from_rist_elem_vec(const RistElemVec &gg) {
    std::vector<unsigned char> gg_bytes(gg.size() * RISTBYTES);
    for (int i = 0; i < gg.size(); i++) {
        memcpy(gg_bytes.data() + i * RISTBYTES, gg[i].element.data(),
               RISTBYTES);
    }
    return gg_bytes;
}

std::vector<unsigned char> bytes_from_rist_scal_vec(const RistScalVec &gg) {
    std::vector<unsigned char> gg_bytes(gg.size() * RISTSCALBYTES);
    for (int i = 0; i < gg.size(); i++) {
        memcpy(gg_bytes.data() + i * RISTSCALBYTES, gg[i].scalar.data(),
               RISTSCALBYTES);
    }
    return gg_bytes;
}

RistElemVec rist_elem_vec_from_bytes(const std::vector<unsigned char> &vv) {
    RistElemVec gg(vv.size() / RISTBYTES);
    for (int i = 0; i < gg.size(); i++) {
        memcpy(gg[i].element.data(), vv.data() + i * RISTBYTES, RISTBYTES);
    }
    return gg;
}

RistScalVec rist_scal_vec_from_bytes(const std::vector<unsigned char> &vv) {
    RistScalVec gg(vv.size() / RISTSCALBYTES);
    for (int i = 0; i < gg.size(); i++) {
        memcpy(gg[i].scalar.data(), vv.data() + i * RISTSCALBYTES,
               RISTSCALBYTES);
    }
    return gg;
}

RistScalVec rist_scal_vec_from_int_vec(const std::vector<int> &s) {
    RistScalVec ret(s.size());
    for (int i = 0; i < s.size(); i++) {
        ret[i] = RistScal(s[i]);
    }
    return ret;
}

RistScalVec rist_scal_vec_from_long_vec(const std::vector<long> &s) {
    RistScalVec ret(s.size());
    for (int i = 0; i < s.size(); i++) {
        ret[i] = RistScal(s[i]);
    }
    return ret;
}

void export_to_bytestream(const RistScalVec &rr, std::vector<unsigned char>::iterator &it) {
    for (auto &&r: rr)
        r.export_to_bytestream(it);
}

void import_from_bytestream(RistScalVec &rr, std::vector<unsigned char>::const_iterator &it) {
    for (auto &&r: rr)
        r.import_from_bytestream(it);
}

void export_to_bytestream(const RistScalMat &rr, std::vector<unsigned char>::iterator &it) {
    for (auto &&r: rr)
        export_to_bytestream(r, it);
}

void import_from_bytestream(RistScalMat &rr, std::vector<unsigned char>::const_iterator &it) {
    for (auto &&r: rr)
        import_from_bytestream(r, it);
}

void rand_init(RistElemVec &pp) {
    std::for_each(
            std::execution::par_unseq,
            pp.begin(),
            pp.end(),
            [](auto &&p) {
                rand_init(p);
            }
    );
}

void hash_init(RistElemVec &pp, const RistHashbytesVec &rr) {
    std::for_each(
            std::execution::par_unseq,
            pp.begin(),
            pp.end(),
            [&pp, &rr](auto &&p) {
                auto i = &p - pp.data();
                hash_init(p, rr[i]);
            }
    );
}

void add(RistElemVec &rr, const RistElemVec &pp, const RistElemVec &qq) {
    std::for_each(
            std::execution::par_unseq,
            rr.begin(),
            rr.end(),
            [&rr, &pp, &qq](auto &&r) {
                auto i = &r - rr.data();
                add(rr[i], pp[i], qq[i]);
            }
    );
}

void scalar_mult(RistElemVec &qq, const RistScalVec &nn, const RistElemVec pp) {
    assert(nn.size() == pp.size());
    std::for_each(
            std::execution::par_unseq,
            qq.begin(),
            qq.end(),
            [&qq, &nn, &pp](auto &&q) {
                auto i = &q - qq.data();
                scalar_mult(qq[i], nn[i], pp[i]);
            }
    );
}

RistElemVec operator*(const RistScalVec &nn, const RistElemVec pp) {
    RistElemVec qq(nn.size());
    scalar_mult(qq, nn, pp);
    return qq;
}

void scalar_mult_base(RistElemVec &qq, const RistScalVec &nn) {
    std::for_each(
            std::execution::par_unseq,
            qq.begin(),
            qq.end(),
            [&qq, &nn](auto &&q) {
                auto i = &q - qq.data();
                scalar_mult_base(qq[i], nn[i]);
            }
    );
}

void rand_init(RistScalVec &pp) {
    std::for_each(
            std::execution::par_unseq,
            pp.begin(),
            pp.end(),
            [](auto &&p) {
                rand_init(p);
            }
    );
}

void scalar_reduce(RistScalVec &rr, const RistNonreducedScalarVec &ss) {
    std::for_each(
            std::execution::par_unseq,
            rr.begin(),
            rr.end(),
            [&rr, &ss](auto &&r) {
                auto i = &r - rr.data();
                scalar_reduce(rr[i], ss[i]);
            }
    );
}

void scalar_reduce(RistScalVec &rr, const RistHashbytesVec &ss) {
    scalar_reduce(rr, *reinterpret_cast<const RistNonreducedScalarVec *>(&ss));
}

void scalar_add(RistScalVec &zz, const RistScalVec &xx, const RistScalVec &yy) {
    assert(xx.size() == yy.size());
    std::for_each(
            std::execution::par_unseq,
            zz.begin(),
            zz.end(),
            [&zz, &xx, &yy](auto &&z) {
                auto i = &z - zz.data();
                scalar_add(zz[i], xx[i], yy[i]);
            }
    );
}

RistScalVec operator+(const RistScalVec &xx, const RistScalVec &yy) {
    RistScalVec zz(xx.size());
    scalar_add(zz, xx, yy);
    return zz;
}

void scalar_sub(RistScalVec &zz, const RistScalVec &xx, const RistScalVec &yy) {
    std::for_each(
            std::execution::par_unseq,
            zz.begin(),
            zz.end(),
            [&zz, &xx, &yy](auto &&z) {
                auto i = &z - zz.data();
                scalar_sub(zz[i], xx[i], yy[i]);
            }
    );
}

RistScalVec operator-(const RistScalVec &xx, const RistScalVec &yy) {
    RistScalVec zz(xx.size());
    scalar_sub(zz, xx, yy);
    return zz;
}


RistScalVec operator-(const RistScalVec &xx) {
    RistScalVec zz(xx.size());
    for (int i = 0; i < xx.size(); i++) {
        zz[i] = -xx[i];
    }
    return zz;
}

void scalar_mul(RistScalVec &zz, const RistScalVec &xx, const RistScalVec &yy) {
    assert(xx.size() == yy.size());
    std::for_each(
            std::execution::par_unseq,
            zz.begin(),
            zz.end(),
            [&zz, &xx, &yy](auto &&z) {
                auto i = &z - zz.data();
                scalar_mul(zz[i], xx[i], yy[i]);
            }
    );
}

RistScalVec operator*(const RistScalVec &xx, const RistScalVec &yy) {
    RistScalVec zz{xx.size()};
    scalar_mul(zz, xx, yy);
    return zz;
}

RistScalVec operator*(const RistScal &x, const RistScalVec &yy) {
    RistScalVec zz{yy.size()};
    std::for_each(
            std::execution::par_unseq,
            zz.begin(),
            zz.end(),
            [&zz, &x, &yy](auto &&z) {
                auto i = &z - zz.data();
                scalar_mul(zz[i], x, yy[i]);
            }
    );
    return zz;
}

RistElem sum(const RistElemVec &pp) {
    return std::reduce(std::execution::par_unseq,
                       pp.begin(), pp.end());
}

RistScal sum(const RistScalVec &pp) {
    return std::reduce(std::execution::par_unseq,
                       pp.begin(), pp.end());
}

RistElem linear_comb(const RistScalVec &nn, const RistElemVec &pp) {
    return sum(nn * pp);
}


RistScalMat operator*(const RistScal &r, const RistScalMat &B) {
    RistScalMat ret;
    ret.reserve(B.size());
    for (int i = 0; i < B.size(); i++) {
        ret.emplace_back(r * B[i]);
    }
    return ret;
}


RistScalMat operator+(const RistScalMat &A, const RistScalMat &B) {
    assert(A.size() == B.size());
    RistScalMat ret;
    ret.reserve(B.size());
    for (int i = 0; i < B.size(); i++) {
        ret.emplace_back(A[i] + B[i]);
    }
    return ret;
}

RistScalMat operator-(const RistScalMat &A) {
    RistScalMat ret;
    ret.reserve(A.size());
    for (int i = 0; i < A.size(); i++) {
        ret.emplace_back(-A[i]);
    }
    return ret;
}
