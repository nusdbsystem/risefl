//
// Created by yizheng on 26/3/23.
//

#include <vector>
#include <algorithm>
#include <execution>

#include <thread>

#include "include/sodium/ed25519_ref10.h"
#include "include/ristretto.h"
#include "include/rist_fast_computation.h"
#include "include/zkp.h"


bool RistElemP3::operator==(const RistElemP3 &other) const {
    // different p3 forms could represent the same Ristretto group element
    // need to convert to bytes to check equality
    RistElem a, b;
    bytes_from_p3(a.element, other);
    bytes_from_p3(b.element, *this);
    return a == b;
}

void swap(RistElemP3 &p, RistElemP3 &q) {
    std::swap(p.element.X, q.element.X);
    std::swap(p.element.Y, q.element.Y);
    std::swap(p.element.Z, q.element.Z);
    std::swap(p.element.T, q.element.T);
}

RistElemP3 &RistElemP3::operator=(RistElemP3 other) {
    swap(*this, other);
    return *this;
}

RistElemP3::RistElemP3(const RistElemP3 &other) {
    p3_copy(*this, other);
}

RistElemP3::RistElemP3(RistElemP3 &&other) noexcept {
    swap(*this, other);
}


void p3_copy(RistElemP3 &r, const RistElemP3 &p) {
    fe25519_copy(r.element.X, p.element.X);
    fe25519_copy(r.element.Y, p.element.Y);
    fe25519_copy(r.element.Z, p.element.Z);
    fe25519_copy(r.element.T, p.element.T);
}
//void p3_copy(RistElemP3Vec& r, const RistElemP3Vec& p) {
//    for (int i = 0; i < p.size(); i++){
//        p3_copy(r[i], p[i]);
//    }
//}

void bytestream_from_p3_vec(std::vector<unsigned char>::iterator &it_bytes, const RistElemP3Vec &pp) {
    std::for_each(
            std::execution::par_unseq,
            pp.begin(),
            pp.end(),
            [&pp, &it_bytes](auto &&p) {
                auto i = &p - pp.data();
                auto it = it_bytes + i * RISTBYTES;
                bytestream_from_p3(it, p);
            }
    );
    it_bytes += pp.size() * RISTBYTES;
//    for (int i = 0; i < pp.size(); i++) {
//        bytestream_from_p3(it_bytes, pp[i]);
//    }
}

std::vector<unsigned char> bytes_from_p3_vec(const RistElemP3Vec &pp) {
    std::vector<unsigned char> bytes(pp.size() * RISTBYTES);
    auto it = bytes.begin();
    bytestream_from_p3_vec(it, pp);
    return bytes;
}

void scalar_mult(RistElemP3Vec &rr, const RistScal &n, const RistElemP3Vec &pp) {
    std::for_each(
            std::execution::par_unseq,
            rr.begin(),
            rr.end(),
            [&rr, &n, &pp](auto &&r) {
                auto i = &r - rr.data();
                scalar_mult(rr[i], n, pp[i]);
            }
    );
}


void pedersen_commit(RistElemP3 &y, const RistScal &x, const RistElemP3 &h, const RistScal &r) {
    RistElemP2 p2;
    ge25519_double_scalarmult_vartime(&p2.element, r.scalar.data(), &h.element, x.scalar.data());
    p2_to_p3(y, p2);
}

void pedersen_commit(RistElemP3 &y, long x, int x_bits, const RistElemP3 &h, const RistScal &r) {
    RistElemP2 p2;
    ge25519_double_scalarmult_vartime_int(&p2.element, r.scalar.data(), &h.element, x, x_bits);
    p2_to_p3(y, p2);
}


void pedersen_commit(RistElemP3Vec &yy, const RistScalVec &xx, const RistElemP3Vec &hh, const RistScal &r) {
    std::for_each(
            std::execution::par_unseq,
            yy.begin(),
            yy.end(),
            [&yy, &xx, &hh, &r](auto &&y) {
                auto i = &y - yy.data();
                pedersen_commit(y, xx[i], hh[i], r);
            }
    );
}

void pedersen_commit(RistElemP3Vec &yy, const RistScalVec &xx, const RistElemP3 &h, const RistScalVec &rr) {
    std::for_each(
            std::execution::par_unseq,
            yy.begin(),
            yy.end(),
            [&yy, &xx, &h, &rr](auto &&y) {
                auto i = &y - yy.data();
                pedersen_commit(y, xx[i], h, rr[i]);
            }
    );
}

void pedersen_commit_precomp(RistElemP3 &y, long x, int x_bits, const RistElemPrecompTable &table, const RistScal &r) {
    scalar_double_mult_precomp(y, x, x_bits, r, table);
}

void pedersen_commit_precomp(RistElemP3 &y, const RistScal &x, const RistElemPrecompTable &table, const RistScal &r) {
    RistElemP3 from_table;
    scalar_mult_precomp(from_table, r, table);
    y = pedersen_zero_commit_p3(x) + from_table;
}


void pedersen_zero_commit(RistElemP3Vec &yy, const RistScalVec &xx) {
    std::for_each(
            std::execution::par_unseq,
            yy.begin(),
            yy.end(),
            [&yy, &xx](auto &&y) {
                auto i = &y - yy.data();
                pedersen_zero_commit(y, xx[i]);
            }
    );
}
//std::vector<unsigned char> pedersen_commit_p3_to_bytes_from_precomp(const RistScalVec& xx, const RistElemP3Vec& hh, const RistScal& r) {
//    std::vector<unsigned char> bytes(xx.size() * RISTBYTES);
//    pedersen_commit_p3_to_bytes_from_precomp(bytes, xx, hh, r);
//    return bytes;
//}

void pedersen_commit_p3_to_bytes(std::vector<unsigned char>::iterator &it_bytes, const RistScalVec &xx,
                                 const RistElemP3Vec &hh, const RistScal &r) {
    std::for_each(
            std::execution::par_unseq,
            xx.begin(),
            xx.end(),
            [&it_bytes, &xx, &hh, &r](auto &&x) {
                auto i = &x - xx.data();
                RistElemP3 y;
                pedersen_commit(y, xx[i], hh[i], r);
                ristretto255_p3_tobytes(&*it_bytes + i * RISTBYTES, &y.element);
            }
    );
    it_bytes += xx.size() * RISTBYTES;
}

void
pedersen_commit_p3_to_bytes(std::vector<unsigned char>::iterator &it_bytes, const std::vector<long> &xx, int x_bits,
                            const RistElemP3Vec &hh, const RistScal &r) {
    std::for_each(
            std::execution::par_unseq,
            xx.begin(),
            xx.end(),
            [&it_bytes, &xx, &hh, &r, &x_bits](auto &&x) {
                auto i = &x - xx.data();
                RistElemP3 y;
                pedersen_commit(y, xx[i], x_bits, hh[i], r);
                ristretto255_p3_tobytes(&*it_bytes + i * RISTBYTES, &y.element);
            }
    );
    it_bytes += xx.size() * RISTBYTES;
}


void
pedersen_commit_p3_to_bytes_from_precomp(std::vector<unsigned char>::iterator &it_bytes, const std::vector<long> &xx,
                                         int xx_bits,
                                         const std::vector<RistElemPrecompTable> &tables,
                                         const RistScalVec &rr) {
    std::for_each(
            std::execution::par_unseq,
            tables.begin(),
            tables.end(),
            [&it_bytes, &xx, &xx_bits, &tables, &rr](auto &&table) {
                auto i = &table - tables.data();
                for (int j = i * rr.size(); j < (i + 1) * rr.size() && j < xx.size(); j++) {
                    RistElemP3 y;
                    pedersen_commit_precomp(y, xx[j], xx_bits, table, rr[j - i * rr.size()]);
                    ristretto255_p3_tobytes(&*it_bytes + j * RISTBYTES, &y.element);
                }
            }
    );
    it_bytes += xx.size() * RISTBYTES;
}


void sum_single_thread(RistElemP3 &s, const RistElemP3Vec &vv) {
    p3_0(s);
    for (auto &&v: vv) {
        p3_add(s, s, v);
    }
}

void sum_single_thread(RistElemP3 &s, const RistElemP3Mat &vv) {
    p3_0(s);
    for (auto &&v: vv) {
        sum_single_thread(s, v);
    }
}


//
//RistElemP3 sum_parallel(const RistElemP3Vec& vv) {
//    return std::reduce(std::execution::par_unseq,std::begin(vv), std::end(vv));
//}


unsigned long get_index(const RistScal &r, int start, int b) {
    unsigned long index = 0;
    for (int i = 0; i < b; i++) {
        if (b_get_bit(r, start + i)) {
            index += (1L << i);
        }
    }
    return index;
}

int get_index(int r, int start, int b) {
    if (r >= 0) {
        return (r >> start) & ((1 << b) - 1);
    } else {
        return -get_index(-r, start, b);
    }
}

void compute_accumulation(RistElemP3Vec &accumulation, const RistScalVec &rr, const RistElemP3Vec &pp,
                          int start, int b) {
    for (int j = 0; j < accumulation.size(); j++)
        p3_0(accumulation[j]);
    for (int i = 0; i < rr.size(); i++) {
        auto index = get_index(rr[i], start, b);
        if (index > 0)
            accumulation[index - 1] = accumulation[index - 1] + pp[i];
    }
}

void compute_accumulation(RistElemP3Vec &accumulation_pos, const std::vector<int> &rr, const RistElemP3Vec &pp,
                          int start, int b) {
    for (int j = 0; j < accumulation_pos.size(); j++) {
        p3_0(accumulation_pos[j]);
    }
    for (int i = 0; i < rr.size(); i++) {
        auto index = get_index(rr[i], start, b);
        if (index > 0)
            accumulation_pos[index - 1] = accumulation_pos[index - 1] + pp[i];
        else if (index < 0) {
            int index_abs = -index;
            accumulation_pos[index_abs - 1] = accumulation_pos[index_abs - 1] - pp[i];
        }
    }
}

void compute_step_sum(RistElemP3 &step_sum, const RistElemP3Vec &acc) {
    RistElemP3 sum;
    auto size = acc.size();
    sum = acc[size - 1];
    step_sum = acc[size - 1];
    for (int i = size - 2; i >= 0; i--) {
        sum += acc[i];
        step_sum += sum;
    }
}

void linear_comb(RistElemP3 &s, const RistScalVec &rr, const RistElemP3Vec &pp) {
    assert(rr.size() == pp.size());

    int n = rr.size();
    if (n == 0) {
        p3_0(s);
        return;
    }
    if (n == 1) {
        scalar_mult(s, rr[0], pp[0]);
        return;
    }

    int b = get_power_two_bound(RistScal(n)) / 2 + 1;
    // Home pc experiments suggest that splitting bits into intervals of size b = log(n) / 2 + 1 is fastest

    int num_intervals = (255 % b == 0) ? 255 / b : 255 / b + 1;
    std::vector<RistElemP3Vec> accumulations(num_intervals, RistElemP3Vec((1 << b) - 1));
    RistElemP3Vec prods(num_intervals);

    std::for_each(
            std::execution::par_unseq,
            accumulations.begin(),
            accumulations.end(),
            [&accumulations, &prods, &rr, &pp, &b](auto &&accumulation) {
                auto i = &accumulation - accumulations.data();
                compute_accumulation(accumulation, rr, pp, b * i, b);
                compute_step_sum(prods[i], accumulation);
            });

    s = prods[num_intervals - 1];
    for (int i = num_intervals - 2; i >= 0; i--) {
        power_2_add(s, b, prods[i]);
    }
}

void linear_comb(RistElemP3 &s, std::vector<int> &rr, int bit_bound, const RistElemP3Vec &pp) {
    assert(rr.size() == pp.size());
    assert(bit_bound < 32);

    int n = rr.size();
    if (n == 1) {
        scalar_mult(s, RistScal(rr[0]), pp[0]);
        return;
    }

    int b = get_power_two_bound(RistScal(n)) / 2;

    int num_intervals = (bit_bound % b == 0) ? bit_bound / b : bit_bound / b + 1;
    std::vector<RistElemP3Vec> accumulations_pos(num_intervals, RistElemP3Vec((1 << b) - 1));
    RistElemP3Vec prods(num_intervals);

    std::for_each(
            std::execution::par_unseq,
            accumulations_pos.begin(),
            accumulations_pos.end(),
            [&accumulations_pos, &prods, &rr, &pp, &b](auto &&accumulation_pos) {
                auto i = &accumulation_pos - accumulations_pos.data();
                compute_accumulation(accumulation_pos, rr, pp, b * i, b);
                compute_step_sum(prods[i], accumulation_pos);
            });

    s = prods[num_intervals - 1];
    for (int i = num_intervals - 2; i >= 0; i--) {
        power_2_add(s, b, prods[i]);
    }
}


void linear_comb(RistElemP3 &s, const RistScalMat &rr, const RistElemP3Mat &pp) {
    assert(rr.size() == pp.size());
    assert(rr[0].size() == pp[0].size());
    linear_comb(s, concatenate_double_vec(rr), concatenate_double_vec(pp));
}


void generate_precomp_table(RistElemPrecompTable &table, const RistElemP3 &r) {
    RistScal n(256);
    RistElemP3 s(r), t(r);
    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 8; j++) {
            ge25519_p3_to_precomp(&table[i][j], &s.element);
            p3_add(s, s, t);
        }
        scalar_mult(t, n, t);
        s = t;
    }
}

int p3_from_bytes(RistElemP3 &p, const std::array<unsigned char, RISTBYTES> &s) {
    return ristretto255_frombytes(&p.element, s.data());
}

int p3_from_bytes(RistElemP3 &p, const std::vector<unsigned char> &s) {
    return ristretto255_frombytes(&p.element, s.data());
}

int p3_from_bytestream(RistElemP3 &p, std::vector<unsigned char>::const_iterator &it_bytes) {
    int ret = ristretto255_frombytes(&p.element, &*it_bytes);
    it_bytes += RISTBYTES;
    return ret;
}

void
p3_vec_from_bytestream(RistElemP3Vec &rr, std::vector<unsigned char>::const_iterator &it_bytes) {
    std::for_each(
            std::execution::par_unseq,
            rr.begin(),
            rr.end(),
            [&rr, &it_bytes](auto &&r) {
                auto i = &r - rr.data();
                ristretto255_frombytes(&r.element, &*it_bytes + i * RISTBYTES);
            }
    );
    it_bytes += rr.size() * RISTBYTES;
}

void p3_vec_from_bytes(RistElemP3Vec &rr, const std::vector<unsigned char> &pp) {
    auto it = pp.cbegin();
    p3_vec_from_bytestream(rr, it);
}

void bytes_from_p3(std::array<unsigned char, RISTBYTES> &s, const RistElemP3 &p) {
    ristretto255_p3_tobytes(s.data(), &p.element);
}

void bytes_from_p3(std::vector<unsigned char> &s, const RistElemP3 &p) {
    ristretto255_p3_tobytes(s.data(), &p.element);
}


void bytestream_from_p3(std::vector<unsigned char>::iterator &it_bytes, const RistElemP3 &p) {
    ristretto255_p3_tobytes(&*it_bytes, &p.element);
    it_bytes += RISTBYTES;
}

std::array<unsigned char, RISTBYTES> bytes_from_p3(const RistElemP3 &p) {
    std::array<unsigned char, RISTBYTES> bytes;
    ristretto255_p3_tobytes(bytes.data(), &p.element);
    return bytes;
}

void rist_elem_vec_from_p3_vec(RistElemVec &rr, RistElemP3Vec &pp) {
    std::for_each(
            std::execution::par_unseq,
            rr.begin(),
            rr.end(),
            [&rr, &pp](auto &&r) {
                auto i = &r - rr.data();
                bytes_from_p3(r.element, pp[i]);
            }
    );
}

RistElemP3 pedersen_commit(const RistScal &x, const RistElemP3 &h, const RistScal &r) {
    RistElemP3 y;
    pedersen_commit(y, x, h, r);
    return y;
}

void pedersen_zero_commit(RistElemP3 &y, const RistScal &x) {
    ge25519_scalarmult_base(&y.element, x.scalar.data());
}

void pedersen_zero_commit(RistElemP3 &y, long x) {
    pedersen_zero_commit(y, RistScal(x));
}

RistElemP3 pedersen_zero_commit_p3(const RistScal &x) {
    RistElemP3 y;
    pedersen_zero_commit(y, x);
    return y;
}

RistElemP3 pedersen_zero_commit_p3(long x) {
    RistElemP3 y;
    pedersen_zero_commit(y, x);
    return y;
}

RistElemP3Vec pedersen_commit(const RistScalVec &xx, const RistElemP3Vec &hh, const RistScal &r) {
    RistElemP3Vec yy(xx.size());
    pedersen_commit(yy, xx, hh, r);
    return yy;
}

void scalar_mult(RistElemP3 &r, const RistScal &n, const RistElemP3 &p) {
    ge25519_scalarmult(&r.element, n.scalar.data(), &p.element);
}

RistElemP3 operator*(const RistScal &n, const RistElemP3 &p) {
    RistElemP3 r;
    scalar_mult(r, n, p);
    return r;
}

RistElemP3Vec operator*(const RistScal &n, const RistElemP3Vec &pp) {
    RistElemP3Vec rr(pp.size());
    scalar_mult(rr, n, pp);
    return rr;
}

void p3_add(RistElemP3Vec &rr, const RistElemP3Vec &pp, const RistElemP3Vec &qq) {
    std::for_each(
            std::execution::par_unseq,
            rr.begin(),
            rr.end(),
            [&rr, &pp, &qq](auto &&r) {
                auto i = &r - rr.data();
                p3_add(rr[i], pp[i], qq[i]);
            }
    );
}

RistElemP3Vec operator+(const RistElemP3Vec &pp, const RistElemP3Vec &qq) {
    RistElemP3Vec rr(pp.size());
    p3_add(rr, pp, qq);
    return rr;
}

RistElemP3 linear_comb(const RistScalVec &rr, const RistElemP3Vec &pp) {
    RistElemP3 s;
    linear_comb(s, rr, pp);
    return s;
}

RistElemP3 linear_comb(const RistScalMat &rr, const RistElemP3Mat &pp) {
    RistElemP3 s;
    linear_comb(s, rr, pp);
    return s;
}

// ss[0] = rr[0] * pp[0] + rr[n] * pp[1] + rr[2 * n] * pp[2] + ...
// ss[1] = rr[1] * pp[0] + rr[n + 1] * pp[1] + rr[2 * n + 1] * pp[2] + ...
// ...
// ss[n - 1] = rr[n - 1] * pp[0] + rr[2 * n - 1] * pp[1] + rr[3 * n - 1] * pp[2] + ...
void linear_comb_block(RistElemP3Vec &ss, const RistScalVec &rr, const RistElemP3Vec &pp) {
    int dim = rr.size();
    int n = ss.size();
    int m = pp.size();
    std::for_each(std::execution::par_unseq, ss.begin(), ss.end(),
                  [&ss, &rr, &pp, &dim, &m, &n](auto &&s) {
                      auto i = &s - ss.data();
                      RistScalVec rr_thread(m);
                      int j = 0;
                      for (; j * n + i < dim; j++) {
                          rr_thread[j] = rr[j * n + i];
                      }
                      for (; j < m; j++) {
                          rr_thread[j] = c_scal_zero;
                      }
                      linear_comb(s, rr_thread, pp);
                  });
}

void linear_comb_block(RistElemP3Vec &ss, const std::vector<int> &rr, int bit_bound, const RistElemP3Vec &pp) {
    int dim = rr.size();
    int n = ss.size();
    int m = pp.size();
    std::for_each(std::execution::par_unseq, ss.begin(), ss.end(),
                  [&ss, &rr, &bit_bound, &pp, &dim, &m, &n](auto &&s) {
                      auto i = &s - ss.data();
                      std::vector<int> rr_thread(m);
                      int j = 0;
                      for (; j * n + i < dim; j++) {
                          rr_thread[j] = rr[j * n + i];
                      }
                      for (; j < m; j++) {
                          rr_thread[j] = 0;
                      }
                      linear_comb(s, rr_thread, bit_bound, pp);
                  });
}

void p3_from_hash(RistElemP3 &r, const RistHashbytes &h) {
    p3_from_hash(&r.element, h.hashbytes.data());
}

RistElemP3 p3_from_hash(const RistHashbytes &h) {
    RistElemP3 r;
    p3_from_hash(r, h);
    return r;
}

void scalar_mult_precomp(RistElemP3 &r, const RistScal &n, const RistElemPrecompTable &table) {
    ge25519_scalarmult_precomp(&r.element, n.scalar.data(), table.data());
}


void
scalar_double_mult_precomp(RistElemP3 &r, long x, int x_bits, const RistScal &n, const RistElemPrecompTable &table) {
    ge25519_double_scalarmult_precomp(&r.element, x, x_bits, n.scalar.data(), table.data());
}

RistP3AndBytes zero_rist_p3_bytes() {
    RistP3AndBytes ret;
    p3_0(ret.elem);
    ret.fill_bytes();
    return ret;
}

RistP3AndBytes one_rist_p3_bytes() {
    RistP3AndBytes ret;
    pedersen_zero_commit(ret.elem, 1);
    ret.fill_bytes();
    return ret;
}

RistElemP3 linear_comb_trim(const RistScalVec &rr, const RistElemP3Vec &ss) {
    int n = rr.size();
    int m = ss.size();
    assert(n >= m);
    return linear_comb(RistScalVec(rr.begin(), rr.begin() + m), ss);
}

RistElemP3Vec operator*(const RistScalVec &nn, const RistElemP3Vec &pp) {
    assert(nn.size() == pp.size());
    RistElemP3Vec ret(nn.size());
    std::for_each(std::execution::par_unseq, ret.begin(), ret.end(), [&ret, &nn, &pp](auto &&r) {
        auto i = &r - ret.data();
        ret[i] = nn[i] * pp[i];
    });
    return ret;
}

RistElemP3 sum_single_thread(const RistElemP3Vec &vv) {
    RistElemP3 s;
    sum_single_thread(s, vv);
    return s;
}

void rand_init(RistP3AndBytes &h) {
    RistHashbytes hash;
    rand_init(hash);
    p3_from_hash(h.elem, hash);
    h.fill_bytes();
}