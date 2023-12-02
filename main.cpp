#include <iostream>
#include "test/test.cpp"
#include "bench/bench_no_mal.cpp"
#include "bench/bench_random.cpp"
#include "bench/bench_inner_prod_mod_p.cpp"
#include "bench/bench_precomp.cpp"
#include "bench/bench_batch_add.cpp"

int main() {

    test();
//    bench_random();
//    bench_inner_prod_mod_p();
//bench_precomp();
//    bench_batch_add();
//    bench_no_mal();
    return 0;
}
