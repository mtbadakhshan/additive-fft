// Tiny standalone driver for perf-stat measurements of the LCH radix-2k FFT.
// Usage: ./perf_driver <variant> <m> <iters>
//   variant in {r2, r4, r2k1, r2k2, r2k3, r2k4, r2k5}
//
// We sample a single random input outside the timed region, then run `iters`
// calls to the chosen FFT variant. Wrap with `perf stat -e ...` to read out
// cache-miss / cycles / instruction counts attributable to the FFT loops
// (basis-conversion + butterfly + replication). The same harness is used for
// all variants so per-call overheads are identical across variants.

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <libff/algebra/fields/binary/gf256.hpp>
#include "LCH/fft.hpp"

template<typename FieldT>
static std::vector<FieldT> random_vec(size_t n) {
    std::vector<FieldT> v(n);
    for (size_t i = 0; i < n; ++i) v[i] = FieldT::random_element();
    return v;
}

template<typename FieldT>
static void run(const std::string& which, size_t m, size_t iters,
                const std::vector<FieldT>& input)
{
    std::vector<FieldT> out;
    auto t0 = std::chrono::steady_clock::now();
    for (size_t k = 0; k < iters; ++k) {
        if      (which == "r2")    out = lch::additive_FFT<FieldT>(input, m, m);
        else if (which == "r4")    out = lch::additive_FFT_radix4<FieldT>(input, m, m);
        else if (which == "r2k1")  out = lch::additive_FFT_radix2k<1, FieldT>(input, m, m);
        else if (which == "r2k2")  out = lch::additive_FFT_radix2k<2, FieldT>(input, m, m);
        else if (which == "r2k3")  out = lch::additive_FFT_radix2k<3, FieldT>(input, m, m);
        else if (which == "r2k4")  out = lch::additive_FFT_radix2k<4, FieldT>(input, m, m);
        else if (which == "r2k5")  out = lch::additive_FFT_radix2k<5, FieldT>(input, m, m);
        else {
            std::fprintf(stderr, "unknown variant: %s\n", which.c_str());
            std::exit(2);
        }
        // Defeat dead-code elimination.
        asm volatile("" : : "r"(out.data()) : "memory");
    }
    auto t1 = std::chrono::steady_clock::now();
    double s = std::chrono::duration<double>(t1 - t0).count();
    std::fprintf(stderr,
        "[perf_driver] variant=%s m=%zu iters=%zu total=%.3fs per_call=%.3fms\n",
        which.c_str(), m, iters, s, 1e3 * s / static_cast<double>(iters));
}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::fprintf(stderr, "usage: %s <variant> <m> <iters>\n", argv[0]);
        return 1;
    }
    const std::string which = argv[1];
    const size_t m     = static_cast<size_t>(std::atoi(argv[2]));
    const size_t iters = static_cast<size_t>(std::atoi(argv[3]));

    typedef libff::gf256 FieldT;

    auto input = random_vec<FieldT>(1ull << m);

    run<FieldT>(which, m, iters, input);
    return 0;
}
