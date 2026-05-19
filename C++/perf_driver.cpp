// Standalone driver for `perf stat` measurements of FFT variants.
//
// Usage:
//   ./perf_driver <variant> <m> <iters>
//
// LCH variants:
//   r2, r4, r2k1, r2k2, r2k3, r2k4, r2k5
//   r2_par, r4_par, r2k1_par..r2k5_par  — OpenMP Strategy A (module parallelism)
//
// Cantor (affine subspace) variants:
//   cantor_r2          — cantor::additive_FFT (radix-2)
//   cantor_r2_par      — cantor::additive_FFT_parallel (OpenMP Strategy A)
//   cantor_r2kK        — cantor::additive_FFT_radix2k<K> for K in {2,3,4}
//   cantor_r2kK_par    — cantor::additive_FFT_radix2k_parallel<K>
//
// Example (pin threads and aggregate perf counters):
//   OMP_NUM_THREADS=4 OMP_PLACES=cores OMP_PROC_BIND=close \
//   taskset -c 0-3 perf stat -e cycles,instructions,cache-references,cache-misses -r 5 \
//     ./perf_driver cantor_r2_par 18 40
//
// Compare serial vs parallel at fixed m:
//   perf stat ... ./perf_driver cantor_r2 18 40
//   perf stat ... ./perf_driver cantor_r2_par 18 40

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <libff/algebra/fields/binary/gf256.hpp>

#include "Cantor/cantor_basis.hpp"
#include "Cantor/fft.hpp"
#include "libiop/algebra/subspace.hpp"
#include "LCH/fft.hpp"

template<typename FieldT>
static std::vector<FieldT> random_vec(size_t n)
{
    std::vector<FieldT> v(n);
    for (size_t i = 0; i < n; ++i)
        v[i] = FieldT::random_element();
    return v;
}

template<typename FieldT>
static void run_lch(const std::string &which, size_t m, size_t iters,
                    const std::vector<FieldT> &input)
{
    std::vector<FieldT> out;
    auto t0 = std::chrono::steady_clock::now();
    for (size_t k = 0; k < iters; ++k) {
        if (which == "r2")
            out = lch::additive_FFT<FieldT>(input, m, m);
        else if (which == "r4")
            out = lch::additive_FFT_radix4<FieldT>(input, m, m);
        else if (which == "r2k1")
            out = lch::additive_FFT_radix2k<1, FieldT>(input, m, m);
        else if (which == "r2k2")
            out = lch::additive_FFT_radix2k<2, FieldT>(input, m, m);
        else if (which == "r2k3")
            out = lch::additive_FFT_radix2k<3, FieldT>(input, m, m);
        else if (which == "r2k4")
            out = lch::additive_FFT_radix2k<4, FieldT>(input, m, m);
        else if (which == "r2k5")
            out = lch::additive_FFT_radix2k<5, FieldT>(input, m, m);
        else if (which == "r2_par")
            out = lch::additive_FFT_parallel<FieldT>(input, m, m);
        else if (which == "r4_par")
            out = lch::additive_FFT_radix4_parallel<FieldT>(input, m, m);
        else if (which == "r2k1_par")
            out = lch::additive_FFT_radix2k_parallel<1, FieldT>(input, m, m);
        else if (which == "r2k2_par")
            out = lch::additive_FFT_radix2k_parallel<2, FieldT>(input, m, m);
        else if (which == "r2k3_par")
            out = lch::additive_FFT_radix2k_parallel<3, FieldT>(input, m, m);
        else if (which == "r2k4_par")
            out = lch::additive_FFT_radix2k_parallel<4, FieldT>(input, m, m);
        else if (which == "r2k5_par")
            out = lch::additive_FFT_radix2k_parallel<5, FieldT>(input, m, m);
        else {
            std::fprintf(stderr, "unknown LCH variant: %s\n", which.c_str());
            std::exit(2);
        }
        asm volatile("" : : "r"(out.data()) : "memory");
    }
    auto t1 = std::chrono::steady_clock::now();
    double s = std::chrono::duration<double>(t1 - t0).count();
    std::fprintf(stderr,
                 "[perf_driver] variant=%s m=%zu iters=%zu total=%.3fs per_call=%.3fms\n",
                 which.c_str(), m, iters, s, 1e3 * s / static_cast<double>(iters));
}

template<typename FieldT>
static void run_cantor(const std::string &which, size_t m, size_t iters,
                       const std::vector<FieldT> &poly_coeffs,
                       const libiop::affine_subspace<FieldT> &domain)
{
    std::vector<FieldT> out;
    auto t0 = std::chrono::steady_clock::now();
    for (size_t k = 0; k < iters; ++k) {
        if (which == "cantor_r2")
            out = cantor::additive_FFT<FieldT>(poly_coeffs, domain);
        else if (which == "cantor_r2_par")
            out = cantor::additive_FFT_parallel<FieldT>(poly_coeffs, domain);
        else if (which == "cantor_r2k2")
            out = cantor::additive_FFT_radix2k<2, FieldT>(poly_coeffs, domain);
        else if (which == "cantor_r2k3")
            out = cantor::additive_FFT_radix2k<3, FieldT>(poly_coeffs, domain);
        else if (which == "cantor_r2k4")
            out = cantor::additive_FFT_radix2k<4, FieldT>(poly_coeffs, domain);
        else if (which == "cantor_r2k2_par")
            out = cantor::additive_FFT_radix2k_parallel<2, FieldT>(poly_coeffs, domain);
        else if (which == "cantor_r2k3_par")
            out = cantor::additive_FFT_radix2k_parallel<3, FieldT>(poly_coeffs, domain);
        else if (which == "cantor_r2k4_par")
            out = cantor::additive_FFT_radix2k_parallel<4, FieldT>(poly_coeffs, domain);
        else {
            std::fprintf(stderr, "unknown Cantor variant: %s\n", which.c_str());
            std::exit(2);
        }
        asm volatile("" : : "r"(out.data()) : "memory");
    }
    auto t1 = std::chrono::steady_clock::now();
    double s = std::chrono::duration<double>(t1 - t0).count();
    std::fprintf(stderr,
                 "[perf_driver] variant=%s m=%zu iters=%zu total=%.3fs per_call=%.3fms\n",
                 which.c_str(), m, iters, s, 1e3 * s / static_cast<double>(iters));
}

int main(int argc, char **argv)
{
    if (argc != 4) {
        std::fprintf(stderr,
                     "usage: %s <variant> <m> <iters>\n"
                     "  LCH: r2, r4, r2k1..r2k5, r2_par, r4_par, r2k1_par..r2k5_par\n"
                     "  Cantor: cantor_r2, cantor_r2_par, cantor_r2k2|3|4, cantor_r2k2_par|3_par|4_par\n",
                     argv[0]);
        return 1;
    }
    const std::string which = argv[1];
    const size_t m          = static_cast<size_t>(std::atoi(argv[2]));
    const size_t iters      = static_cast<size_t>(std::atoi(argv[3]));

    typedef libff::gf256 FieldT;

    auto input = random_vec<FieldT>(1ull << m);

    const bool cantor_affine = which == "cantor_r2" || which == "cantor_r2_par" || which == "cantor_r2k2"
                               || which == "cantor_r2k3" || which == "cantor_r2k4" || which == "cantor_r2k2_par"
                               || which == "cantor_r2k3_par" || which == "cantor_r2k4_par";

    if (cantor_affine) {
        std::vector<FieldT> basis(cantor_basis<FieldT>(m));
        libiop::affine_subspace<FieldT> domain(basis, FieldT::random_element());
        run_cantor<FieldT>(which, m, iters, input, domain);
    } else {
        run_lch<FieldT>(which, m, iters, input);
    }

    return 0;
}
