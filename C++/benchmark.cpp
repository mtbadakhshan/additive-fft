// Google Benchmark suite for additive FFT variants (Cantor, LCH, IFFT, radix-2^k).
//
// Build:  cmake --build C++/build --target run_benchmark
// Run:    ./scripts/bench.sh          (see scripts/README.md for deps and env vars)
// Manual: export BM_MIN_RANGE=14 BM_MAX_RANGE=22 BM_STEP=2
//         C++/build/run_benchmark --benchmark_filter='BM_cantor_additive_fft'
//
// Requires BM_MIN_RANGE, BM_MAX_RANGE, BM_STEP at process start (set by bench.sh).
// For cache/IPC counters use perf_driver + scripts/perf_stat.sh (not this binary).

#include <benchmark/benchmark.h>
#include <cstdlib>
#include <sstream>
#include <string>
#ifdef __linux__
#include <unistd.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
#include <libff/algebra/fields/binary/gf256.hpp>
#include "libiop/algebra/subspace.hpp"
#include "libiop/algebra/utils.hpp"
#include "libiop/fft.hpp"
#include "Cantor/fft.hpp"
#include "Gao/fft.hpp"
#include "LCH/fft.hpp"

// Benchmark for libiop::naive_FFT -----------------------------------------------------------------------------------------------------------------------
static void BM_libiop_naive_fft(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain = libiop::field_subset<FieldT>(libiop::affine_subspace<FieldT>::random_affine_subspace(m));
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(domain);
        benchmark::DoNotOptimize(result = libiop::naive_FFT<FieldT>(poly_coeffs, domain));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

// Benchmark for libiop::additive_FFT -----------------------------------------------------------------------------------------------------------------------
static void BM_libiop_additive_fft(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain = libiop::field_subset<FieldT>(libiop::affine_subspace<FieldT>::random_affine_subspace(m));
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(domain);
        benchmark::DoNotOptimize(result = libiop::additive_FFT<FieldT>(poly_coeffs, domain.subspace()));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

// Benchmark for gao::additive_FFT (lvl1) -----------------------------------------------------------------------------------------------------------------------
static void BM_gao_additive_fft_lvl1(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain = libiop::field_subset<FieldT>(libiop::affine_subspace<FieldT>::random_affine_subspace(m));
    gao::PreComputedValues_Level1<FieldT> values_lvl1 = gao::pre_computation_lvl1<FieldT>(domain.subspace());
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(domain);
        benchmark::DoNotOptimize(result = gao::additive_FFT<FieldT>(poly_coeffs, values_lvl1));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

// Benchmark for gao::additive_FFT (lvl2) -----------------------------------------------------------------------------------------------------------------------
static void BM_gao_additive_fft_lvl2(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain = libiop::field_subset<FieldT>(libiop::affine_subspace<FieldT>::random_affine_subspace(m));
    gao::PreComputedValues_Level2<FieldT> values_lvl2 = gao::pre_computation_lvl2<FieldT>(domain.subspace());
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(domain);
        benchmark::DoNotOptimize(result = gao::additive_FFT<FieldT>(poly_coeffs, values_lvl2));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

// Benchmark for gao::additive_FFT_CO -----------------------------------------------------------------------------------------------------------------------
static void BM_gao_additive_fft_co(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    std::vector<FieldT> basis(cantor_basis<FieldT>(m));
    std::reverse(basis.begin(), basis.end());

    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(domain);
        benchmark::DoNotOptimize(result = gao::additive_FFT_CO<FieldT>(poly_coeffs, domain.subspace()));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

// Benchmark for gao::additive_FFT_CO (lvl2) -----------------------------------------------------------------------------------------------------------------------
static void BM_gao_additive_fft_co_lvl2(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    std::vector<FieldT> basis(cantor_basis<FieldT>(m));
    std::reverse(basis.begin(), basis.end());

    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
    gao::PreComputedValues_CO<FieldT> values = gao::pre_computation_co<FieldT>(domain.subspace());
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(domain);
        benchmark::DoNotOptimize(result = gao::additive_FFT_CO<FieldT>(poly_coeffs, values));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

// Benchmark for cantor::additive_FFT -----------------------------------------------------------------------------------------------------------------------
static void BM_cantor_additive_fft(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    std::vector<FieldT> basis(cantor_basis<FieldT>(m));

    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(domain);
        benchmark::DoNotOptimize(result = cantor::additive_FFT<FieldT>(poly_coeffs, domain.subspace()));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

static void BM_cantor_additive_fft_parallel(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    std::vector<FieldT> basis(cantor_basis<FieldT>(m));

    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(domain);
        benchmark::DoNotOptimize(result = cantor::additive_FFT_parallel<FieldT>(poly_coeffs, domain.subspace()));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

// Templated radix-2^K benchmark; instantiated for K = 2, 3, 4.
template<size_t K>
static void BM_cantor_additive_fft_radix2k(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    std::vector<FieldT> basis(cantor_basis<FieldT>(m));

    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(domain);
        benchmark::DoNotOptimize(result = cantor::additive_FFT_radix2k<K, FieldT>(poly_coeffs, domain.subspace()));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

template<size_t K>
static void BM_cantor_additive_fft_radix2k_parallel(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    std::vector<FieldT> basis(cantor_basis<FieldT>(m));

    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(domain);
        benchmark::DoNotOptimize(
            result = cantor::additive_FFT_radix2k_parallel<K, FieldT>(poly_coeffs, domain.subspace()));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

// IFFT benchmarks (affine chart: random evals as input; same pattern as FFT using random coeffs).
static void BM_cantor_additive_ifft(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    std::vector<FieldT> basis(cantor_basis<FieldT>(m));

    const std::vector<FieldT> evals = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(evals);
        benchmark::DoNotOptimize(domain);
        benchmark::DoNotOptimize(result = cantor::additive_IFFT<FieldT>(evals, domain.subspace()));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

static void BM_cantor_additive_ifft_parallel(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    std::vector<FieldT> basis(cantor_basis<FieldT>(m));

    const std::vector<FieldT> evals = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(evals);
        benchmark::DoNotOptimize(domain);
        benchmark::DoNotOptimize(
            result = cantor::additive_IFFT_parallel<FieldT>(evals, domain.subspace()));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

template<size_t K>
static void BM_cantor_additive_ifft_radix2k(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    std::vector<FieldT> basis(cantor_basis<FieldT>(m));

    const std::vector<FieldT> evals = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(evals);
        benchmark::DoNotOptimize(domain);
        benchmark::DoNotOptimize(
            result = cantor::additive_IFFT_radix2k<K, FieldT>(evals, domain.subspace()));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

template<size_t K>
static void BM_cantor_additive_ifft_radix2k_parallel(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    std::vector<FieldT> basis(cantor_basis<FieldT>(m));

    const std::vector<FieldT> evals = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(evals);
        benchmark::DoNotOptimize(domain);
        benchmark::DoNotOptimize(
            result = cantor::additive_IFFT_radix2k_parallel<K, FieldT>(evals, domain.subspace()));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

// Benchmark for cantor::additive_FFT (precmp-basis)-----------------------------------------------------------------------------------------------------------------------
static void BM_cantor_additive_fft_precmp_basis(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);

    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(result = cantor::additive_FFT<FieldT>(poly_coeffs, m, 31));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

// Cantor table path (domain_dim, shift_dim) — OpenMP parallel (s_i / cantor_combinations).
static void BM_cantor_additive_fft_precmp_basis_parallel(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    constexpr size_t shift_dim = 31;

    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(
            result = cantor::additive_FFT_parallel<FieldT>(poly_coeffs, m, shift_dim));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

template<size_t K>
static void BM_cantor_additive_fft_precmp_basis_radix2k(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);

    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(
            result = cantor::additive_FFT_radix2k<K, FieldT>(poly_coeffs, m, 31));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

template<size_t K>
static void BM_cantor_additive_fft_precmp_basis_radix2k_parallel(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);

    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(
            result = cantor::additive_FFT_radix2k_parallel<K, FieldT>(poly_coeffs, m, 31));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

// Benchmark for cantor::additive_FFT (precmp)-----------------------------------------------------------------------------------------------------------------------
static void BM_cantor_additive_fft_precmp(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    std::vector<FieldT> basis(cantor_basis<FieldT>(m));

    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
    cantor::PreComputedValues<FieldT> values = cantor::pre_computation(domain.subspace());
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(domain);
        benchmark::DoNotOptimize(result = cantor::additive_FFT<FieldT>(poly_coeffs, values));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

// Benchmark for lch::additive_FFT (precmp-basis)-----------------------------------------------------------------------------------------------------------------------
static void BM_lch_additive_fft_precmp_basis(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);

    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(result = lch::additive_FFT<FieldT>(poly_coeffs, m, 31));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

// Benchmark for lch::additive_FFT_radix2k<K> (precmp-basis) ---------------
template<size_t K>
static void BM_lch_additive_fft_radix2k_precmp_basis(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);

    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(result = lch::additive_FFT_radix2k<K, FieldT>(poly_coeffs, m, 31));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

static void BM_lch_additive_fft_parallel_precmp_basis(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);

    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(result = lch::additive_FFT_parallel<FieldT>(poly_coeffs, m, 31));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

template<size_t K>
static void BM_lch_additive_fft_radix2k_parallel_precmp_basis(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);

    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(
            result = lch::additive_FFT_radix2k_parallel<K, FieldT>(poly_coeffs, m, 31));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

// IFFT benchmarks (precmp-basis: random evals as input).
static void BM_cantor_additive_ifft_precmp_basis(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    constexpr size_t shift_dim = 31;

    const std::vector<FieldT> evals = libiop::random_vector<FieldT>(1ull << m);
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(evals);
        benchmark::DoNotOptimize(result = cantor::additive_IFFT<FieldT>(evals, m, shift_dim));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

static void BM_cantor_additive_ifft_precmp_basis_parallel(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    constexpr size_t shift_dim = 31;

    const std::vector<FieldT> evals = libiop::random_vector<FieldT>(1ull << m);
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(evals);
        benchmark::DoNotOptimize(result = cantor::additive_IFFT_parallel<FieldT>(evals, m, shift_dim));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

template<size_t K>
static void BM_cantor_additive_ifft_precmp_basis_radix2k(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    constexpr size_t shift_dim = 31;

    const std::vector<FieldT> evals = libiop::random_vector<FieldT>(1ull << m);
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(evals);
        benchmark::DoNotOptimize(
            result = cantor::additive_IFFT_radix2k<K, FieldT>(evals, m, shift_dim));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

template<size_t K>
static void BM_cantor_additive_ifft_precmp_basis_radix2k_parallel(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    constexpr size_t shift_dim = 31;

    const std::vector<FieldT> evals = libiop::random_vector<FieldT>(1ull << m);
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(evals);
        benchmark::DoNotOptimize(
            result = cantor::additive_IFFT_radix2k_parallel<K, FieldT>(evals, m, shift_dim));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

static void BM_lch_additive_ifft_precmp_basis(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    constexpr size_t shift_dim = 31;

    const std::vector<FieldT> evals = libiop::random_vector<FieldT>(1ull << m);
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(evals);
        benchmark::DoNotOptimize(result = lch::additive_IFFT<FieldT>(evals, m, shift_dim));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

static void BM_lch_additive_ifft_parallel_precmp_basis(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    constexpr size_t shift_dim = 31;

    const std::vector<FieldT> evals = libiop::random_vector<FieldT>(1ull << m);
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(evals);
        benchmark::DoNotOptimize(result = lch::additive_IFFT_parallel<FieldT>(evals, m, shift_dim));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

template<size_t K>
static void BM_lch_additive_ifft_radix2k_precmp_basis(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    constexpr size_t shift_dim = 31;

    const std::vector<FieldT> evals = libiop::random_vector<FieldT>(1ull << m);
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(evals);
        benchmark::DoNotOptimize(result = lch::additive_IFFT_radix2k<K, FieldT>(evals, m, shift_dim));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

template<size_t K>
static void BM_lch_additive_ifft_radix2k_parallel_precmp_basis(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);
    constexpr size_t shift_dim = 31;

    const std::vector<FieldT> evals = libiop::random_vector<FieldT>(1ull << m);
    std::vector<FieldT> result;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(evals);
        benchmark::DoNotOptimize(
            result = lch::additive_IFFT_radix2k_parallel<K, FieldT>(evals, m, shift_dim));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}


const int MIN_RANGE = std::stoi(std::getenv("BM_MIN_RANGE"));
const int MAX_RANGE = std::stoi(std::getenv("BM_MAX_RANGE"));
const int STEP = std::stoi(std::getenv("BM_STEP"));

#define REGISTER_BENCH(fn) \
    BENCHMARK(fn)->DenseRange(MIN_RANGE, MAX_RANGE, STEP)->Unit(benchmark::kMicrosecond)->ReportAggregatesOnly(true)
#define REGISTER_BENCH_K(fn, k) \
    BENCHMARK_TEMPLATE(fn, k)->DenseRange(MIN_RANGE, MAX_RANGE, STEP)->Unit(benchmark::kMicrosecond)->ReportAggregatesOnly(true)

// =============================================================================
// PAPER SUITE (active) — only the forward-FFT table-path (precmp-basis) families
// behind the paper's runtime tables (see scripts/measure.conf.paper). These use
// the *_parallel* routines swept over OMP_NUM_THREADS, where the P=1 pass is the
// serial baseline. To restore the full suite (affine chart, IFFT, and serial
// baselines), flip the `#if 0` block below to `#if 1`.
// =============================================================================

// Cantor table-path forward FFT — OpenMP-free serial baselines
REGISTER_BENCH(BM_cantor_additive_fft_precmp_basis);
REGISTER_BENCH_K(BM_cantor_additive_fft_precmp_basis_radix2k, 2);
REGISTER_BENCH_K(BM_cantor_additive_fft_precmp_basis_radix2k, 3);
REGISTER_BENCH_K(BM_cantor_additive_fft_precmp_basis_radix2k, 4);

// Cantor table-path forward FFT (radix-2 + radix-2^K, parallel)
REGISTER_BENCH(BM_cantor_additive_fft_precmp_basis_parallel);
REGISTER_BENCH_K(BM_cantor_additive_fft_precmp_basis_radix2k_parallel, 2);
REGISTER_BENCH_K(BM_cantor_additive_fft_precmp_basis_radix2k_parallel, 3);
REGISTER_BENCH_K(BM_cantor_additive_fft_precmp_basis_radix2k_parallel, 4);

// LCH table-path forward FFT — OpenMP-free serial baselines
REGISTER_BENCH(BM_lch_additive_fft_precmp_basis);
REGISTER_BENCH_K(BM_lch_additive_fft_radix2k_precmp_basis, 2);
REGISTER_BENCH_K(BM_lch_additive_fft_radix2k_precmp_basis, 3);
REGISTER_BENCH_K(BM_lch_additive_fft_radix2k_precmp_basis, 4);
REGISTER_BENCH_K(BM_lch_additive_fft_radix2k_precmp_basis, 5);

// LCH table-path forward FFT (radix-2 + radix-2^K, parallel)
REGISTER_BENCH(BM_lch_additive_fft_parallel_precmp_basis);
REGISTER_BENCH_K(BM_lch_additive_fft_radix2k_parallel_precmp_basis, 2);
REGISTER_BENCH_K(BM_lch_additive_fft_radix2k_parallel_precmp_basis, 3);
REGISTER_BENCH_K(BM_lch_additive_fft_radix2k_parallel_precmp_basis, 4);
REGISTER_BENCH_K(BM_lch_additive_fft_radix2k_parallel_precmp_basis, 5);

#if 0  // ---- Full suite (affine chart + IFFT): flip to `#if 1` to re-enable ----

// Cantor affine-chart FFT
REGISTER_BENCH(BM_cantor_additive_fft);
REGISTER_BENCH(BM_cantor_additive_fft_parallel);
REGISTER_BENCH_K(BM_cantor_additive_fft_radix2k, 2);
REGISTER_BENCH_K(BM_cantor_additive_fft_radix2k, 3);
REGISTER_BENCH_K(BM_cantor_additive_fft_radix2k, 4);
REGISTER_BENCH_K(BM_cantor_additive_fft_radix2k_parallel, 2);
REGISTER_BENCH_K(BM_cantor_additive_fft_radix2k_parallel, 3);
REGISTER_BENCH_K(BM_cantor_additive_fft_radix2k_parallel, 4);

// Cantor affine-chart IFFT
REGISTER_BENCH(BM_cantor_additive_ifft);
REGISTER_BENCH(BM_cantor_additive_ifft_parallel);
REGISTER_BENCH_K(BM_cantor_additive_ifft_radix2k, 2);
REGISTER_BENCH_K(BM_cantor_additive_ifft_radix2k, 3);
REGISTER_BENCH_K(BM_cantor_additive_ifft_radix2k, 4);
REGISTER_BENCH_K(BM_cantor_additive_ifft_radix2k_parallel, 2);
REGISTER_BENCH_K(BM_cantor_additive_ifft_radix2k_parallel, 3);
REGISTER_BENCH_K(BM_cantor_additive_ifft_radix2k_parallel, 4);

// Cantor table-path inverse FFT (parallel + serial)
REGISTER_BENCH(BM_cantor_additive_ifft_precmp_basis_parallel);
REGISTER_BENCH_K(BM_cantor_additive_ifft_precmp_basis_radix2k_parallel, 2);
REGISTER_BENCH_K(BM_cantor_additive_ifft_precmp_basis_radix2k_parallel, 3);
REGISTER_BENCH_K(BM_cantor_additive_ifft_precmp_basis_radix2k_parallel, 4);
REGISTER_BENCH(BM_cantor_additive_ifft_precmp_basis);
REGISTER_BENCH_K(BM_cantor_additive_ifft_precmp_basis_radix2k, 2);
REGISTER_BENCH_K(BM_cantor_additive_ifft_precmp_basis_radix2k, 3);
REGISTER_BENCH_K(BM_cantor_additive_ifft_precmp_basis_radix2k, 4);

// LCH table-path inverse FFT (parallel + serial)
REGISTER_BENCH(BM_lch_additive_ifft_parallel_precmp_basis);
REGISTER_BENCH_K(BM_lch_additive_ifft_radix2k_parallel_precmp_basis, 2);
REGISTER_BENCH_K(BM_lch_additive_ifft_radix2k_parallel_precmp_basis, 3);
REGISTER_BENCH_K(BM_lch_additive_ifft_radix2k_parallel_precmp_basis, 4);
REGISTER_BENCH_K(BM_lch_additive_ifft_radix2k_parallel_precmp_basis, 5);
REGISTER_BENCH(BM_lch_additive_ifft_precmp_basis);
REGISTER_BENCH_K(BM_lch_additive_ifft_radix2k_precmp_basis, 2);
REGISTER_BENCH_K(BM_lch_additive_ifft_radix2k_precmp_basis, 3);
REGISTER_BENCH_K(BM_lch_additive_ifft_radix2k_precmp_basis, 4);
REGISTER_BENCH_K(BM_lch_additive_ifft_radix2k_precmp_basis, 5);
#endif

#undef REGISTER_BENCH
#undef REGISTER_BENCH_K

// Other families (libiop, Gao, PreComputedValues Cantor path) — enable manually when needed.
// BENCHMARK(BM_libiop_additive_fft)->DenseRange(MIN_RANGE, MAX_RANGE, STEP)->Unit(benchmark::kMicrosecond)->ReportAggregatesOnly(true);
// BENCHMARK(BM_cantor_additive_fft_precmp)->DenseRange(MIN_RANGE, MAX_RANGE, STEP)->Unit(benchmark::kMicrosecond)->ReportAggregatesOnly(true);


int main(int argc, char **argv)
{
    char arg0_default[] = "benchmark";
    char *args_default = arg0_default;
    if (!argv) {
        argc = 1;
        argv = &args_default;
    }

    const char *omp_env = std::getenv("OMP_NUM_THREADS");
    benchmark::AddCustomContext("OMP_NUM_THREADS_env", omp_env ? omp_env : "(unset)");
#ifdef __linux__
    const long nproc = sysconf(_SC_NPROCESSORS_ONLN);
    benchmark::AddCustomContext("processors_online",
                                nproc > 0 ? std::to_string(nproc) : "(unknown)");
#else
    benchmark::AddCustomContext("processors_online", "(sysconf unavailable)");
#endif
#ifdef _OPENMP
    benchmark::AddCustomContext("openmp_max_threads", std::to_string(omp_get_max_threads()));
#else
    benchmark::AddCustomContext("openmp_max_threads", "1 (no OpenMP)");
#endif

    benchmark::Initialize(&argc, argv);
    if (benchmark::ReportUnrecognizedArguments(argc, argv))
        return 1;
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
