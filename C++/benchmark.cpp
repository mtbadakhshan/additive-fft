#include <benchmark/benchmark.h>
#include <libff/algebra/fields/binary/gf256.hpp>
#include "libiop/algebra/subspace.hpp"
#include "libiop/algebra/utils.hpp"
#include "libiop/fft.hpp"

//Benchmark for libiop::additive_fft
static void BM_libiop_additive_fft(benchmark::State &state){
    typedef libff::gf256 FieldT;
    const size_t m = state.range(0);

    for (auto _ : state)
    {
        state.PauseTiming();
        std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
        libiop::field_subset<FieldT> domain = libiop::field_subset<FieldT>(libiop::affine_subspace<FieldT>::random_affine_subspace(m));
        state.ResumeTiming();

        benchmark::DoNotOptimize(poly_coeffs);
        benchmark::DoNotOptimize(domain);

        const std::vector<FieldT> result = libiop::additive_FFT<FieldT>(poly_coeffs, domain.subspace());
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

BENCHMARK(BM_libiop_additive_fft)->Range(4, 15)->Unit(benchmark::kMicrosecond);

BENCHMARK_MAIN();
