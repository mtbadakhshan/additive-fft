// Sweeps the dyadic FFT cache-chunk size on the current machine and prints
// mean timings, with LCH as a baseline. Use the best column to pick a value
// for dyadic_fft_set_cache_block() (or to adjust the auto-detection rule).
//
// Build:  make tune.out
// Run:    ./tune.out
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <immintrin.h>

#include "dyadic/dyadic_fft.h"
#include "utils/lch_api.h"
#include "utils/utils.h"

static double now_ms(void){
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1e6;
}

static double bench_dyadic(unsigned n, unsigned iters){
    double total = 0;
    for (unsigned it = 0; it < iters; ++it){
        __m128i* fx = random_polynomial_gf2128(n);
        double t0 = now_ms();
        __m128i* evals = dyadic_fft_gf2128(fx, n);
        total += now_ms() - t0;
        free(fx);
        free(evals);
    }
    return total / iters;
}

static double bench_lch(unsigned n, unsigned iters){
    double total = 0;
    for (unsigned it = 0; it < iters; ++it){
        __m128i* fx = random_polynomial_gf2128(n);
        double t0 = now_ms();
        __m128i* evals = lch_fft_gf2128(fx, n);
        total += now_ms() - t0;
        free(fx);
        free(evals);
    }
    return total / iters;
}

int main(){
    srand(time(NULL));
    const unsigned block_log[] = {13, 14, 15, 16, 17, 18, 19, 20, 21};
    const unsigned n_blocks = sizeof(block_log) / sizeof(block_log[0]);

    dyadic_fft_set_cache_block(0);  // ensure auto mode
    unsigned l2 = dyadic_fft_get_cache_block();
    unsigned big = dyadic_fft_cache_block_for(1u << 24);
    printf("auto L2-scale chunk: %u elements (%u KB)\n", l2, l2 / 64);
    printf("auto L3-scale chunk: %u elements (%u KB)  [used when n > that size]\n",
           big, big / 64);
    printf("(setter overrides both for the sweep below)\n\n");

    printf("m\tLCH");
    for (unsigned b = 0; b < n_blocks; ++b) printf("\t2^%u", block_log[b]);
    printf("\t(ms, mean; best dyadic chunk marked *)\n");

    for (unsigned m = 17; m <= 24; ++m){
        unsigned n = 1u << m;
        unsigned iters = (m <= 19) ? 20 : (m <= 22 ? 8 : 4);

        double lch = bench_lch(n, iters);
        double dy[sizeof(block_log) / sizeof(block_log[0])];
        unsigned best = 0;
        for (unsigned b = 0; b < n_blocks; ++b){
            dyadic_fft_set_cache_block(1u << block_log[b]);
            dy[b] = bench_dyadic(n, iters);
            if (dy[b] < dy[best]) best = b;
        }

        printf("%u\t%.2f", m, lch);
        for (unsigned b = 0; b < n_blocks; ++b)
            printf("\t%.2f%s", dy[b], b == best ? "*" : "");
        printf("\n");
    }

    dyadic_fft_set_cache_block(0);  // restore auto-detection
    return 0;
}
