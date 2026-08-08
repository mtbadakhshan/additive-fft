#ifndef _DYADIC_FFT_H_
#define _DYADIC_FFT_H_
#include <immintrin.h>
    // Recursive variant
    __m128i* dyadic_fft_gf2128(__m128i* poly, unsigned n_term);

    // Iterative variant
    __m128i* dyadic_fft_gf2128_iter(__m128i* poly, unsigned n_term);

    // Cache chunk B (16-byte elements) for the blocked schedule.
    // Priority: (1) dyadic_fft_set_cache_block override, (2) value from
    // dyadic_tune_params.h after `make tune`, (3) round_pow2(L2/16),
    // (4) fallback 2^16. get returns B; cache_block_source returns a short
    // label ("forced", "tuned", "L2", "fallback"). Not thread-safe w.r.t.
    // concurrent FFT calls.
    void dyadic_fft_set_cache_block(unsigned n_elem);
    unsigned dyadic_fft_get_cache_block(void);
    unsigned dyadic_fft_cache_block_for(unsigned n_term);
    const char* dyadic_fft_cache_block_source(void);

#endif
