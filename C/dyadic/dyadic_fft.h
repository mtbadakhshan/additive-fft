#ifndef _DYADIC_FFT_H_
#define _DYADIC_FFT_H_
#include <immintrin.h>
    // Recursive variant
    __m128i* dyadic_fft_gf2128(__m128i* poly, unsigned n_term);

    // Iterative variant
    __m128i* dyadic_fft_gf2128_iter(__m128i* poly, unsigned n_term);

    // Cache chunk control (sizes in 16-byte elements). By default each FFT
    // call auto-selects an L2-scale chunk when the array fits in L3, and an
    // L3-scale chunk when it spills to DRAM. The setter forces one size for
    // all subsequent calls (rounded down to a power of two; n_elem = 0
    // restores auto-selection). get returns the forced size, or the L2-scale
    // auto default; cache_block_for(n) returns the size a call of length n
    // would use. Not thread-safe w.r.t. concurrent FFT calls.
    void dyadic_fft_set_cache_block(unsigned n_elem);
    unsigned dyadic_fft_get_cache_block(void);
    unsigned dyadic_fft_cache_block_for(unsigned n_term);

#endif
