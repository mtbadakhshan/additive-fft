#ifndef _LCH_API_H_
#define _LCH_API_H_
#include <immintrin.h>
    /*
    The bitpolymul doesn't have butterfly_net_clmul function, since it is for multiplying polynomials. 
    Therefore, they have butterfly_net_half_inp_clmul instead.
    */
    // butterfly_net_clmul( __m128i *poly , unsigned n_fx ); 
    __m128i* lch_fft_gf2128(__m128i* poly, unsigned n_term);
#endif
