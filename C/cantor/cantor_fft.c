#include "cantor_fft.h"
#include <stdio.h>
#include <string.h>    // Header for memcpy

#include "s_i.h"
#include "utils/utils.h"
#include "bitpolymul/bitmat_prod.h"
#include "bitpolymul/gf2128_cantor_iso.h"
#include "bitpolymul/gfext_aesni.h"



#define LOG2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))

__m128i* cantor_fft_gf2128(__m128i* fx, unsigned n_term){
    #ifdef NCOPY_POLY //For small polynomials copy to a new array to keep the input values unchanged
    __m128i* poly = fx;
    #else
    __m128i* poly = (__m128i*)aligned_alloc( 32 , sizeof(__m128i)*n_term );
    poly = memcpy( poly, fx, sizeof(__m128i)*n_term );
    #endif
    unsigned m = LOG2(n_term);
    unsigned S_index = m, input_size = n_term, n_modules = 1;
    unsigned j, t, ii, rr, offset, offset2, half_input_size, half_half_input_size;
    __m128i mult_factor, poly_k;
    __m256i* poly256;
    // unsigned nz_S[15];
    const unsigned* nz_S;
    for (unsigned r = 0; r < m-1; ++r){
        // Computing S_i
        S_index--; 
        nz_S = s_i[S_index]; t=n_terms[S_index];
        offset = 0;
        half_input_size = input_size >> 1;
        half_half_input_size = half_input_size>>1;
        for (unsigned module = 0; module < n_modules; ++module){
            mult_factor = bitmat_prod_accu_64x128_M8R_sse(_mm_setzero_si128(),  gfCantorto2128_8R, module<<1);
            offset2 = offset + half_input_size;
            for (unsigned k = offset2 + half_input_size - 1; k >= offset2; --k){
                poly_k = poly[k];
                for (unsigned i = 0; i < t; ++i){
                    poly[k - nz_S[i]] ^= poly_k;
                }
                // poly[k-half_input_size] ^= _gf2ext128_mul_sse(poly_k, mult_factor);
            }
            poly256 =(__m256i*) (poly+offset);
            for (j = 0; j < half_half_input_size; ++j) { // we use half_half_input_size since the steps are 256-bits 
                poly256[j + half_half_input_size] = poly256[j] ^= _gf2ext128_mul_2x1_avx2( poly256[j + half_half_input_size] , mult_factor );
                // poly256[j + half_half_input_size] ^= poly256[j];  
            }      
            offset += input_size;
        }
        input_size = half_input_size;
        n_modules <<= 1;
    }
    
    offset = 0;
    // r = m (last layer)
    for (unsigned module = 0; module < n_modules; ++module){
        mult_factor = bitmat_prod_accu_64x128_M8R_sse(_mm_setzero_si128(),  gfCantorto2128_8R, module<<1);
        offset2 = offset + half_input_size;
        poly[offset] ^= _gf2ext128_mul_sse(poly[offset+1], mult_factor);
        poly[offset+1] ^= poly[offset];            
        offset += 2;
    }

    return poly;
}
