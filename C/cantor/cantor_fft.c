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
    for (unsigned r = 0; r < m; ++r){
        // Computing S_i
        S_index--; 
        // t=0;
        // if (S_index > 0){
        //     nz_S[t] = (1<<S_index) - 1; t++;
        //     for ( unsigned i = 1; i < S_index; ++i){
        //         ii = i; rr = S_index;
        //         while (((rr & 1) | (~ii & 1)) && ii>0) {ii >>= 1; rr >>= 1;}
        //         if (ii == 0) {nz_S[t] = (1<<S_index) - (1<<i); t++;}
        //     }
        // }

        nz_S = s_i[S_index]; t=n_terms[S_index];

        // printf("{%u", nz_S[0]); for(unsigned i = 1; i < t; ++i) printf(", %u", nz_S[i]); printf("}\t\t\t\t\t \\\\S_%u\n", S_index);

        offset = 0;
        half_input_size = input_size >> 1;
        half_half_input_size = half_input_size>>1;
        for (unsigned module = 0; module < n_modules; ++module){
            mult_factor = bitmat_prod_accu_64x128_M8R_sse(_mm_setzero_si128(),  gfCantorto2128_8R, module << 1);
            offset2 = offset + half_input_size;
            for (unsigned k = offset2 + half_input_size - 1; k >= offset2; --k){
                poly_k = poly[k];
                for (unsigned i = 0; i < t; ++i){
                    poly[k - nz_S[i]] ^= poly_k;
                }
                poly[k-half_input_size] ^= _gf2ext128_mul_sse(poly_k, mult_factor);
            }
            j = 0;
            poly256 =(__m256i*) (poly+offset);
            for (; j < half_half_input_size; ++j) poly256[j + half_half_input_size] ^= poly256[j];
            j*=2;
            for (; j < half_input_size; ++j) poly[offset2+j] ^= poly[offset+j];            
            offset += input_size;
        }
        input_size = half_input_size;
        n_modules <<= 1;
    }
    return poly;
}
