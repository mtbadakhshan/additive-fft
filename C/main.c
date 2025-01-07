#include <stdio.h>
#include <stdlib.h>
#include <emmintrin.h>
#include <immintrin.h>
#include "utils/lch_api.h"
#include "utils/utils.h"
#include "bitpolymul/bc.h"
#include "bitpolymul/butterfly_net.h"


#define LOG2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))

__m128i* random_polynomial_gf2128(unsigned n){
    __m128i* fx = (__m128i*)aligned_alloc( 32 , sizeof(__m128i)*n );
	if( NULL == fx ) { printf("alloc fail.\n"); exit(-1); }
	for (unsigned i = 0; i < n; ++i) fx[i] = _mm_set_epi32(rand(), rand(), rand(), rand());
    return fx;
} 

void test_bitpolymul_lch(){
	unsigned m = 3; // degree n = 2^m
	unsigned n = (1ULL) << m; // n is even and it must be even for keeping the 32 byte

    // Generating the random polynomial in GF_2^128
	__m128i* fx = random_polynomial_gf2128(n);
    print_vector_128bit_name(fx, n, "fx");
    __m128i* evals = lch_fft_gf2128(fx, n);
    print_vector_128bit_name(fx, n, "fx");

    
    // _bc_to_lch_128( fx , n , 1 );
    // butterfly_net_clmul( fx ,n );


}

int main(){
    test_bitpolymul_lch();
    return 0;
}