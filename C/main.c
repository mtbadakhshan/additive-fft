#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <immintrin.h>

#include "utils/lch_api.h"
#include "utils/utils.h"
#include "bitpolymul/bc.h"
#include "bitpolymul/butterfly_net.h"


#define LOG2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))

void test_bitpolymul_lch(){
	unsigned m = 5; // degree n = 2^m
	unsigned n = (1ULL) << m; // n is even and it must be even for keeping the 32 byte
    printf("m = %u, n = %u\n", m, n);
    // Generating the random polynomial in GF_2^128
	__m128i* fx = random_polynomial_gf2128(n);
    __m128i* evals = lch_fft_gf2128(fx, n);
    __m128i* evals2 = naive_evaluate(fx, n);
    printf("are equal = %b\n", are_equal_vec_128bit(evals, evals2, n));
}

void test_cantor(){
    
}

int main(){
    srand(time(NULL));
    // validate_cantor_basis();
    test_bitpolymul_lch();
    return 0;
}