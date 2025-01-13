#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <immintrin.h>

#include "cantor/cantor_fft.h"
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
    unsigned m = 10; // degree n = 2^m
	unsigned n = (1ULL) << m; // n is even and it must be even for keeping the 32 byte
    printf("m = %u, n = %u\n", m, n);
    // Generating the random polynomial in GF_2^128
	__m128i* fx = random_polynomial_gf2128(n);
    __m128i* evals = cantor_fft_gf2128(fx, n);
    __m128i* evals2 = naive_evaluate(fx, n);
    printf("are equal = %b\n", are_equal_vec_128bit(evals, evals2, n));
}

#define ITERATIONS 10
void cantor_vs_lch(){
    printf("m\tCantor\t\tLCH\n");
    for (unsigned m = 3; m < 22; m++){
	    unsigned n = (1ULL) << m; 
        clock_t start, end;
        double time_cantor = 0, time_lch = 0;

        for (unsigned iter = 0; iter < ITERATIONS; ++iter){
            __m128i* fx1 = random_polynomial_gf2128(n);
            start = clock();
            __m128i* evals2 = lch_fft_gf2128(fx1, n);
            end = clock();
            time_lch += 1000 * ((double) (end - start)) / CLOCKS_PER_SEC;

            start = clock();
            __m128i* evals1 = cantor_fft_gf2128(fx1, n);
            end = clock();
            time_cantor += 1000 * ((double) (end - start)) / CLOCKS_PER_SEC;

            // printf("are equal = %b\n", are_equal_vec_128bit(evals1, evals2, n));
        }   
        // __m128i* fx2 = random_polynomial_gf2128(n);

        printf("%u\t%fms\t%fms\n", m, time_cantor/ITERATIONS, time_lch/ITERATIONS);
    }
}

int main(){
    srand(time(NULL));
    // validate_cantor_basis();
    // test_bitpolymul_lch();
    // test_cantor();
    cantor_vs_lch();
    return 0;
}