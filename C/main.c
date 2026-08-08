#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <immintrin.h>

#include "cantor/cantor_fft.h"
#include "dyadic/dyadic_fft.h"
#include "utils/lch_api.h"
#include "utils/utils.h"
#include "bitpolymul/bc.h"
#include "bitpolymul/butterfly_net.h"


#define LOG2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))

void test_bitpolymul_lch(){
	unsigned m = 3; // degree n = 2^m
	unsigned n = (1ULL) << m; // n is even and it must be even for keeping the 32 byte
    printf("m = %u, n = %u\n", m, n);
    // Generating the random polynomial in GF_2^128
	__m128i* fx = random_polynomial_gf2128(n);
    __m128i* evals = lch_fft_gf2128(fx, n);
    __m128i* evals2 = naive_evaluate(fx, n);
    printf("are equal = %b\n", are_equal_vec_128bit(evals, evals2, n));
}

void test_cantor(){
    unsigned m = 15; // degree n = 2^m
	unsigned n = (1ULL) << m; // n is even and it must be even for keeping the 32 byte
    printf("m = %u, n = %u\n", m, n);
    // Generating the random polynomial in GF_2^128
	__m128i* fx = random_polynomial_gf2128(n);
    __m128i* evals = cantor_fft_gf2128(fx, n);
    __m128i* evals3 = cantor_fft_hc_circuits_gf2128(fx, n);
    // __m128i* evals2 = naive_evaluate(fx, n);
    printf("are equal = %b\n", are_equal_vec_128bit(evals, evals3, n));
}

void test_cantor_parallel(){
    unsigned m = 20; // degree n = 2^m
	unsigned n = (1ULL) << m; // n is even and it must be even for keeping the 32 byte
    printf("m = %u, n = %u\n", m, n);
    // Generating the random polynomial in GF_2^128
	__m128i* fx = random_polynomial_gf2128(n);
    // __m128i* evals = cantor_fft_gf2128(fx, n);
    #pragma omp parallel
        {
            #pragma omp single
            {
                __m128i* evals3 = cantor_fft_gf2128_parallel(fx, n);
            }
        }
    // __m128i* evals2 = naive_evaluate(fx, n);
    // printf("are equal = %b\n", are_equal_vec_128bit(evals, evals3, n));
}

void test_dyadic(){
    for (unsigned m = 1; m <= 15; ++m){
        unsigned n = (1ULL) << m;
        __m128i* fx = random_polynomial_gf2128(n);
        __m128i* evals = dyadic_fft_gf2128(fx, n);
        __m128i* evals2 = cantor_fft_gf2128(fx, n);
        bool ok = are_equal_vec_128bit(evals, evals2, n);
        if (m <= 10){ // naive evaluation is quadratic; only check small sizes
            __m128i* evals3 = naive_evaluate(fx, n);
            ok &= are_equal_vec_128bit(evals, evals3, n);
            free(evals3);
        }
        printf("m = %u, n = %u, dyadic equals naive = %b\n", m, n, ok);
        free(fx);
        free(evals);
        free(evals2);
    }
}

#define ITERATIONS 100
void cantor_vs_lch(){
    printf("m\tCantor\t\tCantor HC\tCantor PARALLEL\tLCH\n");
    for (unsigned m = 9; m < 22; m++){
	    unsigned n = (1ULL) << m; 
        clock_t start, end;
        double time_cantor = 0, time_cantor_hc = 0, time_lch = 0, time_cantor_hc_cache=0;
        
        for (unsigned iter = 0; iter < ITERATIONS; ++iter){
            __m128i* fx1 = random_polynomial_gf2128(n);
            start = clock();
            __m128i* evals = cantor_fft_gf2128(fx1, n);
            end = clock();
            time_cantor += 1000 * ((double) (end - start)) / CLOCKS_PER_SEC;
            free(fx1);
            free(evals);
        }

        #pragma omp parallel
        {
            #pragma omp single
            {
        for (unsigned iter = 0; iter < ITERATIONS; ++iter){
            __m128i* fx1 = random_polynomial_gf2128(n);
            start = clock();
            __m128i* evals = cantor_fft_gf2128_parallel(fx1, n);
            end = clock();
            time_cantor_hc_cache += 1000 * ((double) (end - start)) / CLOCKS_PER_SEC;
            free(fx1);
            free(evals);
        }
    }       
}

        for (unsigned iter = 0; iter < ITERATIONS; ++iter){
            __m128i* fx1 = random_polynomial_gf2128(n);
            start = clock();
            __m128i* evals = cantor_fft_hc_circuits_gf2128(fx1, n);
            end = clock();
            time_cantor_hc += 1000 * ((double) (end - start)) / CLOCKS_PER_SEC;
            free(fx1);
            free(evals);
        }

        for (unsigned iter = 0; iter < ITERATIONS; ++iter){
            __m128i* fx1 = random_polynomial_gf2128(n);
            start = clock();
            __m128i* evals = lch_fft_gf2128(fx1, n);
            end = clock();
            time_lch += 1000 * ((double) (end - start)) / CLOCKS_PER_SEC;
            free(fx1);
            free(evals);
        }
    

        printf("%u\t%f ms\t%f ms\t%f ms\t%f ms\n", m, time_cantor/ITERATIONS, time_cantor_hc/ITERATIONS, time_cantor_hc_cache/ITERATIONS, time_lch/ITERATIONS);
    }
}

// Monotonic wall-clock time in milliseconds (C11 timespec_get).
static double now_ms(void){
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1e6;
}

static double mean_of(const double* x, unsigned n){
    double s = 0;
    for (unsigned i = 0; i < n; ++i) s += x[i];
    return s / n;
}

// sample standard deviation (n - 1 denominator)
static double stddev_of(const double* x, unsigned n, double mean){
    double s = 0;
    for (unsigned i = 0; i < n; ++i) s += (x[i] - mean) * (x[i] - mean);
    return sqrt(s / (n - 1));
}

void lch_vs_dyadic(){
    static double t_lch[ITERATIONS], t_dyadic[ITERATIONS];

    // Use library B (tuned via `make tune`, else L2 heuristic). Always log
    // B and its source so paper runs are reproducible from the log + tune file.
    dyadic_fft_set_cache_block(0);
    printf("Dyadic cache block B = %u elements (%u KB)  [source=%s]\n",
           dyadic_fft_get_cache_block(), dyadic_fft_get_cache_block() / 64,
           dyadic_fft_cache_block_source());
    printf("iterations = %u (interleaved LCH/Dyadic, mean +- sample stddev)\n",
           (unsigned) ITERATIONS);
    printf("m\tLCH (ms)\t\t\tDyadic (ms)\n");
    for (unsigned m = 9; m < 29; m++){
	    unsigned n = (1ULL) << m; 

        // Interleave the two algorithms inside each iteration so both see
        // the same machine conditions (frequency scaling, thermal drift).
        for (unsigned iter = 0; iter < ITERATIONS; ++iter){
            __m128i* fx1 = random_polynomial_gf2128(n);
            double t0 = now_ms();
            __m128i* evals = lch_fft_gf2128(fx1, n);
            t_lch[iter] = now_ms() - t0;
            free(fx1);
            free(evals);

            __m128i* fx2 = random_polynomial_gf2128(n);
            t0 = now_ms();
            __m128i* evals2 = dyadic_fft_gf2128(fx2, n);
            t_dyadic[iter] = now_ms() - t0;
            free(fx2);
            free(evals2);
        }

        double mean_lch = mean_of(t_lch, ITERATIONS);
        double mean_dyadic = mean_of(t_dyadic, ITERATIONS);
        printf("%u\t%.4f +- %.4f\t\t%.4f +- %.4f\n", m,
               mean_lch, stddev_of(t_lch, ITERATIONS, mean_lch),
               mean_dyadic, stddev_of(t_dyadic, ITERATIONS, mean_dyadic));
    }
}

int main(){
    srand(time(NULL));
    // validate_cantor_basis();
    // test_bitpolymul_lch();
    test_dyadic();
    // test_cantor();
    // cantor_vs_lch();
    lch_vs_dyadic();
    // test_cantor_parallel();
    return 0;
}
