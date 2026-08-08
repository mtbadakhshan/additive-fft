#include "dyadic_fft.h"
#include <stdio.h>
#include <string.h>    // Header for memcpy

#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>    // sysconf: L2 cache size
#endif

#include "bitpolymul/bitmat_prod.h"
#include "bitpolymul/gf2128_cantor_iso.h"
#include "bitpolymul/gfext_aesni.h"


#define LOG2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))

// Cache chunk B (16-byte elements) for blocked Taylor + subtree schedule.
//
// B is a *machine parameter*, not part of the algorithm — same idea as gf2x
// or FFTW "wisdom": the blocked schedule is fixed; B is chosen by an
// explicit tuning phase (`make tune` -> dyadic_tune_params.h) for
// reproducibility. Untuned fallback: round_pow2(L2/16), then 2^16.
#include "dyadic_tune_params.h"

#ifndef DYADIC_FALLBACK_BLOCK
#define DYADIC_FALLBACK_BLOCK (1u << 16)
#endif

#define DYADIC_MIN_BLOCK (1u << 10)
#define DYADIC_MAX_BLOCK (1u << 22)

static unsigned dyadic_l2_block = 0;       // 0 = not resolved yet (untuned path)
static unsigned dyadic_forced_block = 0;   // non-zero = runtime override
static unsigned dyadic_block = DYADIC_FALLBACK_BLOCK;
static const char* dyadic_block_source = "fallback";

static unsigned dyadic_clamp_block(unsigned long n_elem){
    if (n_elem < DYADIC_MIN_BLOCK) return DYADIC_MIN_BLOCK;
    if (n_elem > DYADIC_MAX_BLOCK) return DYADIC_MAX_BLOCK;
    return 1u << LOG2(n_elem);      // round down to a power of two
}

static unsigned dyadic_configured_block(void){
    if (dyadic_forced_block){
        dyadic_block_source = "forced";
        return dyadic_forced_block;
    }
#ifdef DYADIC_TUNED_BLOCK
    dyadic_block_source = "tuned";
    return dyadic_clamp_block(DYADIC_TUNED_BLOCK);
#else
    if (!dyadic_l2_block){
#ifdef _SC_LEVEL2_CACHE_SIZE
        long l2_bytes = sysconf(_SC_LEVEL2_CACHE_SIZE);
        if (l2_bytes > 0){
            dyadic_l2_block = dyadic_clamp_block((unsigned long) l2_bytes / sizeof(__m128i));
            dyadic_block_source = "L2";
            return dyadic_l2_block;
        }
#endif
        dyadic_l2_block = DYADIC_FALLBACK_BLOCK;
        dyadic_block_source = "fallback";
    }
    return dyadic_l2_block;
#endif
}

void dyadic_fft_set_cache_block(unsigned n_elem){
    if (n_elem == 0){
        dyadic_forced_block = 0;
        dyadic_l2_block = 0;
        return;
    }
    dyadic_forced_block = dyadic_clamp_block(n_elem);
}

unsigned dyadic_fft_get_cache_block(void){
    return dyadic_configured_block();
}

unsigned dyadic_fft_cache_block_for(unsigned n_term){
    (void) n_term;
    return dyadic_configured_block();
}

const char* dyadic_fft_cache_block_source(void){
    (void) dyadic_configured_block();
    return dyadic_block_source;
}

// One halving pass: poly[p - shift] ^= poly[p] for p descending over the
// half-length source range starting at top. The descending order is
// semantic: sources in [top, top + half - shift) must be read after they
// were updated from distance shift above. A 2-wide (256-bit) descending
// sweep preserves that order whenever shift >= 2, because any update to a
// chunk's sources comes from a strictly earlier (higher) chunk; the only
// shift == 1 case (half == 2, mu1 == 1, s == 0) stays scalar.
// Source loads are 32B-aligned; targets may be misaligned by shift.
static inline void dyadic_xor_down(__m128i* top, unsigned half, unsigned shift){
    if (shift >= 2){
        __m128i* dst = top - shift;
        for (unsigned i = half; i >= 2; i -= 2){
            __m256i v = _mm256_load_si256((const __m256i*)(top + i - 2));
            __m256i d = _mm256_loadu_si256((const __m256i*)(dst + i - 2));
            _mm256_storeu_si256((__m256i*)(dst + i - 2), d ^ v);
        }
    } else {
        for (unsigned i = half; i-- > 0; ){
            top[(int)i - (int)shift] ^= top[i];
        }
    }
}

// In-place Taylor expansion phase T(mu, s): the array is partitioned into
// independent subproblems of length 2^mu stored with stride 2^s. Each
// subproblem is Taylor-expanded w.r.t. S^{mu1}(x) = x^(2^mu1) + x, where mu1
// is the largest power of two strictly less than mu (valid for a Cantor
// basis). All 2^s interleaved subproblems inside a block of size 2^(s+mu)
// share the same additions, so each block is processed with plain XOR sweeps
// over physical positions.
//
// Halving schedule (Gao-Mateer radix conversion, cf. additive_FFT_CO in
// C++/Gao/fft.tcc): for a piece of logical length L, split
// f = f0 + x^(L/2) f1 and use (x^t + x)^(L/(2t)) = x^(L/2) + x^(L/(2t))
// over GF(2), with t = 2^mu1. One descending pass
// f[p - (L/2 - L/(2t))] ^= f[p] over the top half rewrites f as
// g0 + S^{mu1}(x)^(L/(2t)) g1 (the descending order matters: additions that
// land back inside the top half are folded down again by later iterations);
// then recurse on both halves. This costs mu2 = mu - mu1 passes instead of
// the naive 2^mu2 shrinking re-sweeps.
//
// Cache blocking: running each halving level as one sweep over the whole
// array re-streams the working set from DRAM once per level. Levels whose
// span exceeds the cache chunk are unavoidable full sweeps; all remaining
// levels of a cache-sized block are executed back-to-back while the block
// is hot (cf. the depth-first recursion in bc_to_lch_128), so they cost one
// memory pass in total instead of one per level.
static void dyadic_taylor_phase(__m128i* poly, unsigned n_term, unsigned mu, unsigned s){
    unsigned mu1 = 1;
    while (mu1 * 2 < mu) mu1 *= 2;

    unsigned bottom = (1u << s) << mu1;  // stop when pieces reach logical length t = 2^mu1
    unsigned len = (1u << s) << mu;

    // levels spanning more than a cache block: full sweeps
    for (; len > bottom && len > dyadic_block; len >>= 1){
        unsigned half = len >> 1;
        unsigned shift = half - (half >> mu1);  // physical distance target <- source
        for (unsigned b0 = 0; b0 < n_term; b0 += len){
            dyadic_xor_down(poly + b0 + half, half, shift);
        }
    }
    if (len <= bottom) return;

    // remaining levels: finish each cache-sized block completely before
    // moving to the next one
    for (unsigned b0 = 0; b0 < n_term; b0 += len){
        for (unsigned l = len; l > bottom; l >>= 1){
            unsigned half = l >> 1;
            unsigned shift = half - (half >> mu1);
            for (unsigned c0 = b0; c0 < b0 + len; c0 += l){
                dyadic_xor_down(poly + c0 + half, half, shift);
            }
        }
    }
}

// In-place butterfly stage B(2^s) restricted to [base, base + region): for
// every pair (p, p + 2^s) with bit_s(p) = 0, evaluate the length-2 polynomial
// at {lam, lam + 1}. With zero affine shift, lam has Cantor-basis coordinate
// p >> s = module << 1 (module index is global, hence derived from base),
// mapped to GF(2^128) via the gfCantorto2128_8R bit-matrix (as in
// cantor_fft.c).
static void dyadic_butterfly_stage(__m128i* poly, unsigned base, unsigned region, unsigned s){
    unsigned half = 1u << s;
    unsigned module0 = base >> (s + 1);
    unsigned n_modules = region >> (s + 1);
    unsigned l = 0;

    // module = 0: lam = 0, only the XOR half of the butterfly remains
    if (module0 == 0){
        if (half >= 2){
            __m256i* poly256 = (__m256i*)(poly + base);
            unsigned half_half = half >> 1; // 256-bit steps
            for (unsigned i = 0; i < half_half; ++i){
                poly256[half_half + i] ^= poly256[i];
            }
        } else {
            poly[base + 1] ^= poly[base];
        }
        l = 1;
    }

    for (; l < n_modules; ++l){
        __m128i mult_factor = bitmat_prod_accu_64x128_M8R_sse(_mm_setzero_si128(),  gfCantorto2128_8R, (module0 + l)<<1);
        __m128i* p0 = poly + base + (l << (s + 1));
        __m128i* p1 = p0 + half;
        for (unsigned i = 0; i < half; ++i){
            p0[i] ^= _gf2ext128_mul_sse(p1[i], mult_factor);
            p1[i] ^= p0[i];
        }
    }
}

// Fused pair of butterfly stages B(2^(s+1)) followed by B(2^s), i.e. the whole
// remaining work of a (mu = 2, s) task after its Taylor phase. For each i,
// the 4-tuple (x0[i], x1[i], x2[i], x3[i]) taken at stride sigma = 2^s is one
// 2x2 block (two length-2 column FFTs, then two length-2 row FFTs), evaluated
// in a single memory pass. For group h the stage-1 twiddle is lam = C(2h) and the
// stage-2 twiddles are mu0 = C(4h) and mu1 = C(4h + 2) = mu0 + C(2), where C
// is the Cantor-to-GF(2^128) bit-matrix map; so only 2 bit-matrix products
// (plus one XOR) are needed per group instead of 3, and each 4-tuple is
// loaded and stored once instead of twice.
static void dyadic_butterfly_2x2(__m128i* poly, unsigned base, unsigned region, unsigned s){
    unsigned sigma = 1u << s;              // quarter length
    unsigned group0 = base >> (s + 2);
    unsigned n_groups = region >> (s + 2);
    const __m128i c2 = _mm_load_si128((const __m128i*) gfCantorto2128_8R + 2); // C(2)
    unsigned g = 0;

    // group 0: lam = 0, mu0 = 0, mu1 = C(2); only one multiplication survives
    if (group0 == 0){
        __m128i* x0 = poly + base;
        __m128i* x1 = x0 + sigma;
        __m128i* x2 = x1 + sigma;
        __m128i* x3 = x2 + sigma;
        for (unsigned i = 0; i < sigma; ++i){
            x2[i] ^= x0[i];
            x3[i] ^= x1[i];
            x1[i] ^= x0[i];
            x2[i] ^= _gf2ext128_mul_sse(x3[i], c2);
            x3[i] ^= x2[i];
        }
        g = 1;
    }

    for (; g < n_groups; ++g){
        unsigned h = group0 + g;
        __m128i lam = bitmat_prod_accu_64x128_M8R_sse(_mm_setzero_si128(), gfCantorto2128_8R, h << 1);
        __m128i mu0 = bitmat_prod_accu_64x128_M8R_sse(_mm_setzero_si128(), gfCantorto2128_8R, h << 2);
        __m128i mu1 = mu0 ^ c2;
        __m128i* x0 = poly + base + (g << (s + 2));
        __m128i* x1 = x0 + sigma;
        __m128i* x2 = x1 + sigma;
        __m128i* x3 = x2 + sigma;
        for (unsigned i = 0; i < sigma; ++i){
            __m128i a0 = x0[i], a1 = x1[i], a2 = x2[i], a3 = x3[i];
            a0 ^= _gf2ext128_mul_sse(a2, lam);
            a2 ^= a0;
            a1 ^= _gf2ext128_mul_sse(a3, lam);
            a3 ^= a1;
            a0 ^= _gf2ext128_mul_sse(a1, mu0);
            a1 ^= a0;
            a2 ^= _gf2ext128_mul_sse(a3, mu1);
            a3 ^= a2;
            x0[i] = a0; x1[i] = a1; x2[i] = a2; x3[i] = a3;
        }
    }
}

// Recursive schedule of the dyadic AFFT restricted to [base, base + region):
//   Sched(mu, s) = T(mu, s) + Sched(mu2, s + mu1) + Sched(mu1, s),
//   Sched(1, s) = B(2^s),  Sched(2, s) = T(2, s) + fused 2x2 kernel.
// Each butterfly stride 2^s, s = m-1 .. 0, occurs exactly once overall.
//
// Cache blocking across the whole subtree: every operation of Sched(mu, s)
// stays inside span-aligned blocks of size span = 2^(s+mu), so blocks are
// fully independent. Once span fits in the configured chunk B, the entire
// remaining subtree is executed chunk by chunk while the chunk is hot -
// one trip through DRAM instead of one per phase. Only phases whose span
// exceeds B remain full sweeps.
static void dyadic_sched(__m128i* poly, unsigned base, unsigned region, unsigned mu, unsigned s){
    unsigned span = (1u << s) << mu;
    if (region > dyadic_block && span <= dyadic_block){
        for (unsigned b = base; b < base + region; b += dyadic_block){
            dyadic_sched(poly, b, dyadic_block, mu, s);
        }
        return;
    }
    if (mu == 1){
        dyadic_butterfly_stage(poly, base, region, s);
        return;
    }
    if (mu == 2){
        dyadic_taylor_phase(poly + base, region, 2, s);
        dyadic_butterfly_2x2(poly, base, region, s);
        return;
    }
    unsigned mu1 = 1; // largest power of 2 strictly less than mu
    while (mu1 * 2 < mu) mu1 *= 2;
    dyadic_taylor_phase(poly + base, region, mu, s);
    dyadic_sched(poly, base, region, mu - mu1, s + mu1);  // column FFTs
    dyadic_sched(poly, base, region, mu1, s);             // row FFTs
}

__m128i* dyadic_fft_gf2128(__m128i* fx, unsigned n_term){
    #ifdef NCOPY_POLY // operate on fx in place; otherwise work on a copy so the input stays unchanged
    __m128i* poly = fx;
    #else
    __m128i* poly = (__m128i*)aligned_alloc( 64 , sizeof(__m128i)*n_term );
    poly = memcpy( poly, fx, sizeof(__m128i)*n_term );
    #endif
    unsigned m = LOG2(n_term);

    dyadic_block = dyadic_configured_block();
    dyadic_sched(poly, 0, n_term, m, 0);

    return poly;
}

// Reference variant of the driver: the same schedule flattened into an
// explicit stack of (mu, s) tasks, executing every phase as a full pass over
// the array (no subtree cache blocking). Kept as a readable reference and
// correctness baseline: it performs the same arithmetic as dyadic_fft_gf2128
// (the blocked recursion only reorders operations acting on disjoint
// blocks), so the outputs are bit-identical.
__m128i* dyadic_fft_gf2128_iter(__m128i* fx, unsigned n_term){
    #ifdef NCOPY_POLY // operate on fx in place; otherwise work on a copy so the input stays unchanged
    __m128i* poly = fx;
    #else
    __m128i* poly = (__m128i*)aligned_alloc( 64 , sizeof(__m128i)*n_term );
    poly = memcpy( poly, fx, sizeof(__m128i)*n_term );
    #endif
    unsigned m = LOG2(n_term);

    dyadic_block = dyadic_configured_block();
    unsigned stack_mu[32], stack_s[32];
    unsigned top = 0;
    stack_mu[top] = m; stack_s[top] = 0; top++;
    while (top){
        top--;
        unsigned mu = stack_mu[top], s = stack_s[top];
        // Follow the column-FFT chain; row tasks are stacked for later.
        // Stop at mu = 2: those tasks are handled by the fused 2x2 kernel.
        while (mu > 2){
            unsigned mu1 = 1; // largest power of 2 strictly less than mu
            while (mu1 * 2 < mu) mu1 *= 2;
            unsigned mu2 = mu - mu1;
            dyadic_taylor_phase(poly, n_term, mu, s);
            stack_mu[top] = mu1; stack_s[top] = s; top++;
            mu = mu2; s += mu1;
        }
        if (mu == 2){
            dyadic_taylor_phase(poly, n_term, 2, s);
            dyadic_butterfly_2x2(poly, 0, n_term, s);
        } else {
            dyadic_butterfly_stage(poly, 0, n_term, s);
        }
    }

    return poly;
}
