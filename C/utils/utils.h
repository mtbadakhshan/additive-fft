#ifndef _UTILS_H_
#define _UTILS_H_
#include <immintrin.h>
#include "bitpolymul/byte_inline_func.h"

static inline void print_vector_128bit_name(__m128i * vec, unsigned n, const char * vec_name){printf("%s :", vec_name); xmm_dump(vec, n); puts("");}
static inline void print_vector_128bit(__m128i * vec, unsigned n){print_vector_128bit_name(vec, n, "vec");}


#endif
