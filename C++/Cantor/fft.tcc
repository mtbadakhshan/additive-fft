#include <Cantor/fft.hpp>
#include <array>
#include <cstddef>
#include <libff/algebra/field_utils/algorithms.hpp>

#include <iostream>
#include "utils/utils.hpp"
#include "Cantor/cantor_basis.hpp"

namespace cantor {

// The precomputations are hard coded. This is for the case that the affine shift is zero
template<typename FieldT>   
std::vector<FieldT> additive_FFT(const std::vector<FieldT> &poly_coeffs, const size_t domain_dim, const size_t shift_dim){
    const size_t m = domain_dim;
    std::vector<FieldT> g(poly_coeffs);
    g.resize((1ULL << m), FieldT::zero());
    const size_t n = g.size();
    assert(n == (1ull<<m));

    FieldT* cantor_combinations;
    if(FieldT::extension_degree() == 128) cantor_combinations = (FieldT*) cantor_combinations_8R_in_gf2to128;
    else if(FieldT::extension_degree() == 192) cantor_combinations = (FieldT*) cantor_combinations_8R_in_gf2to192;
    else if(FieldT::extension_degree() == 256) cantor_combinations = (FieldT*) cantor_combinations_8R_in_gf2to256;
    else throw std::invalid_argument("The field size should be either 128, or 256 for using the cantor basis");

    size_t input_size = n;
    size_t n_modules = 1;

    const unsigned* nz_S;
    size_t S_index = m, t;
    for (int r = 0; r < m; ++r)
    {   
        --S_index; 
        nz_S = s_i[S_index]; t=n_terms[S_index];
        size_t offset = 0;
        size_t half_input_size = input_size >> 1;
        size_t shift_bit = shift_dim==0 ? 0 :n_modules << (shift_dim - m);


        for (size_t module = 0; module < n_modules; ++module){
            size_t module_shifted = (module|shift_bit) << 1, i_256 = 0;
            FieldT mult_factor = FieldT::zero();
            while(module_shifted){
                mult_factor += cantor_combinations[i_256 + (module_shifted & 0xff)];
                i_256 += 256;
                module_shifted >>= 8;
            }
            size_t offset2 = offset + half_input_size;
            for (unsigned k = offset2 + half_input_size - 1; k >= offset2; --k){
                FieldT gk = g[k];
                for (unsigned i = 0; i < t; ++i){
                    g[k - nz_S[i]] += gk;
                }
                g[k-half_input_size] += gk * mult_factor;
            }
            for (size_t j = 0; j < half_input_size; ++j)
                g[offset2+j] += g[offset+j];

            offset += input_size;
        }
        input_size = half_input_size;
        n_modules <<= 1;
    }
    return g;
}

template<typename FieldT>   
std::vector<FieldT> additive_IFFT(const std::vector<FieldT> &evals, 
                                 const size_t domain_dim, const size_t shift_dim){
    const size_t m = domain_dim;
    std::vector<FieldT> g(evals);
    // g.resize((1ULL << m), FieldT::zero());
    const size_t n = g.size();
    assert(n == (1ull<<m));

    FieldT* cantor_combinations;
    if(FieldT::extension_degree() == 128) cantor_combinations = (FieldT*) cantor_combinations_8R_in_gf2to128;
    else if(FieldT::extension_degree() == 192) cantor_combinations = (FieldT*) cantor_combinations_8R_in_gf2to192;
    else if(FieldT::extension_degree() == 256) cantor_combinations = (FieldT*) cantor_combinations_8R_in_gf2to256;
    else throw std::invalid_argument("The field size should be either 128, or 256 for using the cantor basis");

    size_t input_size = 1;
    size_t n_modules = n;

    const unsigned* nz_S;
    size_t S_index = 0, t;
    for (int r = m-1; r >=0; --r)
    {   
        nz_S = s_i[S_index]; t=n_terms[S_index];
        n_modules >>= 1;
        size_t half_input_size = input_size;
        input_size <<= 1;
        size_t shift_bit = shift_dim==0 ? 0 :n_modules << (shift_dim - m);
        size_t offset = 0;
        for (size_t module = 0; module < n_modules; ++module){
            size_t module_shifted = (module|shift_bit) << 1, i_256 = 0;
            // std::cout<<"module_shifted: "<< std::bitset<16>(module_shifted) << std::endl;
            FieldT mult_factor = FieldT::zero();
            while(module_shifted){
                mult_factor += cantor_combinations[i_256 + (module_shifted & 0xff)];
                i_256 += 256;
                module_shifted >>= 8;
                // std::cout<<"module_shifted: "<< std::bitset<16>(module_shifted) << std::endl;
            }
            size_t offset2 = offset + half_input_size;
            for (size_t j = 0; j < half_input_size; ++j)
                g[offset2+j] += g[offset+j];

            for (size_t k =  offset2; k < offset2+half_input_size; ++k){
                FieldT gk = g[k];
                g[k-half_input_size] += gk * mult_factor;
                for (size_t i = 0; i < t; ++i){
                    g[k - nz_S[i]] += gk;
                }
            }
            offset += input_size;
        }
        S_index ++;
    }
    return g;
}

template<typename FieldT>
PreComputedValues<FieldT> pre_computation(const libiop::affine_subspace<FieldT> &domain)
{
    const size_t m = domain.dimension();
    FieldT affine_shift = domain.shift();
    std::vector<FieldT> W(domain.basis());

    PreComputedValues<FieldT> values(m);

    size_t cnt = 0;
    size_t S_index = m;
    size_t n_modules = 1;
    for (size_t r = 0; r < m; ++r)
    {
        FieldT affine_shift_round = affine_shift;
        --S_index; 
        std::vector<size_t> nz_S;
        if(S_index > 0){
            nz_S.reserve(S_index);
            nz_S.emplace_back((1<<S_index)-1); // C(x,0) = 1, so we assigned that before loop.
            for (size_t i = 1; i < S_index; ++i){
                size_t ii = i ;
                size_t rr = S_index;
                while (((rr & 1) | (~ii & 1)) && ii>0){
                    ii >>= 1;
                    rr >>= 1;
                }
                if (ii == 0){
                    nz_S.emplace_back((1<<S_index) - (1<<i));
                    affine_shift_round += affine_shift ^ (1<<i);
                }
            }
            affine_shift_round += affine_shift ^ (1<<S_index);
        }  

        for (size_t module = 0; module < n_modules; ++module){
            
                // Computing the multiplication factor
                FieldT mult_factor = affine_shift_round;
                for (size_t i = 0; i < r; ++i){
                    if (module & (1<<i))
                        mult_factor += W[i+1];
                }
                values.mult_factor_vec.emplace_back(mult_factor);
        }

        n_modules <<= 1;
        values.nz_S_vecs.emplace_back(nz_S);
    }
    return values;
}

// The pricomputed values are passed to the algorithm.
template<typename FieldT>
std::vector<FieldT> additive_FFT(const std::vector<FieldT> &poly_coeffs, 
                                        const PreComputedValues<FieldT> &values)
{
    const size_t m = values.dimension;

    std::vector<FieldT> g(poly_coeffs);
    g.resize(1ull<<m, FieldT::zero());
    const size_t n = g.size();
    assert(n == (1ull<<m));


    size_t input_size = n;
    size_t n_modules = 1;
    size_t cnt = 0;
    for (size_t r = 0; r < m; ++r)
    {
        std::vector<size_t> nz_S(values.nz_S_vecs[r]);
        size_t offset = 0;
        size_t half_input_size = input_size >> 1;

        for (size_t module = 0; module < n_modules; ++module){
            FieldT mult_factor = values.mult_factor_vec[cnt++];
            size_t offset2 = offset + half_input_size;
            for (size_t k = offset2+half_input_size -1; k >= offset2; --k){
                FieldT gk = g[k];
                for (const auto& nz : nz_S)
                    g[k - nz] += gk ;
                g[k-half_input_size] += gk * mult_factor;
            }

            for (size_t j = 0; j < half_input_size; ++j)
                g[offset2+j] += g[offset+j];

            offset += input_size;
        }
        input_size = half_input_size;
        n_modules <<= 1;
    }
    return g;
}

template<typename FieldT>
std::vector<FieldT> additive_IFFT(const std::vector<FieldT> &evals, 
                                        const PreComputedValues<FieldT> &values)
{
    const size_t m = values.dimension;

    std::vector<FieldT> g(evals);
    const size_t n = g.size();
    assert(n == (1ull<<m));


    size_t input_size = 1;
    size_t n_modules = n;
    size_t cnt = n-2;
    for (int r = m-1; r >= 0; --r)
    {
        std::vector<size_t> nz_S(values.nz_S_vecs[r]);
        n_modules >>= 1;
        size_t half_input_size = input_size;
        input_size <<= 1;

        size_t offset = n - input_size;
        for (int module = n_modules - 1; module >= 0; --module){
            FieldT mult_factor = values.mult_factor_vec[cnt--];
            size_t offset2 = offset + half_input_size;
            for (size_t j = 0; j < half_input_size; ++j)
                g[offset2+j] += g[offset+j];
            for (size_t k = offset2; k < offset2+half_input_size; ++k){
                FieldT gk = g[k];
                g[k-half_input_size] += gk * mult_factor;
                for (const auto& nz : nz_S)
                    g[k - nz] += gk ;
            }
            offset -= input_size;
        }
    }
    return g;
}


template<typename FieldT>
std::vector<FieldT> additive_FFT(const std::vector<FieldT> &poly_coeffs,
                                        const libiop::affine_subspace<FieldT> &domain)
{
    std::vector<FieldT> g(poly_coeffs);
    g.resize(domain.num_elements(), FieldT::zero());
    
    FieldT affine_shift = domain.shift();
    std::vector<FieldT> W(domain.basis());

    const size_t n = g.size();
    const size_t m = domain.dimension();
    assert(n == (1ull<<m));


    size_t input_size = n;
    size_t n_modules = 1;
    size_t S_index = m;
    for (size_t r = 0; r < m; ++r)
    {
        FieldT affine_shift_round = affine_shift;
        // Computing non-zero indices in S_{m-r}(X) for round r in the reverse order. Becuse in division algorithm we start from 
        --S_index; 
        std::vector<size_t> nz_S;

        if(S_index > 0){
            nz_S.reserve(S_index);
            nz_S.emplace_back((1<<S_index)-1); // C(x,0) = 1, so we assigned that before loop.
            for (size_t i = 1; i < S_index; ++i){
                size_t ii = i ;
                size_t rr = S_index;
                while (((rr & 1) | (~ii & 1)) && ii>0){
                    ii >>= 1;
                    rr >>= 1;
                }
                if (ii == 0){
                    nz_S.emplace_back((1<<S_index) - (1<<i));
                    affine_shift_round += affine_shift ^ (1<<i);
                }
            }
            affine_shift_round += affine_shift ^ (1<<S_index);
        }  

        size_t offset = 0;
        size_t half_input_size = input_size >> 1;

        for (size_t module = 0; module < n_modules; ++module){
            
            // Computing the multiplication factor
            FieldT mult_factor = affine_shift_round;
            for (size_t i = 0; i < r; ++i){
                if (module & (1<<i))
                    mult_factor += W[i+1];
            }
            size_t offset2 = offset + half_input_size;
            for (size_t k = offset2+half_input_size -1; k >= offset2; --k){
                FieldT gk = g[k];
                for (const auto& nz : nz_S)
                    g[k - nz] += gk ;
                g[k-half_input_size] += gk * mult_factor;
            }
            for (size_t j = 0; j < half_input_size; ++j)
                g[offset2+j] += g[offset+j];

            offset += input_size;
        }
        input_size = half_input_size;
        n_modules <<= 1;
    }

    return g;
}

// Radix-4 variant of additive_FFT(poly_coeffs, domain), thin wrapper around
// additive_FFT_radix2k<2>.
//
// Earlier, this function was implemented as a hand-coded "two-sweep" radix-4
// round (one descending sweep that fused level-1 division + T_2 twiddle +
// level-2 division on both halves, plus one ascending butterfly closer). That
// design is mathematically incorrect for m >= 3: the level-2 writes from
// iters in [3L/4, L) land in upper-half positions in [L/2, L-2], which are
// then read by the level-1 division at later iters in [L/2, 3L/4),
// corrupting the round-0 quotient bits.
//
// The correct memory-efficient radix-4 round needs 3 sweeps (one per level
// + one closer), achieved by additive_FFT_radix2k<2>. See RADIX4_CANTOR_FFT.md.
template<typename FieldT>
std::vector<FieldT> additive_FFT_radix4(const std::vector<FieldT> &poly_coeffs,
                                        const libiop::affine_subspace<FieldT> &domain)
{
    return additive_FFT_radix2k<2, FieldT>(poly_coeffs, domain);
}

// Templated radix-2^K Cantor AFFT, generalizing additive_FFT_radix4.
//
// Per round j (0-indexed, doing floor(m/K) rounds total):
//   * module size L = 2^{m - K*j}, n_modules = 2^{K*j}
//   * sub-chunk size = L / 2^K
//
// Each radix-2^K round runs in K+1 array sweeps per module:
//
//   (Sweep ell, ell=1..K) Descending k-loop firing the level-ell division
//             g[k - nz_S[ell]] += g[k] for all nz in nz_S[ell], at all
//             positions k whose bit (logL - ell) is set.
//             Sweep 1 also folds in the stage-1 cross-half twiddle
//                  g[k - L/2] += g[k] * T_{u, 1, 0}    for k in [L/2, L).
//             We do NOT fuse level-ell+1 writes into the level-ell descending
//             pass: level-(ell+1) writes from iter k_iter to upper-half
//             positions can pollute g[k] read by level-ell at a later iter.
//             The K levels MUST be done in separate descending passes.
//
//   (Sweep K+1) Ascending tuple-loop over jj in [0, L/2^K). For each tuple,
//             load 2^K values into registers, run stages 2..K of the Cantor
//             butterfly in-register, write 2^K values back. Stage-1's "+="
//             closer is also performed here.
//
// Per-tuple twiddles for stages 2..K are derived from K base twiddles
// base_T[1..K], where base_T[ell] = S^{m-K*j-ell}(theta_module). Within a
// stage, block b's twiddle is base_T[ell] + sum_{i: b_i=1} W[i+1].
//
// Total memory passes per FFT (treating each sweep as ~one pass over the
// active footprint):
//   * radix-2:           2m  passes
//   * radix-2^K (K>=1):  (K+1) * m / K  passes  (asymptotes to m as K->inf)
//
// When m is not a multiple of K, the function performs (m/K) radix-2^K rounds
// and then does (m mod K) residual radix-2 rounds (2 sweeps each).
template<size_t K, typename FieldT>
std::vector<FieldT> additive_FFT_radix2k(const std::vector<FieldT> &poly_coeffs,
                                          const libiop::affine_subspace<FieldT> &domain)
{
    static_assert(K >= 1, "K must be >= 1");
    static_assert(K <= 5, "K supported up to 5 (radix up to 32)");

    std::vector<FieldT> g(poly_coeffs);
    g.resize(domain.num_elements(), FieldT::zero());

    FieldT affine_shift = domain.shift();
    std::vector<FieldT> W(domain.basis());

    const size_t n = g.size();
    const size_t m = domain.dimension();
    assert(n == (1ull << m));

    constexpr size_t Q      = 1ull << K;     // = 2^K, the radix
    constexpr size_t Q_half = Q >> 1;        // = 2^{K-1}

    // block_offsets[c] = sum of W[i+1] for set bits i of c (c < 2^{K-1}).
    // Used in Sweep K+1 to derive in-stage block twiddles from base twiddles.
    std::array<FieldT, (Q_half == 0 ? 1 : Q_half)> block_offsets;
    block_offsets[0] = FieldT::zero();
    for (size_t c = 1; c < Q_half; ++c) {
        const size_t low_bit = __builtin_ctzll(c);
        block_offsets[c] = block_offsets[c ^ (1ull << low_bit)] + W[low_bit + 1];
    }

    const size_t num_rounds = m / K;
    const size_t residual   = m - K * num_rounds;

    size_t input_size = n;
    size_t n_modules  = 1;

    // ---------- Main: num_rounds rounds of radix-2^K ----------
    for (size_t round_j = 0; round_j < num_rounds; ++round_j) {
        const size_t L        = input_size;
        const size_t logL     = m - K * round_j;
        const size_t chunk    = L >> K;
        const size_t half     = L >> 1;
        const size_t num_bits = K * round_j;

        // Per-level nz_S and asr (shift accumulated through S^{logL - ell}).
        std::array<std::vector<size_t>, K + 1> nz_S;
        std::array<FieldT, K + 1>              asr;
        for (size_t ell = 1; ell <= K; ++ell) {
            const size_t S_idx = logL - ell;
            asr[ell] = affine_shift;
            if (S_idx > 0) {
                nz_S[ell].reserve(S_idx);
                nz_S[ell].emplace_back((1ull << S_idx) - 1);
                for (size_t i = 1; i < S_idx; ++i) {
                    size_t ii = i, rr = S_idx;
                    while (((rr & 1) | (~ii & 1)) && ii > 0) { ii >>= 1; rr >>= 1; }
                    if (ii == 0) {
                        nz_S[ell].emplace_back((1ull << S_idx) - (1ull << i));
                        asr[ell] += affine_shift ^ (1ull << i);
                    }
                }
                asr[ell] += affine_shift ^ (1ull << S_idx);
            }
        }

        size_t offset = 0;
        for (size_t module = 0; module < n_modules; ++module) {
            // Base twiddles base_T[ell] = S^{logL - ell}(theta_module).
            std::array<FieldT, K + 1> base_T;
            for (size_t ell = 1; ell <= K; ++ell)
                base_T[ell] = asr[ell];
            for (size_t i = 0; i < num_bits; ++i) {
                if (module & (1ull << i)) {
                    for (size_t ell = 1; ell <= K; ++ell)
                        base_T[ell] += W[i + ell];
                }
            }

            // ---- Sweep 1: level-1 division + stage-1 cross-half twiddle ----
            //   For k in [L/2, L): apply S^{logL-1} division and the T_{u,1,0}
            //   cross-half twiddle. After this sweep:
            //     g[0..L/2)    = h_0 = r + T_{u,1,0} * q (round-0 lower)
            //     g[L/2..L)    = q   (round-0 quotient, untouched here)
            //   This is exactly one pure radix-2 sweep-A.
            {
                const FieldT T1 = base_T[1];
                for (size_t k = L; k > half; ) {
                    --k;
                    const FieldT gk = g[offset + k];
                    for (size_t nz : nz_S[1])
                        g[offset + k - nz] += gk;
                    g[offset + k - half] += gk * T1;
                }
            }

            // ---- Sweeps 2..K: level-ell divisions ----
            //   At level ell there are 2^{ell-1} disjoint active bands of
            //   size L/2^ell each, total L/2 active positions. We iterate
            //   each band in standard descending order; the bands themselves
            //   are independent so we can visit them in any order. We pick
            //   ascending band order for a sequential memory pass.
            for (size_t ell = 2; ell <= K; ++ell) {
                const size_t span      = 1ull << (logL - ell);     // = L/2^ell
                const size_t step      = span << 1;                // = L/2^{ell-1}
                const size_t num_bands = 1ull << (ell - 1);
                for (size_t b = 0; b < num_bands; ++b) {
                    const size_t band_start = b * step + span;
                    const size_t band_end   = band_start + span;
                    for (size_t k = band_end; k > band_start; ) {
                        --k;
                        const FieldT gk = g[offset + k];
                        for (size_t nz : nz_S[ell])
                            g[offset + k - nz] += gk;
                    }
                }
            }

            // ---- Sweep K+1: in-register K-stage butterfly, per tuple ----
            for (size_t jj = 0; jj < chunk; ++jj) {
                FieldT v[Q];
                for (size_t i = 0; i < Q; ++i)
                    v[i] = g[offset + i * chunk + jj];

                // Stage-1 closer: v[Q/2 + i] += v[i].
                for (size_t i = 0; i < Q_half; ++i)
                    v[Q_half + i] += v[i];

                // Stages 2..K (no-op when K == 1).
                for (size_t ell = 2; ell <= K; ++ell) {
                    const size_t stride     = 1ull << (K - ell);
                    const size_t num_blocks = 1ull << (ell - 1);
                    for (size_t b = 0; b < num_blocks; ++b) {
                        const FieldT T_lb = base_T[ell] + block_offsets[b];
                        for (size_t t = 0; t < stride; ++t) {
                            const size_t i0 = b * (stride << 1) + t;
                            const size_t i1 = i0 + stride;
                            v[i0] += T_lb * v[i1];
                            v[i1] += v[i0];
                        }
                    }
                }

                for (size_t i = 0; i < Q; ++i)
                    g[offset + i * chunk + jj] = v[i];
            }

            offset += L;
        }

        input_size >>= K;
        n_modules  <<= K;
    }

    // ---------- Residual radix-2 rounds (mixed-radix epilogue) ----------
    for (size_t r_res = 0; r_res < residual; ++r_res) {
        const size_t round_idx = K * num_rounds + r_res;
        const size_t S_idx     = m - round_idx - 1;
        const size_t L         = input_size;
        const size_t half      = L >> 1;

        std::vector<size_t> nz_S;
        FieldT asr = affine_shift;
        if (S_idx > 0) {
            nz_S.reserve(S_idx);
            nz_S.emplace_back((1ull << S_idx) - 1);
            for (size_t i = 1; i < S_idx; ++i) {
                size_t ii = i, rr = S_idx;
                while (((rr & 1) | (~ii & 1)) && ii > 0) { ii >>= 1; rr >>= 1; }
                if (ii == 0) {
                    nz_S.emplace_back((1ull << S_idx) - (1ull << i));
                    asr += affine_shift ^ (1ull << i);
                }
            }
            asr += affine_shift ^ (1ull << S_idx);
        }

        size_t offset = 0;
        for (size_t module = 0; module < n_modules; ++module) {
            FieldT mult_factor = asr;
            for (size_t i = 0; i < round_idx; ++i) {
                if (module & (1ull << i))
                    mult_factor += W[i + 1];
            }
            const size_t offset2 = offset + half;
            for (size_t k = offset2 + half; k > offset2; ) {
                --k;
                FieldT gk = g[k];
                for (const auto& nz : nz_S)
                    g[k - nz] += gk;
                g[k - half] += gk * mult_factor;
            }
            for (size_t j = 0; j < half; ++j)
                g[offset2 + j] += g[offset + j];
            offset += L;
        }
        input_size >>= 1;
        n_modules  <<= 1;
    }

    return g;
}

template<typename FieldT>
std::vector<FieldT> additive_IFFT(const std::vector<FieldT> &evals,
                                        const libiop::affine_subspace<FieldT> &domain)
{
    std::vector<FieldT> g(evals);
    
    FieldT affine_shift = domain.shift();
    std::vector<FieldT> W(domain.basis());

    const size_t n = g.size();
    const size_t m = domain.dimension();
    assert(n == (1ull<<m));

    size_t input_size = 1;
    size_t n_modules = n;
    size_t S_index = 0;
    for (int r = m-1; r >=0; --r)
    {
        FieldT affine_shift_round = affine_shift;
        n_modules >>= 1;
        size_t half_input_size = input_size;
        input_size <<= 1;

        // Computing non-zero indices in S_{m-r}(X) for round r in the reverse order. Becuse in division algorithm we start from 
        std::vector<size_t> nz_S;

        if(S_index > 0){
            nz_S.reserve(S_index);
            nz_S.emplace_back((1<<S_index)-1); // C(x,0) = 1, so we assigned that before loop.
            for (size_t i = 1; i < S_index; ++i){
                size_t ii = i ;
                size_t rr = S_index;
                while (((rr & 1) | (~ii & 1)) && ii>0){
                    ii >>= 1;
                    rr >>= 1;
                }
                if (ii == 0){
                    nz_S.emplace_back((1<<S_index) - (1<<i));
                    affine_shift_round += affine_shift ^ (1<<i);
                }
            }
            affine_shift_round += affine_shift ^ (1<<S_index);
        }  
        size_t offset = 0;
        for (size_t module = 0; module < n_modules; ++module){
            // Computing the multiplication factor
            FieldT mult_factor = affine_shift_round;
            for (size_t i = 0; i < r; ++i){
                if (module & (1<<i))
                    mult_factor += W[i+1];
            }
            size_t offset2 = offset + half_input_size;
            for (size_t j = 0; j < half_input_size; ++j)
                g[offset2+j] += g[offset+j];    
            for (size_t k =  offset2; k < offset2+half_input_size; ++k){
                FieldT gk = g[k];
                g[k-half_input_size] += gk * mult_factor;
                for (const auto& nz : nz_S)
                    g[k - nz] += gk ;
            }
            offset += input_size;
        }
        S_index ++;
    }
    return g;
}

} // namespace cantor
