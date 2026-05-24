namespace cantor {

/// One radix-2 butterfly module for \ref additive_FFT(poly, domain_dim, shift_dim):
/// \c s_i / \c n_terms sparse indices and twiddle from \c cantor_combinations via
/// \c (module|shift_bit)<<1 (same inner loops as the serial table-path implementation).
template<typename FieldT>
void cantor_additive_fft_combination_radix2_module(std::vector<FieldT> &g,
                                                   const size_t offset,
                                                   const size_t half_input_size,
                                                   const size_t module,
                                                   const size_t m,
                                                   const size_t shift_dim,
                                                   const size_t n_modules,
                                                   const unsigned *nz_S,
                                                   const unsigned t,
                                                   FieldT *cantor_combinations)
{
    const size_t shift_bit      = shift_dim == 0 ? 0 : n_modules << (shift_dim - m);
    size_t       module_shifted = (module | shift_bit) << 1;
    size_t       i_256          = 0;
    FieldT       mult_factor    = FieldT::zero();
    while (module_shifted) {
        mult_factor += cantor_combinations[i_256 + (module_shifted & 0xff)];
        i_256 += 256;
        module_shifted >>= 8;
    }
    const size_t offset2 = offset + half_input_size;
    for (unsigned k = offset2 + half_input_size - 1; k >= offset2; --k) {
        FieldT gk = g[k];
        for (unsigned i = 0; i < t; ++i) {
            g[k - nz_S[i]] += gk;
        }
        g[k - half_input_size] += gk * mult_factor;
    }
    for (size_t j = 0; j < half_input_size; ++j) {
        g[offset2 + j] += g[offset + j];
    }
}

/// Decode Cantor-coordinate bit pattern into a field element (same 8-bit
/// chunked table layout as \ref additive_FFT(poly, domain_dim, shift_dim)).
template<typename FieldT>
static inline FieldT cantor_element_from_cantor_bits(size_t bits, FieldT *cantor_combinations)
{
    FieldT acc   = FieldT::zero();
    size_t i_256 = 0;
    while (bits) {
        acc += cantor_combinations[i_256 + (bits & 0xff)];
        i_256 += 256;
        bits >>= 8;
    }
    return acc;
}

/// One fused radix-2^K module for the Cantor-combination table path
/// (\c s_i / \c n_terms + \c cantor_combinations), mirroring
/// \ref cantor_additive_fft_radix2k_process_module but with LCH-style twiddles
/// \c base_T[ell] = element_from_cantor_bits((module<<1 | (shift_bit_1<<1)) << (ell-1)).
template<size_t K, typename FieldT>
void cantor_additive_fft_combination_radix2k_process_module(
    std::vector<FieldT> &g,
    const size_t offset,
    const size_t L,
    const size_t chunk,
    const size_t half,
    const size_t logL,
    const size_t module,
    const size_t m,
    const size_t shift_dim,
    const size_t n_modules,
    FieldT *cantor_combinations,
    const std::array<const unsigned *, K + 1> &nz_ptr,
    const std::array<unsigned, K + 1> &t_arr)
{
    constexpr size_t Q      = 1ull << K;
    constexpr size_t Q_half = Q >> 1;

    std::array<FieldT, (Q_half == 0 ? 1 : Q_half)> block_offsets;
    block_offsets[0] = FieldT::zero();
    for (size_t c = 1; c < Q_half; ++c) {
        const size_t low_bit = __builtin_ctzll(c);
        block_offsets[c] =
            block_offsets[c ^ (1ull << low_bit)] + cantor_combinations[1ull << (low_bit + 1)];
    }

    const size_t shift_bit_1 = shift_dim == 0 ? 0 : (n_modules << (shift_dim - m));
    const size_t mshift_1    = (module << 1) | (shift_bit_1 << 1);

    std::array<FieldT, K + 1> base_T;
    for (size_t ell = 1; ell <= K; ++ell)
        base_T[ell] = cantor_element_from_cantor_bits<FieldT>(mshift_1 << (ell - 1), cantor_combinations);

    {
        const FieldT    T1  = base_T[1];
        const unsigned *nz1 = nz_ptr[1];
        const unsigned  t1  = t_arr[1];
        for (size_t k = L; k > half;) {
            --k;
            const FieldT gk = g[offset + k];
            for (unsigned i = 0; i < t1; ++i)
                g[offset + k - nz1[i]] += gk;
            g[offset + k - half] += gk * T1;
        }
    }

    for (size_t ell = 2; ell <= K; ++ell) {
        const size_t    span      = 1ull << (logL - ell);
        const size_t    step      = span << 1;
        const size_t    num_bands = 1ull << (ell - 1);
        const unsigned *nzel      = nz_ptr[ell];
        const unsigned  te        = t_arr[ell];
        for (size_t b = 0; b < num_bands; ++b) {
            const size_t band_start = b * step + span;
            const size_t band_end   = band_start + span;
            for (size_t k = band_end; k > band_start;) {
                --k;
                const FieldT gk = g[offset + k];
                for (unsigned i = 0; i < te; ++i)
                    g[offset + k - nzel[i]] += gk;
            }
        }
    }

    for (size_t jj = 0; jj < chunk; ++jj) {
        FieldT v[Q];
        for (size_t i = 0; i < Q; ++i)
            v[i] = g[offset + i * chunk + jj];

        for (size_t i = 0; i < Q_half; ++i)
            v[Q_half + i] += v[i];

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
}

/// One radix-2 butterfly module for the affine-subspace chart
/// (\ref additive_FFT(poly, domain) / \ref additive_FFT_radix2k(poly, domain)).
template<typename FieldT>
void cantor_additive_fft_radix2_process_module(std::vector<FieldT> &g,
                                               const size_t offset,
                                               const size_t half_input_size,
                                               const size_t module,
                                               const size_t r,
                                               const FieldT &affine_shift_round,
                                               const std::vector<FieldT> &W,
                                               const std::vector<size_t> &nz_S)
{
    FieldT mult_factor = affine_shift_round;
    for (size_t i = 0; i < r; ++i) {
        if (module & (1ull << i))
            mult_factor += W[i + 1];
    }
    const size_t offset2 = offset + half_input_size;
    for (size_t k = offset2 + half_input_size - 1; k >= offset2; --k) {
        FieldT gk = g[k];
        for (const auto &nz : nz_S)
            g[k - nz] += gk;
        g[k - half_input_size] += gk * mult_factor;
    }
    for (size_t j = 0; j < half_input_size; ++j)
        g[offset2 + j] += g[offset + j];
}

/// One fused radix-2^K module for the affine-subspace chart
/// (\ref additive_FFT_radix2k(poly, domain)).
template<size_t K, typename FieldT>
void cantor_additive_fft_radix2k_process_module(
    std::vector<FieldT> &g,
    const size_t offset,
    const size_t L,
    const size_t chunk,
    const size_t half,
    const size_t logL,
    const size_t num_bits,
    const size_t module,
    const std::array<std::vector<size_t>, K + 1> &nz_S,
    const std::array<FieldT, K + 1> &asr,
    const std::vector<FieldT> &W,
    const std::array<FieldT, (((1ull << K) >> 1) == 0 ? 1 : ((1ull << K) >> 1))> &block_offsets)
{
    constexpr size_t Q      = 1ull << K;
    constexpr size_t Q_half = Q >> 1;

    std::array<FieldT, K + 1> base_T;
    for (size_t ell = 1; ell <= K; ++ell)
        base_T[ell] = asr[ell];
    for (size_t i = 0; i < num_bits; ++i) {
        if (module & (1ull << i)) {
            for (size_t ell = 1; ell <= K; ++ell)
                base_T[ell] += W[i + ell];
        }
    }

    {
        const FieldT T1 = base_T[1];
        for (size_t k = L; k > half;) {
            --k;
            const FieldT gk = g[offset + k];
            for (size_t nz : nz_S[1])
                g[offset + k - nz] += gk;
            g[offset + k - half] += gk * T1;
        }
    }

    for (size_t ell = 2; ell <= K; ++ell) {
        const size_t span      = 1ull << (logL - ell);
        const size_t step      = span << 1;
        const size_t num_bands = 1ull << (ell - 1);
        for (size_t b = 0; b < num_bands; ++b) {
            const size_t band_start = b * step + span;
            const size_t band_end   = band_start + span;
            for (size_t k = band_end; k > band_start;) {
                --k;
                const FieldT gk = g[offset + k];
                for (size_t nz : nz_S[ell])
                    g[offset + k - nz] += gk;
            }
        }
    }

    for (size_t jj = 0; jj < chunk; ++jj) {
        FieldT v[Q];
        for (size_t i = 0; i < Q; ++i)
            v[i] = g[offset + i * chunk + jj];

        for (size_t i = 0; i < Q_half; ++i)
            v[Q_half + i] += v[i];

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
}

/// Inverse of \ref cantor_additive_fft_combination_radix2_module (radix-2 table path).
template<typename FieldT>
void cantor_additive_ifft_combination_radix2_module(std::vector<FieldT> &g,
                                                  const size_t offset,
                                                  const size_t half_input_size,
                                                  const size_t module,
                                                  const size_t m,
                                                  const size_t shift_dim,
                                                  const size_t n_modules,
                                                  const unsigned *nz_S,
                                                  const unsigned t,
                                                  FieldT *cantor_combinations)
{
    const size_t offset2 = offset + half_input_size;
    for (size_t j = 0; j < half_input_size; ++j)
        g[offset2 + j] += g[offset + j];

    const size_t shift_bit      = shift_dim == 0 ? 0 : n_modules << (shift_dim - m);
    size_t       module_shifted = (module | shift_bit) << 1;
    size_t       i_256          = 0;
    FieldT       mult_factor    = FieldT::zero();
    while (module_shifted) {
        mult_factor += cantor_combinations[i_256 + (module_shifted & 0xff)];
        i_256 += 256;
        module_shifted >>= 8;
    }
    for (size_t k = offset2; k < offset2 + half_input_size; ++k) {
        FieldT gk = g[k];
        g[k - half_input_size] += gk * mult_factor;
        for (unsigned i = 0; i < t; ++i)
            g[k - nz_S[i]] += gk;
    }
}

/// Inverse of \ref cantor_additive_fft_combination_radix2k_process_module.
template<size_t K, typename FieldT>
void cantor_additive_ifft_combination_radix2k_process_module(
    std::vector<FieldT> &g,
    const size_t offset,
    const size_t L,
    const size_t chunk,
    const size_t half,
    const size_t logL,
    const size_t module,
    const size_t m,
    const size_t shift_dim,
    const size_t n_modules,
    FieldT *cantor_combinations,
    const std::array<const unsigned *, K + 1> &nz_ptr,
    const std::array<unsigned, K + 1> &t_arr)
{
    constexpr size_t Q      = 1ull << K;
    constexpr size_t Q_half = Q >> 1;

    std::array<FieldT, (Q_half == 0 ? 1 : Q_half)> block_offsets;
    block_offsets[0] = FieldT::zero();
    for (size_t c = 1; c < Q_half; ++c) {
        const size_t low_bit = __builtin_ctzll(c);
        block_offsets[c] =
            block_offsets[c ^ (1ull << low_bit)] + cantor_combinations[1ull << (low_bit + 1)];
    }

    const size_t shift_bit_1 = shift_dim == 0 ? 0 : (n_modules << (shift_dim - m));
    const size_t mshift_1    = (module << 1) | (shift_bit_1 << 1);

    std::array<FieldT, K + 1> base_T;
    for (size_t ell = 1; ell <= K; ++ell)
        base_T[ell] = cantor_element_from_cantor_bits<FieldT>(mshift_1 << (ell - 1), cantor_combinations);

    for (size_t jj = 0; jj < chunk; ++jj) {
        FieldT v[Q];
        for (size_t i = 0; i < Q; ++i)
            v[i] = g[offset + i * chunk + jj];

        for (size_t ell = K; ell >= 2; --ell) {
            const size_t stride     = 1ull << (K - ell);
            const size_t num_blocks = 1ull << (ell - 1);
            for (size_t b = num_blocks; b-- > 0;) {
                const FieldT T_lb = base_T[ell] + block_offsets[b];
                for (size_t t = stride; t-- > 0;) {
                    const size_t i0 = b * (stride << 1) + t;
                    const size_t i1 = i0 + stride;
                    v[i1] += v[i0];
                    v[i0] += T_lb * v[i1];
                }
            }
        }

        for (size_t i = 0; i < Q_half; ++i)
            v[Q_half + i] += v[i];

        for (size_t i = 0; i < Q; ++i)
            g[offset + i * chunk + jj] = v[i];
    }

    for (size_t ell = K; ell >= 2; --ell) {
        const size_t    span      = 1ull << (logL - ell);
        const size_t    step      = span << 1;
        const size_t    num_bands = 1ull << (ell - 1);
        const unsigned *nzel      = nz_ptr[ell];
        const unsigned  te        = t_arr[ell];
        for (size_t b = 0; b < num_bands; ++b) {
            const size_t band_start = b * step + span;
            const size_t band_end   = band_start + span;
            for (size_t k = band_start; k < band_end; ++k) {
                FieldT gk = g[offset + k];
                for (unsigned i = 0; i < te; ++i)
                    g[offset + k - nzel[i]] += gk;
            }
        }
    }

    {
        const FieldT    T1  = base_T[1];
        const unsigned *nz1 = nz_ptr[1];
        const unsigned  t1  = t_arr[1];
        for (size_t k = half; k < L; ++k) {
            FieldT gk = g[offset + k];
            g[offset + k - half] += gk * T1;
            for (unsigned i = 0; i < t1; ++i)
                g[offset + k - nz1[i]] += gk;
        }
    }
}

/// Inverse of \ref cantor_additive_fft_radix2_process_module (affine-subspace chart).
template<typename FieldT>
void cantor_additive_ifft_radix2_process_module(std::vector<FieldT> &g,
                                              const size_t offset,
                                              const size_t half_input_size,
                                              const size_t module,
                                              const size_t r,
                                              const FieldT &affine_shift_round,
                                              const std::vector<FieldT> &W,
                                              const std::vector<size_t> &nz_S)
{
    FieldT mult_factor = affine_shift_round;
    for (size_t i = 0; i < r; ++i) {
        if (module & (1ull << i))
            mult_factor += W[i + 1];
    }
    const size_t offset2 = offset + half_input_size;
    for (size_t j = 0; j < half_input_size; ++j)
        g[offset2 + j] += g[offset + j];

    for (size_t k = offset2; k < offset2 + half_input_size; ++k) {
        FieldT gk = g[k];
        g[k - half_input_size] += gk * mult_factor;
        for (const auto &nz : nz_S)
            g[k - nz] += gk;
    }
}

/// Inverse of \ref cantor_additive_fft_radix2k_process_module (affine-subspace chart).
template<size_t K, typename FieldT>
void cantor_additive_ifft_radix2k_process_module(
    std::vector<FieldT> &g,
    const size_t offset,
    const size_t L,
    const size_t chunk,
    const size_t half,
    const size_t logL,
    const size_t num_bits,
    const size_t module,
    const std::array<std::vector<size_t>, K + 1> &nz_S,
    const std::array<FieldT, K + 1> &asr,
    const std::vector<FieldT> &W,
    const std::array<FieldT, (((1ull << K) >> 1) == 0 ? 1 : ((1ull << K) >> 1))> &block_offsets)
{
    constexpr size_t Q      = 1ull << K;
    constexpr size_t Q_half = Q >> 1;

    std::array<FieldT, K + 1> base_T;
    for (size_t ell = 1; ell <= K; ++ell)
        base_T[ell] = asr[ell];
    for (size_t i = 0; i < num_bits; ++i) {
        if (module & (1ull << i)) {
            for (size_t ell = 1; ell <= K; ++ell)
                base_T[ell] += W[i + ell];
        }
    }

    for (size_t jj = 0; jj < chunk; ++jj) {
        FieldT v[Q];
        for (size_t i = 0; i < Q; ++i)
            v[i] = g[offset + i * chunk + jj];

        for (size_t ell = K; ell >= 2; --ell) {
            const size_t stride     = 1ull << (K - ell);
            const size_t num_blocks = 1ull << (ell - 1);
            for (size_t b = num_blocks; b-- > 0;) {
                const FieldT T_lb = base_T[ell] + block_offsets[b];
                for (size_t t = stride; t-- > 0;) {
                    const size_t i0 = b * (stride << 1) + t;
                    const size_t i1 = i0 + stride;
                    v[i1] += v[i0];
                    v[i0] += T_lb * v[i1];
                }
            }
        }

        for (size_t i = 0; i < Q_half; ++i)
            v[Q_half + i] += v[i];

        for (size_t i = 0; i < Q; ++i)
            g[offset + i * chunk + jj] = v[i];
    }

    for (size_t ell = K; ell >= 2; --ell) {
        const size_t span      = 1ull << (logL - ell);
        const size_t step      = span << 1;
        const size_t num_bands = 1ull << (ell - 1);
        for (size_t b = 0; b < num_bands; ++b) {
            const size_t band_start = b * step + span;
            const size_t band_end   = band_start + span;
            for (size_t k = band_start; k < band_end; ++k) {
                FieldT gk = g[offset + k];
                for (size_t nz : nz_S[ell])
                    g[offset + k - nz] += gk;
            }
        }
    }

    {
        const FieldT T1 = base_T[1];
        for (size_t k = half; k < L; ++k) {
            FieldT gk = g[offset + k];
            g[offset + k - half] += gk * T1;
            for (size_t nz : nz_S[1])
                g[offset + k - nz] += gk;
        }
    }
}

} // namespace cantor
