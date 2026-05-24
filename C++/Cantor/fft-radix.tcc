
namespace cantor{

    /// Fused radix-2^K Cantor FFT on the \b combination-table path
    /// (\c s_i / \c n_terms / \c cantor_combinations), same result as
    /// \ref additive_FFT(poly, domain_dim, shift_dim). Twiddle bits follow the same
    /// rule as LCH \c butterfly_radix2k (see \c cantor_additive_fft_combination_radix2k_process_module).
    template<size_t K, typename FieldT>
    std::vector<FieldT> additive_FFT_radix2k(const std::vector<FieldT> &poly_coeffs,
                                             const size_t domain_dim,
                                             const size_t shift_dim)
    {
        static_assert(K >= 1, "K must be >= 1");
        static_assert(K <= 5, "K supported up to 5 (radix up to 32)");

        const size_t m = domain_dim;
        std::vector<FieldT> g(poly_coeffs);
        g.resize((1ULL << m), FieldT::zero());
        const size_t n = g.size();
        assert(n == (1ull << m));

        FieldT *cantor_combinations;
        if (FieldT::extension_degree() == 128)
            cantor_combinations = (FieldT *)cantor_combinations_8R_in_gf2to128;
        else if (FieldT::extension_degree() == 192)
            cantor_combinations = (FieldT *)cantor_combinations_8R_in_gf2to192;
        else if (FieldT::extension_degree() == 256)
            cantor_combinations = (FieldT *)cantor_combinations_8R_in_gf2to256;
        else
            throw std::invalid_argument(
                "The field size should be either 128, 192, or 256 for using the cantor basis");

        const size_t num_rounds = m / K;
        const size_t residual   = m - K * num_rounds;

        size_t input_size = n;
        size_t n_modules  = 1;

        for (size_t round_j = 0; round_j < num_rounds; ++round_j) {
            const size_t L     = input_size;
            const size_t logL  = m - K * round_j;
            const size_t chunk = L >> K;
            const size_t half  = L >> 1;

            std::array<const unsigned *, K + 1> nz_ptr{};
            std::array<unsigned, K + 1>        t_arr{};
            for (size_t ell = 1; ell <= K; ++ell) {
                const size_t S_idx = logL - ell;
                nz_ptr[ell]        = s_i[S_idx];
                t_arr[ell]         = n_terms[S_idx];
            }

            for (size_t module = 0; module < n_modules; ++module) {
                const size_t offset = module * L;
                cantor_additive_fft_combination_radix2k_process_module<K, FieldT>(
                    g, offset, L, chunk, half, logL, module, m, shift_dim, n_modules, cantor_combinations,
                    nz_ptr, t_arr);
            }

            input_size >>= K;
            n_modules <<= K;
        }

        for (size_t r_res = 0; r_res < residual; ++r_res) {
            const size_t round_idx = K * num_rounds + r_res;
            const size_t S_index   = m - round_idx - 1;
            const unsigned *nz_S   = s_i[S_index];
            const unsigned  t      = n_terms[S_index];
            const size_t L         = input_size;
            const size_t half      = L >> 1;

            for (size_t module = 0; module < n_modules; ++module) {
                const size_t offset = module * L;
                cantor_additive_fft_combination_radix2_module<FieldT>(
                    g, offset, half, module, m, shift_dim, n_modules, nz_S, t, cantor_combinations);
            }
            input_size >>= 1;
            n_modules <<= 1;
        }

        return g;
    }

    // Templated radix-2^K Cantor AFFT (K=2 is radix-4).
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

        constexpr size_t Q      = 1ull << K;
        constexpr size_t Q_half = Q >> 1;

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

        for (size_t round_j = 0; round_j < num_rounds; ++round_j) {
            const size_t L        = input_size;
            const size_t logL     = m - K * round_j;
            const size_t chunk    = L >> K;
            const size_t half     = L >> 1;
            const size_t num_bits = K * round_j;

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

            for (size_t module = 0; module < n_modules; ++module) {
                cantor_additive_fft_radix2k_process_module<K, FieldT>(
                    g, module * L, L, chunk, half, logL, num_bits, module, nz_S, asr, W, block_offsets);
            }

            input_size >>= K;
            n_modules  <<= K;
        }

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

            for (size_t module = 0; module < n_modules; ++module) {
                cantor_additive_fft_radix2_process_module(g, module * L, half, module, round_idx, asr, W,
                                                          nz_S);
            }
            input_size >>= 1;
            n_modules  <<= 1;
        }

        return g;
    }

    template<size_t K, typename FieldT>
    std::vector<FieldT> additive_IFFT_radix2k(const std::vector<FieldT> &evals,
                                              const size_t domain_dim,
                                              const size_t shift_dim)
    {
        static_assert(K >= 1, "K must be >= 1");
        static_assert(K <= 5, "K supported up to 5 (radix up to 32)");

        const size_t        m = domain_dim;
        std::vector<FieldT> g(evals);
        const size_t        n = g.size();
        assert(n == (1ull << m));

        FieldT *cantor_combinations;
        if (FieldT::extension_degree() == 128)
            cantor_combinations = (FieldT *)cantor_combinations_8R_in_gf2to128;
        else if (FieldT::extension_degree() == 192)
            cantor_combinations = (FieldT *)cantor_combinations_8R_in_gf2to192;
        else if (FieldT::extension_degree() == 256)
            cantor_combinations = (FieldT *)cantor_combinations_8R_in_gf2to256;
        else
            throw std::invalid_argument(
                "The field size should be either 128, 192, or 256 for using the cantor basis");

        const size_t num_rounds = m / K;
        const size_t residual   = m - K * num_rounds;

        size_t input_size = 1;
        size_t n_modules  = n;

        for (size_t r_res = residual; r_res-- > 0;) {
            const size_t round_idx = K * num_rounds + r_res;
            n_modules >>= 1;
            input_size <<= 1;
            const size_t    S_index = m - round_idx - 1;
            const unsigned *nz_S    = s_i[S_index];
            const unsigned  t       = n_terms[S_index];
            const size_t L          = input_size;
            const size_t half       = L >> 1;

            for (size_t module = 0; module < n_modules; ++module) {
                cantor_additive_ifft_combination_radix2_module<FieldT>(
                    g, module * L, half, module, m, shift_dim, n_modules, nz_S, t, cantor_combinations);
            }
        }

        for (size_t round_j = num_rounds; round_j-- > 0;) {
            n_modules >>= K;
            input_size <<= K;
            const size_t L     = input_size;
            const size_t logL  = m - K * round_j;
            const size_t chunk = L >> K;
            const size_t half  = L >> 1;

            std::array<const unsigned *, K + 1> nz_ptr{};
            std::array<unsigned, K + 1>        t_arr{};
            for (size_t ell = 1; ell <= K; ++ell) {
                const size_t S_idx = logL - ell;
                nz_ptr[ell]        = s_i[S_idx];
                t_arr[ell]         = n_terms[S_idx];
            }

            for (size_t module = 0; module < n_modules; ++module) {
                cantor_additive_ifft_combination_radix2k_process_module<K, FieldT>(
                    g, module * L, L, chunk, half, logL, module, m, shift_dim, n_modules, cantor_combinations,
                    nz_ptr, t_arr);
            }
        }

        return g;
    }

    template<size_t K, typename FieldT>
    std::vector<FieldT> additive_IFFT_radix2k(const std::vector<FieldT> &evals,
                                              const libiop::affine_subspace<FieldT> &domain)
    {
        static_assert(K >= 1, "K must be >= 1");
        static_assert(K <= 5, "K supported up to 5 (radix up to 32)");

        std::vector<FieldT> g(evals);
        g.resize(domain.num_elements(), FieldT::zero());

        const FieldT              affine_shift = domain.shift();
        const std::vector<FieldT> W(domain.basis());

        const size_t n = g.size();
        const size_t m = domain.dimension();
        assert(n == (1ull << m));

        constexpr size_t Q      = 1ull << K;
        constexpr size_t Q_half = Q >> 1;

        std::array<FieldT, (Q_half == 0 ? 1 : Q_half)> block_offsets;
        block_offsets[0] = FieldT::zero();
        for (size_t c = 1; c < Q_half; ++c) {
            const size_t low_bit = __builtin_ctzll(c);
            block_offsets[c] = block_offsets[c ^ (1ull << low_bit)] + W[low_bit + 1];
        }

        const size_t num_rounds = m / K;
        const size_t residual   = m - K * num_rounds;

        size_t input_size = 1;
        size_t n_modules  = n;

        for (size_t r_res = residual; r_res-- > 0;) {
            const size_t round_idx = K * num_rounds + r_res;
            n_modules >>= 1;
            input_size <<= 1;
            const size_t S_idx = m - round_idx - 1;
            const size_t L     = input_size;
            const size_t half  = L >> 1;

            std::vector<size_t> nz_S;
            FieldT              asr = affine_shift;
            if (S_idx > 0) {
                nz_S.reserve(S_idx);
                nz_S.emplace_back((1ull << S_idx) - 1);
                for (size_t i = 1; i < S_idx; ++i) {
                    size_t ii = i, rr = S_idx;
                    while (((rr & 1) | (~ii & 1)) && ii > 0) {
                        ii >>= 1;
                        rr >>= 1;
                    }
                    if (ii == 0) {
                        nz_S.emplace_back((1ull << S_idx) - (1ull << i));
                        asr += affine_shift ^ (1ull << i);
                    }
                }
                asr += affine_shift ^ (1ull << S_idx);
            }

            for (size_t module = 0; module < n_modules; ++module) {
                cantor_additive_ifft_radix2_process_module(g, module * L, half, module, round_idx, asr, W,
                                                           nz_S);
            }
        }

        for (size_t round_j = num_rounds; round_j-- > 0;) {
            n_modules >>= K;
            input_size <<= K;
            const size_t L        = input_size;
            const size_t logL     = m - K * round_j;
            const size_t chunk    = L >> K;
            const size_t half     = L >> 1;
            const size_t num_bits = K * round_j;

            std::array<std::vector<size_t>, K + 1> nz_S;
            std::array<FieldT, K + 1>              asr;
            for (size_t ell = 1; ell <= K; ++ell) {
                const size_t S_idx = logL - ell;
                asr[ell]           = affine_shift;
                if (S_idx > 0) {
                    nz_S[ell].reserve(S_idx);
                    nz_S[ell].emplace_back((1ull << S_idx) - 1);
                    for (size_t i = 1; i < S_idx; ++i) {
                        size_t ii = i, rr = S_idx;
                        while (((rr & 1) | (~ii & 1)) && ii > 0) {
                            ii >>= 1;
                            rr >>= 1;
                        }
                        if (ii == 0) {
                            nz_S[ell].emplace_back((1ull << S_idx) - (1ull << i));
                            asr[ell] += affine_shift ^ (1ull << i);
                        }
                    }
                    asr[ell] += affine_shift ^ (1ull << S_idx);
                }
            }

            for (size_t module = 0; module < n_modules; ++module) {
                cantor_additive_ifft_radix2k_process_module<K, FieldT>(
                    g, module * L, L, chunk, half, logL, num_bits, module, nz_S, asr, W, block_offsets);
            }
        }

        return g;
    }

} // namespace cantor
