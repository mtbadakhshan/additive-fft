    
namespace cantor{    
    /// OpenMP parallel Cantor-combination table path
    /// (\c s_i / \c n_terms, same as \ref additive_FFT(poly, domain_dim, shift_dim)).
    /// Preconditions match the serial routine: for \c shift_dim != 0, use \c shift_dim >= domain_dim
    /// so that \c (shift_dim - m) is well-defined for \c n_modules << (...); otherwise behaviour
    /// matches the legacy serial code (which has the same constraint).
    template<typename FieldT>
    std::vector<FieldT> additive_FFT_parallel(const std::vector<FieldT> &poly_coeffs,
                                              const size_t domain_dim,
                                              const size_t shift_dim)
    {
    #ifndef _OPENMP
        return additive_FFT<FieldT>(poly_coeffs, domain_dim, shift_dim);
    #else
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
    
        size_t        input_size = n;
        size_t        n_modules  = 1;
        size_t        S_index    = m;
        const unsigned *nz_S     = nullptr;
        unsigned      t          = 0;
        size_t        half_input_size = 0;
    
    #pragma omp parallel shared(g, m, input_size, n_modules, S_index, nz_S, t, half_input_size, shift_dim,         \
                                    cantor_combinations)
        {
            for (size_t r = 0; r < m; ++r) {
    #pragma omp single
                {
                    --S_index;
                    nz_S            = s_i[S_index];
                    t               = n_terms[S_index];
                    half_input_size = input_size >> 1;
                }
    #pragma omp barrier
    
                if (n_modules >= 2) {
    #pragma omp for schedule(static)
                    for (size_t module = 0; module < n_modules; ++module) {
                        const size_t offset = module * input_size;
                        cantor_additive_fft_combination_radix2_module<FieldT>(
                            g, offset, half_input_size, module, m, shift_dim, n_modules, nz_S, t,
                            cantor_combinations);
                    }
                } else {
    #pragma omp single
                    {
                        cantor_additive_fft_combination_radix2_module<FieldT>(
                            g, 0, half_input_size, 0, m, shift_dim, n_modules, nz_S, t, cantor_combinations);
                    }
                }
    
    #pragma omp barrier
    #pragma omp single
                {
                    input_size = half_input_size;
                    n_modules <<= 1;
                }
    #pragma omp barrier
            }
        }
    
        return g;
    #endif
    }


template<typename FieldT>
    std::vector<FieldT> additive_FFT_parallel(const std::vector<FieldT> &poly_coeffs,
                                              const libiop::affine_subspace<FieldT> &domain)
    {
    #ifndef _OPENMP
        return additive_FFT(poly_coeffs, domain);
    #else
        std::vector<FieldT> g(poly_coeffs);
        g.resize(domain.num_elements(), FieldT::zero());
    
        const FieldT          affine_shift = domain.shift();
        const std::vector<FieldT> W(domain.basis());
    
        const size_t n = g.size();
        const size_t m = domain.dimension();
        assert(n == (1ull << m));
    
        size_t              input_size         = n;
        size_t              n_modules          = 1;
        size_t              S_index            = m;
        FieldT              affine_shift_round = FieldT::zero();
        size_t              half_input_size    = 0;
        std::vector<size_t> nz_S;
    
    #pragma omp parallel shared(g, affine_shift, W, m, input_size, n_modules, S_index, affine_shift_round, \
                                    half_input_size, nz_S)
        {
            for (size_t r = 0; r < m; ++r) {
    #pragma omp single
                {
                    affine_shift_round = affine_shift;
                    --S_index;
                    nz_S.clear();
                    if (S_index > 0) {
                        nz_S.reserve(S_index);
                        nz_S.emplace_back((1ull << S_index) - 1);
                        for (size_t i = 1; i < S_index; ++i) {
                            size_t ii = i;
                            size_t rr = S_index;
                            while (((rr & 1) | (~ii & 1)) && ii > 0) {
                                ii >>= 1;
                                rr >>= 1;
                            }
                            if (ii == 0) {
                                nz_S.emplace_back((1ull << S_index) - (1ull << i));
                                affine_shift_round += affine_shift ^ (1ull << i);
                            }
                        }
                        affine_shift_round += affine_shift ^ (1ull << S_index);
                    }
                    half_input_size = input_size >> 1;
                }
    #pragma omp barrier
    
                if (n_modules >= 2) {
    #pragma omp for schedule(static)
                    for (size_t module = 0; module < n_modules; ++module) {
                        const size_t offset = module * input_size;
                        cantor_additive_fft_radix2_process_module(g, offset, half_input_size, module,
                                                                  r, affine_shift_round, W, nz_S);
                    }
                } else {
    #pragma omp single
                    {
                        cantor_additive_fft_radix2_process_module(g, 0, half_input_size, 0, r,
                                                                  affine_shift_round, W, nz_S);
                    }
                }
    
    #pragma omp single
                {
                    input_size = half_input_size;
                    n_modules <<= 1;
                }
    #pragma omp barrier
            }
        }
    
        return g;
    #endif
    }

    template<size_t K, typename FieldT>
    std::vector<FieldT> additive_FFT_radix2k_parallel(const std::vector<FieldT> &poly_coeffs,
                                                      const size_t domain_dim,
                                                      const size_t shift_dim)
    {
    #ifndef _OPENMP
        return additive_FFT_radix2k<K, FieldT>(poly_coeffs, domain_dim, shift_dim);
    #else
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
    
        std::array<const unsigned *, K + 1> nz_ptr{};
        std::array<unsigned, K + 1>        t_arr{};
        size_t                             L     = 0;
        size_t                             logL  = 0;
        size_t                             chunk = 0;
        size_t                             half  = 0;
    
        size_t       round_idx = 0;
        const unsigned *nz_res = nullptr;
        unsigned     t_res     = 0;
        size_t       L_res     = 0;
        size_t       half_res  = 0;
    
    #pragma omp parallel shared(g, m, shift_dim, cantor_combinations, num_rounds, residual, input_size, n_modules, \
                                    nz_ptr, t_arr, L, logL, chunk, half, round_idx, nz_res, t_res, L_res, half_res)
        {
            for (size_t round_j = 0; round_j < num_rounds; ++round_j) {
    #pragma omp single
                {
                    L     = input_size;
                    logL  = m - K * round_j;
                    chunk = L >> K;
                    half  = L >> 1;
                    for (size_t ell = 1; ell <= K; ++ell) {
                        const size_t S_idx = logL - ell;
                        nz_ptr[ell]        = s_i[S_idx];
                        t_arr[ell]         = n_terms[S_idx];
                    }
                }
    #pragma omp barrier
    
                if (n_modules >= 2) {
    #pragma omp for schedule(static)
                    for (size_t module = 0; module < n_modules; ++module) {
                        cantor_additive_fft_combination_radix2k_process_module<K, FieldT>(
                            g, module * L, L, chunk, half, logL, module, m, shift_dim, n_modules,
                            cantor_combinations, nz_ptr, t_arr);
                    }
                } else {
    #pragma omp single
                    {
                        cantor_additive_fft_combination_radix2k_process_module<K, FieldT>(
                            g, 0, L, chunk, half, logL, 0, m, shift_dim, n_modules, cantor_combinations, nz_ptr,
                            t_arr);
                    }
                }
    
    #pragma omp single
                {
                    input_size >>= K;
                    n_modules <<= K;
                }
    #pragma omp barrier
            }
    
            for (size_t r_res = 0; r_res < residual; ++r_res) {
    #pragma omp single
                {
                    round_idx = K * num_rounds + r_res;
                    const size_t S_index = m - round_idx - 1;
                    nz_res               = s_i[S_index];
                    t_res                = n_terms[S_index];
                    L_res                = input_size;
                    half_res             = L_res >> 1;
                }
    #pragma omp barrier
    
                if (n_modules >= 2) {
    #pragma omp for schedule(static)
                    for (size_t module = 0; module < n_modules; ++module) {
                        cantor_additive_fft_combination_radix2_module<FieldT>(
                            g, module * L_res, half_res, module, m, shift_dim, n_modules, nz_res, t_res,
                            cantor_combinations);
                    }
                } else {
    #pragma omp single
                    {
                        cantor_additive_fft_combination_radix2_module<FieldT>(
                            g, 0, half_res, 0, m, shift_dim, n_modules, nz_res, t_res, cantor_combinations);
                    }
                }
    
    #pragma omp single
                {
                    input_size >>= 1;
                    n_modules <<= 1;
                }
    #pragma omp barrier
            }
        }
    
        return g;
    #endif
    }

    template<size_t K, typename FieldT>
    std::vector<FieldT> additive_FFT_radix2k_parallel(const std::vector<FieldT> &poly_coeffs,
                                                        const libiop::affine_subspace<FieldT> &domain)
    {
    #ifndef _OPENMP
        return additive_FFT_radix2k<K, FieldT>(poly_coeffs, domain);
    #else
        static_assert(K >= 1, "K must be >= 1");
        static_assert(K <= 5, "K supported up to 5 (radix up to 32)");
    
        std::vector<FieldT> g(poly_coeffs);
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
    
        size_t input_size = n;
        size_t n_modules  = 1;
    
        std::array<std::vector<size_t>, K + 1> nz_S_main;
        std::array<FieldT, K + 1>              asr_main;
        size_t                                 L        = 0;
        size_t                                 logL     = 0;
        size_t                                 chunk    = 0;
        size_t                                 half     = 0;
        size_t                                 num_bits = 0;
    
        std::vector<size_t> nz_S_res;
        FieldT              asr_res   = FieldT::zero();
        size_t              L_res     = 0;
        size_t              half_res  = 0;
        size_t              round_idx = 0;
    
    #pragma omp parallel shared(g, affine_shift, W, m, num_rounds, residual, input_size, n_modules, block_offsets, \
                                    nz_S_main, asr_main, L, logL, chunk, half, num_bits, nz_S_res, asr_res, L_res,      \
                                    half_res, round_idx)
        {
            for (size_t round_j = 0; round_j < num_rounds; ++round_j) {
    #pragma omp single
                {
                    L        = input_size;
                    logL     = m - K * round_j;
                    chunk    = L >> K;
                    half     = L >> 1;
                    num_bits = K * round_j;
                    for (size_t ell = 1; ell <= K; ++ell) {
                        const size_t S_idx = logL - ell;
                        asr_main[ell] = affine_shift;
                        nz_S_main[ell].clear();
                        if (S_idx > 0) {
                            nz_S_main[ell].reserve(S_idx);
                            nz_S_main[ell].emplace_back((1ull << S_idx) - 1);
                            for (size_t i = 1; i < S_idx; ++i) {
                                size_t ii = i, rr = S_idx;
                                while (((rr & 1) | (~ii & 1)) && ii > 0) {
                                    ii >>= 1;
                                    rr >>= 1;
                                }
                                if (ii == 0) {
                                    nz_S_main[ell].emplace_back((1ull << S_idx) - (1ull << i));
                                    asr_main[ell] += affine_shift ^ (1ull << i);
                                }
                            }
                            asr_main[ell] += affine_shift ^ (1ull << S_idx);
                        }
                    }
                }
    #pragma omp barrier
    
                if (n_modules >= 2) {
    #pragma omp for schedule(static)
                    for (size_t module = 0; module < n_modules; ++module) {
                        cantor_additive_fft_radix2k_process_module<K>(
                            g, module * L, L, chunk, half, logL, num_bits, module, nz_S_main, asr_main, W,
                            block_offsets);
                    }
                } else {
    #pragma omp single
                    {
                        cantor_additive_fft_radix2k_process_module<K>(
                            g, 0, L, chunk, half, logL, num_bits, 0, nz_S_main, asr_main, W, block_offsets);
                    }
                }
    
    #pragma omp single
                {
                    input_size >>= K;
                    n_modules <<= K;
                }
    #pragma omp barrier
            }
    
            for (size_t r_res = 0; r_res < residual; ++r_res) {
    #pragma omp single
                {
                    round_idx = K * num_rounds + r_res;
                    const size_t S_idx = m - round_idx - 1;
                    L_res    = input_size;
                    half_res = L_res >> 1;
                    nz_S_res.clear();
                    asr_res = affine_shift;
                    if (S_idx > 0) {
                        nz_S_res.reserve(S_idx);
                        nz_S_res.emplace_back((1ull << S_idx) - 1);
                        for (size_t i = 1; i < S_idx; ++i) {
                            size_t ii = i, rr = S_idx;
                            while (((rr & 1) | (~ii & 1)) && ii > 0) {
                                ii >>= 1;
                                rr >>= 1;
                            }
                            if (ii == 0) {
                                nz_S_res.emplace_back((1ull << S_idx) - (1ull << i));
                                asr_res += affine_shift ^ (1ull << i);
                            }
                        }
                        asr_res += affine_shift ^ (1ull << S_idx);
                    }
                }
    #pragma omp barrier
    
                if (n_modules >= 2) {
    #pragma omp for schedule(static)
                    for (size_t module = 0; module < n_modules; ++module) {
                        cantor_additive_fft_radix2_process_module(g, module * L_res, half_res, module,
                                                                  round_idx, asr_res, W, nz_S_res);
                    }
                } else {
    #pragma omp single
                    {
                        cantor_additive_fft_radix2_process_module(g, 0, half_res, 0, round_idx, asr_res,
                                                                  W, nz_S_res);
                    }
                }
    
    #pragma omp single
                {
                    input_size >>= 1;
                    n_modules <<= 1;
                }
    #pragma omp barrier
            }
        }
    
        return g;
    #endif
    }

    template<typename FieldT>
    std::vector<FieldT> additive_IFFT_parallel(const std::vector<FieldT> &evals,
                                               const size_t domain_dim,
                                               const size_t shift_dim)
    {
    #ifndef _OPENMP
        return additive_IFFT<FieldT>(evals, domain_dim, shift_dim);
    #else
        const size_t m = domain_dim;
        std::vector<FieldT> g(evals);
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

        size_t        input_size      = 1;
        size_t        n_modules       = n;
        size_t        S_index         = 0;
        const unsigned *nz_S          = nullptr;
        unsigned      t               = 0;
        size_t        half_input_size = 0;

    #pragma omp parallel shared(g, m, input_size, n_modules, S_index, nz_S, t, half_input_size, shift_dim, \
                                    cantor_combinations)
        {
            for (int r = static_cast<int>(m) - 1; r >= 0; --r) {
    #pragma omp single
                {
                    nz_S            = s_i[S_index];
                    t               = n_terms[S_index];
                    n_modules >>= 1;
                    input_size <<= 1;
                    half_input_size = input_size >> 1;
                }
    #pragma omp barrier

                if (n_modules >= 2) {
    #pragma omp for schedule(static)
                    for (size_t module = 0; module < n_modules; ++module) {
                        cantor_additive_ifft_combination_radix2_module<FieldT>(
                            g, module * input_size, half_input_size, module, m, shift_dim, n_modules, nz_S, t,
                            cantor_combinations);
                    }
                } else {
    #pragma omp single
                    {
                        cantor_additive_ifft_combination_radix2_module<FieldT>(
                            g, 0, half_input_size, 0, m, shift_dim, n_modules, nz_S, t, cantor_combinations);
                    }
                }

    #pragma omp single
                {
                    ++S_index;
                }
    #pragma omp barrier
            }
        }

        return g;
    #endif
    }

    template<typename FieldT>
    std::vector<FieldT> additive_IFFT_parallel(const std::vector<FieldT> &evals,
                                               const libiop::affine_subspace<FieldT> &domain)
    {
    #ifndef _OPENMP
        return additive_IFFT<FieldT>(evals, domain);
    #else
        std::vector<FieldT> g(evals);
        g.resize(domain.num_elements(), FieldT::zero());

        const FieldT              affine_shift = domain.shift();
        const std::vector<FieldT> W(domain.basis());

        const size_t n = g.size();
        const size_t m = domain.dimension();
        assert(n == (1ull << m));

        size_t              input_size         = 1;
        size_t              n_modules          = n;
        size_t              S_index            = 0;
        FieldT              affine_shift_round = FieldT::zero();
        size_t              half_input_size    = 0;
        std::vector<size_t> nz_S;

    #pragma omp parallel shared(g, affine_shift, W, m, input_size, n_modules, S_index, affine_shift_round, \
                                    half_input_size, nz_S)
        {
            for (int r = static_cast<int>(m) - 1; r >= 0; --r) {
    #pragma omp single
                {
                    affine_shift_round = affine_shift;
                    nz_S.clear();
                    if (S_index > 0) {
                        nz_S.reserve(S_index);
                        nz_S.emplace_back((1ull << S_index) - 1);
                        for (size_t i = 1; i < S_index; ++i) {
                            size_t ii = i;
                            size_t rr = S_index;
                            while (((rr & 1) | (~ii & 1)) && ii > 0) {
                                ii >>= 1;
                                rr >>= 1;
                            }
                            if (ii == 0) {
                                nz_S.emplace_back((1ull << S_index) - (1ull << i));
                                affine_shift_round += affine_shift ^ (1ull << i);
                            }
                        }
                        affine_shift_round += affine_shift ^ (1ull << S_index);
                    }
                    n_modules >>= 1;
                    input_size <<= 1;
                    half_input_size = input_size >> 1;
                }
    #pragma omp barrier

                if (n_modules >= 2) {
    #pragma omp for schedule(static)
                    for (size_t module = 0; module < n_modules; ++module) {
                        cantor_additive_ifft_radix2_process_module(
                            g, module * input_size, half_input_size, module, static_cast<size_t>(r),
                            affine_shift_round, W, nz_S);
                    }
                } else {
    #pragma omp single
                    {
                        cantor_additive_ifft_radix2_process_module(g, 0, half_input_size, 0, static_cast<size_t>(r),
                                                                   affine_shift_round, W, nz_S);
                    }
                }

    #pragma omp single
                {
                    ++S_index;
                }
    #pragma omp barrier
            }
        }

        return g;
    #endif
    }

    template<size_t K, typename FieldT>
    std::vector<FieldT> additive_IFFT_radix2k_parallel(const std::vector<FieldT> &evals,
                                                       const size_t domain_dim,
                                                       const size_t shift_dim)
    {
    #ifndef _OPENMP
        return additive_IFFT_radix2k<K, FieldT>(evals, domain_dim, shift_dim);
    #else
        static_assert(K >= 1, "K must be >= 1");
        static_assert(K <= 5, "K supported up to 5 (radix up to 32)");

        const size_t m = domain_dim;
        std::vector<FieldT> g(evals);
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

        size_t input_size = 1;
        size_t n_modules  = n;

        std::array<const unsigned *, K + 1> nz_ptr{};
        std::array<unsigned, K + 1>        t_arr{};
        size_t                             L     = 0;
        size_t                             logL  = 0;
        size_t                             chunk = 0;
        size_t                             half  = 0;

        size_t         round_idx = 0;
        const unsigned *nz_res   = nullptr;
        unsigned       t_res     = 0;
        size_t         L_res     = 0;
        size_t         half_res  = 0;

    #pragma omp parallel shared(g, m, shift_dim, cantor_combinations, num_rounds, residual, input_size, n_modules, \
                                    nz_ptr, t_arr, L, logL, chunk, half, round_idx, nz_res, t_res, L_res, half_res)
        {
            for (size_t r_res = residual; r_res-- > 0;) {
    #pragma omp single
                {
                    round_idx = K * num_rounds + r_res;
                    n_modules >>= 1;
                    input_size <<= 1;
                    const size_t S_index = m - round_idx - 1;
                    nz_res               = s_i[S_index];
                    t_res                = n_terms[S_index];
                    L_res                = input_size;
                    half_res             = L_res >> 1;
                }
    #pragma omp barrier

                if (n_modules >= 2) {
    #pragma omp for schedule(static)
                    for (size_t module = 0; module < n_modules; ++module) {
                        cantor_additive_ifft_combination_radix2_module<FieldT>(
                            g, module * L_res, half_res, module, m, shift_dim, n_modules, nz_res, t_res,
                            cantor_combinations);
                    }
                } else {
    #pragma omp single
                    {
                        cantor_additive_ifft_combination_radix2_module<FieldT>(
                            g, 0, half_res, 0, m, shift_dim, n_modules, nz_res, t_res, cantor_combinations);
                    }
                }
    #pragma omp barrier
            }

            for (size_t round_j = num_rounds; round_j-- > 0;) {
    #pragma omp single
                {
                    n_modules >>= K;
                    input_size <<= K;
                    L     = input_size;
                    logL  = m - K * round_j;
                    chunk = L >> K;
                    half  = L >> 1;
                    for (size_t ell = 1; ell <= K; ++ell) {
                        const size_t S_idx = logL - ell;
                        nz_ptr[ell]        = s_i[S_idx];
                        t_arr[ell]         = n_terms[S_idx];
                    }
                }
    #pragma omp barrier

                if (n_modules >= 2) {
    #pragma omp for schedule(static)
                    for (size_t module = 0; module < n_modules; ++module) {
                        cantor_additive_ifft_combination_radix2k_process_module<K, FieldT>(
                            g, module * L, L, chunk, half, logL, module, m, shift_dim, n_modules,
                            cantor_combinations, nz_ptr, t_arr);
                    }
                } else {
    #pragma omp single
                    {
                        cantor_additive_ifft_combination_radix2k_process_module<K, FieldT>(
                            g, 0, L, chunk, half, logL, 0, m, shift_dim, n_modules, cantor_combinations, nz_ptr,
                            t_arr);
                    }
                }
    #pragma omp barrier
            }
        }

        return g;
    #endif
    }

    template<size_t K, typename FieldT>
    std::vector<FieldT> additive_IFFT_radix2k_parallel(const std::vector<FieldT> &evals,
                                                       const libiop::affine_subspace<FieldT> &domain)
    {
    #ifndef _OPENMP
        return additive_IFFT_radix2k<K, FieldT>(evals, domain);
    #else
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

        std::array<std::vector<size_t>, K + 1> nz_S_main;
        std::array<FieldT, K + 1>              asr_main;
        size_t                                 L        = 0;
        size_t                                 logL     = 0;
        size_t                                 chunk    = 0;
        size_t                                 half     = 0;
        size_t                                 num_bits = 0;

        std::vector<size_t> nz_S_res;
        FieldT              asr_res   = FieldT::zero();
        size_t              L_res     = 0;
        size_t              half_res  = 0;
        size_t              round_idx = 0;

    #pragma omp parallel shared(g, affine_shift, W, m, num_rounds, residual, input_size, n_modules, block_offsets, \
                                    nz_S_main, asr_main, L, logL, chunk, half, num_bits, nz_S_res, asr_res, L_res,      \
                                    half_res, round_idx)
        {
            for (size_t r_res = residual; r_res-- > 0;) {
    #pragma omp single
                {
                    round_idx = K * num_rounds + r_res;
                    n_modules >>= 1;
                    input_size <<= 1;
                    const size_t S_idx = m - round_idx - 1;
                    L_res              = input_size;
                    half_res           = L_res >> 1;
                    nz_S_res.clear();
                    asr_res = affine_shift;
                    if (S_idx > 0) {
                        nz_S_res.reserve(S_idx);
                        nz_S_res.emplace_back((1ull << S_idx) - 1);
                        for (size_t i = 1; i < S_idx; ++i) {
                            size_t ii = i, rr = S_idx;
                            while (((rr & 1) | (~ii & 1)) && ii > 0) {
                                ii >>= 1;
                                rr >>= 1;
                            }
                            if (ii == 0) {
                                nz_S_res.emplace_back((1ull << S_idx) - (1ull << i));
                                asr_res += affine_shift ^ (1ull << i);
                            }
                        }
                        asr_res += affine_shift ^ (1ull << S_idx);
                    }
                }
    #pragma omp barrier

                if (n_modules >= 2) {
    #pragma omp for schedule(static)
                    for (size_t module = 0; module < n_modules; ++module) {
                        cantor_additive_ifft_radix2_process_module(g, module * L_res, half_res, module, round_idx,
                                                                   asr_res, W, nz_S_res);
                    }
                } else {
    #pragma omp single
                    {
                        cantor_additive_ifft_radix2_process_module(g, 0, half_res, 0, round_idx, asr_res, W, nz_S_res);
                    }
                }
    #pragma omp barrier
            }

            for (size_t round_j = num_rounds; round_j-- > 0;) {
    #pragma omp single
                {
                    n_modules >>= K;
                    input_size <<= K;
                    L        = input_size;
                    logL     = m - K * round_j;
                    chunk    = L >> K;
                    half     = L >> 1;
                    num_bits = K * round_j;
                    for (size_t ell = 1; ell <= K; ++ell) {
                        const size_t S_idx = logL - ell;
                        asr_main[ell]      = affine_shift;
                        nz_S_main[ell].clear();
                        if (S_idx > 0) {
                            nz_S_main[ell].reserve(S_idx);
                            nz_S_main[ell].emplace_back((1ull << S_idx) - 1);
                            for (size_t i = 1; i < S_idx; ++i) {
                                size_t ii = i, rr = S_idx;
                                while (((rr & 1) | (~ii & 1)) && ii > 0) {
                                    ii >>= 1;
                                    rr >>= 1;
                                }
                                if (ii == 0) {
                                    nz_S_main[ell].emplace_back((1ull << S_idx) - (1ull << i));
                                    asr_main[ell] += affine_shift ^ (1ull << i);
                                }
                            }
                            asr_main[ell] += affine_shift ^ (1ull << S_idx);
                        }
                    }
                }
    #pragma omp barrier

                if (n_modules >= 2) {
    #pragma omp for schedule(static)
                    for (size_t module = 0; module < n_modules; ++module) {
                        cantor_additive_ifft_radix2k_process_module<K>(
                            g, module * L, L, chunk, half, logL, num_bits, module, nz_S_main, asr_main, W,
                            block_offsets);
                    }
                } else {
    #pragma omp single
                    {
                        cantor_additive_ifft_radix2k_process_module<K>(
                            g, 0, L, chunk, half, logL, num_bits, 0, nz_S_main, asr_main, W, block_offsets);
                    }
                }
    #pragma omp barrier
            }
        }

        return g;
    #endif
    }

} // namespace cantor