
namespace lch {

    static inline
    unsigned deg_si( unsigned si ) {
	return (1<<si);
    }

    static inline 
    unsigned get_num_blocks( unsigned poly_len , unsigned blk_size ) {
        return poly_len/blk_size;
    }
    

    static inline
    unsigned get_si_2_pow( unsigned si , unsigned deg ) {
        unsigned si_deg = (1<<si);
        unsigned r=1;
        while( (si_deg<<r) < deg ) {
            r += 1;
        }
        return (1<<(r-1));
    }

    static inline
    unsigned get_max_si( unsigned deg ) {
        unsigned si = 0;
        unsigned si_attempt = 1;
        uint64_t deg64 = deg;
        while( deg64 > ((1ULL)<<si_attempt) ) {
            si = si_attempt;
            si_attempt <<= 1;
        }
        return si;
    }

    // DIRECT
    template<typename FieldT>   
    static inline
    void xor_down_128(  std::vector<FieldT> &poly , unsigned offset, unsigned st , unsigned len , unsigned diff )
    {
        for( unsigned i=0;i<len;i++) {
            poly[st-i-1 + offset] += poly[st-i-1+diff+offset]; // + stands for ^ in libff
        }
    }


    template<typename FieldT>   
    static inline
    void poly_div_128(  std::vector<FieldT> &poly, unsigned offset, unsigned n_terms , unsigned blk_size , unsigned si , unsigned pow )
    {
        if( 0 == si ) return;
        unsigned si_degree = deg_si(si)*pow;
        unsigned deg_diff = si_degree - pow;
        unsigned deg_blk = get_num_blocks( n_terms , blk_size ) -1;

        xor_down_128( poly , offset, (deg_blk-deg_diff+1)*blk_size , (deg_blk-si_degree+1)*blk_size , deg_diff*blk_size );
    }

    template<typename FieldT>   
    static inline
    void represent_in_si_128(  std::vector<FieldT> &poly , size_t offset_in, unsigned n_terms , unsigned blk_size , unsigned si )
    {
        if( 0 == si ) return;
        unsigned num_blocks = get_num_blocks( n_terms , blk_size );
        if( 2 >= num_blocks ) return;
        unsigned degree_in_blocks = num_blocks - 1;
        unsigned degree_basic_form_si = deg_si(si);
        if( degree_basic_form_si > degree_in_blocks ) return;

        unsigned pow = get_si_2_pow( si , degree_in_blocks );

        while( 0 < pow ) {
            for(unsigned offset=0; offset<n_terms; offset+= blk_size*2*pow*deg_si(si) ) {
                poly_div_128( poly, offset+offset_in , blk_size*2*pow*deg_si(si) , blk_size , si , pow );
            }
            pow >>= 1;
        }
    }


    template<typename FieldT>   
    void basis_conversion_recursive( std::vector<FieldT> &poly_coeffs, 
                                            size_t offset_in,
                                            const size_t n_terms, 
                                            const size_t blk_size){
        unsigned num_blocks =  get_num_blocks( n_terms , blk_size );
        if( 2 >= num_blocks ) return;
        unsigned degree_in_blocks = num_blocks - 1;
        unsigned si = get_max_si( degree_in_blocks );
        represent_in_si_128( poly_coeffs , offset_in, n_terms , blk_size , si );
        unsigned new_blk_size = deg_si(si)*blk_size;
        basis_conversion_recursive( poly_coeffs, offset_in, n_terms , new_blk_size );
        for(unsigned offset=0; offset<n_terms; offset+= new_blk_size ) {
            basis_conversion_recursive( poly_coeffs, offset+offset_in , new_blk_size , blk_size );
        }
    }

    template<typename FieldT>   
    void basis_conversion( std::vector<FieldT> &poly_coeffs, const size_t n_terms){
        basis_conversion_recursive(poly_coeffs, 0, n_terms, 1);
        
    }

    static inline
    unsigned get_s_k_a_cantor( unsigned k , unsigned a ) { return (a>>k); }


    template<typename FieldT>   
    static
    void butterfly_0( std::vector<FieldT> &poly , unsigned unit )
    {
        unsigned unit_2= unit/2;
        for(unsigned i=0;i<unit_2;i++) {
            poly[unit_2+i] += poly[i];
        }
    }

    template<typename FieldT>   
    static
    void butterfly_op( std::vector<FieldT> &poly, unsigned offset, unsigned unit , unsigned ska, size_t shift_bit, FieldT* cantor_combinations)
    {
        size_t module_shifted = (ska) | (shift_bit<<1), i_256 = 0;
        FieldT mult_factor = FieldT::zero();
        while(module_shifted){
            mult_factor += cantor_combinations[i_256 + (module_shifted & 0xff)];
            i_256 += 256;
            module_shifted >>= 8;
        }

        // std::cout<<"mult_factor: "<< mult_factor << " offset: "<<offset << " unit: " << unit << std::endl;

        unsigned unit_2= unit/2;
        for(unsigned i=0;i<unit_2;i++) {
            poly[offset + i] += poly[offset + unit_2+i] * mult_factor;
            poly[offset + unit_2+i] += poly[offset + i];
        }
    }


    template<typename FieldT>   
    void butterfly( std::vector<FieldT> &poly_coeffs,
                    const size_t n_terms, const size_t shift_dim, 
                    FieldT* cantor_combinations){
        if( 1 >= n_terms ) return;

	    unsigned log_n = __builtin_ctz( n_terms );
	    unsigned m     = __builtin_ctz( poly_coeffs.size());

        // std::cout<<"m: " << m << std::endl;
        // std::cout<<"log_n: " << log_n << std::endl ;

	    for(unsigned i=log_n; i > 0; i--) {
		    unsigned unit = (1<<i);
		    unsigned num = poly_coeffs.size() / unit;
            // std::cout<<"i: " << i <<" num: " << num << " unit: " << unit << std::endl ;
            size_t shift_bit = shift_dim==0 ? 0 :(num) << (shift_dim - m);

		    // butterfly_0( poly_coeffs , unit );
		for(unsigned j=0;j<num;j++) {
			butterfly_op( poly_coeffs, j*unit , unit , get_s_k_a_cantor( i-1 , j*unit ), shift_bit, cantor_combinations);
            // std::cout<<"get_s_k_a_cantor( i-1 , j*unit ): " << get_s_k_a_cantor( i-1 , j*unit ) << std::endl ;
            }
            // my_print_vector(poly_coeffs);
        }
    }

// RADIX-4 -------------------------------------------------------------------------------------------------------

    // Decode a bit-pattern into the field element whose Cantor-basis
    // coordinates equal those bits. The 8-bits-at-a-time table layout is
    //   cantor_combinations[256 * p + b] = element with Cantor coords b << (8p),
    // hence a 32-bit input is decomposed into four bytes and summed.
    template<typename FieldT>
    static inline
    FieldT element_from_cantor_bits(size_t bits, FieldT* cantor_combinations) {
        FieldT acc = FieldT::zero();
        size_t i_256 = 0;
        while (bits) {
            acc += cantor_combinations[i_256 + (bits & 0xff)];
            i_256 += 256;
            bits >>= 8;
        }
        return acc;
    }

    // Radix-4 LCH butterfly. Mathematically equivalent to two consecutive
    // radix-2 stages of `butterfly`, fused into a single ascending sweep over
    // the four quarters of each module so that each radix-4 round touches the
    // module exactly once. When log_n is odd we run (log_n-1)/2 radix-4 rounds
    // followed by a single residual radix-2 stage at i = 1.
    //
    // Inside one radix-4 round the four quarters undergo
    //
    //     | h_0 |   | 1   T_1   T_2   T_1 T_2       | | f_0 |
    //     | h_1 | = | 1  T_1+1  T_2  (T_1+1) T_2    | | f_1 |
    //     | h_2 |   | 1  T_1+b  T_2+1 (T_1+b)(T_2+1)| | f_2 |
    //     | h_3 |   | 1 T_1+b+1 T_2+1 (T_1+b+1)(T_2+1)| | f_3 |
    //
    // with b = beta_1 (so b^2 + b = 1, hence T_3 = T_1 + b is the twiddle on
    // the upper half during the inner stage). T_2 is the radix-2 twiddle of
    // the outer stage on this module; T_1 is the radix-2 twiddle of the inner
    // stage on the lower-half sub-module. The bit-pattern that encodes T_1
    // is the bit-pattern that encodes T_2 shifted up by one (because
    // S^{i-2}(beta_k) = beta_{k-i+2} = beta_{(k-i+1)+1}).
    template<typename FieldT>
    void butterfly_radix4(std::vector<FieldT> &poly_coeffs,
                          const size_t n_terms, const size_t shift_dim,
                          FieldT* cantor_combinations) {
        if (1 >= n_terms) return;

        const unsigned log_n = __builtin_ctz(n_terms);
        const unsigned m     = __builtin_ctz(poly_coeffs.size());

        // beta_1 lives at cantor_combinations[2] (bits-of-2 = bit 1 set, so
        // the encoded element has Cantor coordinate beta_1 only).
        const FieldT beta_1 = cantor_combinations[2];

        const unsigned num_radix4 = log_n / 2;

        for (unsigned outer = 0; outer < num_radix4; ++outer) {
            const unsigned i_high  = log_n - 2 * outer;
            const unsigned unit    = 1u << i_high;
            const unsigned half    = unit >> 1;
            const unsigned quarter = unit >> 2;
            const size_t   num_modules = poly_coeffs.size() / unit;

            const size_t shift_bit_high = (shift_dim == 0) ? 0
                : (num_modules << (shift_dim - m));

            for (size_t mod = 0; mod < num_modules; ++mod) {
                const size_t offset = mod * unit;

                const size_t mshift_T2 = (mod << 1) | (shift_bit_high << 1);
                const size_t mshift_T1 = mshift_T2 << 1;

                const FieldT T_2 = element_from_cantor_bits<FieldT>(mshift_T2, cantor_combinations);
                const FieldT T_1 = element_from_cantor_bits<FieldT>(mshift_T1, cantor_combinations);
                const FieldT T_3 = T_1 + beta_1;

                // Single ascending sweep over the four quarters; 4 reads and
                // 4 writes per inner iteration, all intermediates kept in
                // registers.
                for (unsigned jj = 0; jj < quarter; ++jj) {
                    const FieldT q0 = poly_coeffs[offset + jj];
                    const FieldT q1 = poly_coeffs[offset + quarter + jj];
                    const FieldT q2 = poly_coeffs[offset + half + jj];
                    const FieldT q3 = poly_coeffs[offset + half + quarter + jj];

                    const FieldT u0 = q0 + T_2 * q2;
                    const FieldT u1 = q1 + T_2 * q3;
                    const FieldT u2 = u0 + q2;
                    const FieldT u3 = u1 + q3;

                    const FieldT h0 = u0 + T_1 * u1;
                    const FieldT h2 = u2 + T_3 * u3;
                    const FieldT h1 = h0 + u1;
                    const FieldT h3 = h2 + u3;

                    poly_coeffs[offset + jj]                  = h0;
                    poly_coeffs[offset + quarter + jj]        = h1;
                    poly_coeffs[offset + half + jj]           = h2;
                    poly_coeffs[offset + half + quarter + jj] = h3;
                }
            }
        }

        // Mixed-radix epilogue: when log_n is odd, do the final radix-2 stage
        // at i = 1 (unit = 2). This is exactly butterfly_op with unit = 2.
        if (log_n & 1u) {
            const unsigned unit = 2;
            const size_t   num_modules = poly_coeffs.size() / unit;
            const size_t   shift_bit = (shift_dim == 0) ? 0
                : (num_modules << (shift_dim - m));

            for (size_t mod = 0; mod < num_modules; ++mod) {
                butterfly_op(poly_coeffs, mod * unit, unit,
                             get_s_k_a_cantor(0, mod * unit),
                             shift_bit, cantor_combinations);
            }
        }
    }

    // Templated radix-2^K LCH butterfly. Generalises butterfly_radix4
    // (which is the K = 2 specialisation). Per outer round, processes a module
    // of size L = 2^{logL} as Q = 2^K quarters of size L/Q in a single
    // ascending sweep, executing K consecutive radix-2 stages in registers.
    //
    // Twiddle structure: for round_j and module index `mod`, level `ell` in
    // the in-register K-stage butterfly uses the radix-2 stage twiddle for
    // i_ell = logL - ell + 1. The bit-pattern for the radix-2 twiddle is
    //
    //     mshift_ell = (mod << ell) | (shift_bit_ell << 1),
    //
    // and shift_bit_ell = shift_bit_1 << (ell - 1), so we just shift mshift_1
    // up by (ell-1) bits to get mshift_ell. Adding `block_offsets[b]` (the
    // element with cantor coords (b << 1)) yields the per-block twiddle for
    // the b-th radix-2 sub-module inside the radix-2^K module at level ell.
    //
    // When log_n is not a multiple of K, the function performs (log_n / K)
    // radix-2^K rounds and then (log_n mod K) residual radix-2 rounds.
    template<size_t K, typename FieldT>
    void butterfly_radix2k(std::vector<FieldT> &poly_coeffs,
                           const size_t n_terms, const size_t shift_dim,
                           FieldT* cantor_combinations) {
        static_assert(K >= 1, "K must be >= 1");
        static_assert(K <= 5, "K supported up to 5 (radix up to 32)");

        if (1 >= n_terms) return;

        const unsigned log_n = __builtin_ctz(n_terms);
        const unsigned m     = __builtin_ctz(poly_coeffs.size());

        constexpr size_t Q      = 1ull << K;
        constexpr size_t Q_half = Q >> 1;

        // block_offsets[c] = element with Cantor coords (c << 1) for
        // c in [0, 2^{K-1}). Equivalently, block_offsets[c] = sum of
        // beta_{i+1} over set bits i of c. cantor_combinations[1<<k]
        // is the Cantor basis element beta_k.
        std::array<FieldT, (Q_half == 0 ? 1 : Q_half)> block_offsets;
        block_offsets[0] = FieldT::zero();
        for (size_t c = 1; c < Q_half; ++c) {
            const size_t low_bit = __builtin_ctzll(c);
            block_offsets[c] = block_offsets[c ^ (1ull << low_bit)]
                             + cantor_combinations[1ull << (low_bit + 1)];
        }

        const unsigned num_rounds = log_n / K;
        const unsigned residual   = log_n - K * num_rounds;

        // ---- Main: num_rounds radix-2^K rounds ----
        for (unsigned round_j = 0; round_j < num_rounds; ++round_j) {
            const unsigned logL = log_n - K * round_j;
            const size_t   L     = 1ull << logL;
            const size_t   chunk = L >> K;
            const size_t   num_modules = poly_coeffs.size() >> logL;

            // shift_bit for stage i_1 = logL: matches radix-2 butterfly's
            // (poly_coeffs.size() / 2^logL) << (shift_dim - m).
            const size_t shift_bit_1 = (shift_dim == 0) ? 0
                : (num_modules << (shift_dim - m));

            for (size_t mod = 0; mod < num_modules; ++mod) {
                const size_t offset = mod * L;

                // mshift_ell = mshift_1 << (ell - 1).
                const size_t mshift_1 = (mod << 1) | (shift_bit_1 << 1);
                std::array<FieldT, K + 1> base_T;
                for (size_t ell = 1; ell <= K; ++ell) {
                    base_T[ell] = element_from_cantor_bits<FieldT>(
                        mshift_1 << (ell - 1), cantor_combinations);
                }

                // Single ascending sweep: K-stage in-register butterfly per
                // tuple. 2^K reads, 2^K writes, K * 2^{K-1} mults+adds in
                // registers.
                for (size_t jj = 0; jj < chunk; ++jj) {
                    FieldT v[Q];
                    for (size_t i = 0; i < Q; ++i)
                        v[i] = poly_coeffs[offset + i * chunk + jj];

                    for (size_t ell = 1; ell <= K; ++ell) {
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
                        poly_coeffs[offset + i * chunk + jj] = v[i];
                }
            }
        }

        // ---- Mixed-radix epilogue: (log_n mod K) residual radix-2 rounds ----
        for (unsigned r_res = 0; r_res < residual; ++r_res) {
            const unsigned i_stage = log_n - K * num_rounds - r_res;
            const unsigned unit    = 1u << i_stage;
            const size_t   num     = poly_coeffs.size() >> i_stage;
            const size_t   shift_bit = (shift_dim == 0) ? 0
                : (num << (shift_dim - m));

            for (size_t mod = 0; mod < num; ++mod) {
                butterfly_op(poly_coeffs, mod * unit, unit,
                             get_s_k_a_cantor(i_stage - 1, mod * unit),
                             shift_bit, cantor_combinations);
            }
        }
    }

// INVERSE ---------------------------------------------------------------------------------------------------------

    template<typename FieldT>   
    static inline
    void xor_up_128( std::vector<FieldT> &poly , unsigned offset, unsigned st , unsigned len , unsigned diff )
    {
        for( unsigned i=0;i<len;i++) {
            poly[st+i + offset] += poly[st+i+diff + offset];
        }
    }

    template<typename FieldT> 
    static inline
    void i_poly_div_128( std::vector<FieldT> &poly,  unsigned offset, unsigned n_terms , unsigned blk_size , unsigned si , unsigned pow )
    {
        if( 0 == si ) return;
        unsigned si_degree = deg_si(si)*pow;
        unsigned deg_diff = si_degree - pow;
        unsigned deg_blk = get_num_blocks( n_terms , blk_size ) -1;

        xor_up_128( poly , offset, blk_size*(si_degree-deg_diff) , (deg_blk-si_degree+1)*blk_size , deg_diff*blk_size );

    }

    template<typename FieldT> 
    static inline
    void i_represent_in_si_128( std::vector<FieldT> &poly, size_t offset_in, unsigned n_terms , unsigned blk_size , unsigned si )
    {
        if( 0 == si ) return;
        unsigned num_blocks = get_num_blocks( n_terms , blk_size );
        if( 2 >= num_blocks ) return;
        unsigned degree_in_blocks = num_blocks - 1;
        unsigned degree_basic_form_si = deg_si(si);
        if( degree_basic_form_si > degree_in_blocks ) return;

        unsigned pow = 1;
        while( pow*deg_si(si) <= degree_in_blocks ) {
            for(unsigned offset=0; offset<n_terms; offset+= blk_size*2*pow*deg_si(si) ) {
                i_poly_div_128( poly, offset+offset_in, blk_size*2*pow*deg_si(si) , blk_size , si , pow );
            }
            pow *= 2;
        }
    }

    template<typename FieldT> 
    void inv_basis_conversion_recursive( std::vector<FieldT> &poly_coeffs, 
                                         size_t offset_in,
                                         const size_t n_terms, 
                                         const size_t blk_size )
    {

        unsigned num_blocks = get_num_blocks( n_terms , blk_size );
        if( 2 >= num_blocks ) return;
        unsigned degree_in_blocks = num_blocks - 1;
        unsigned si = get_max_si( degree_in_blocks );
        unsigned new_blk_size = deg_si(si)*blk_size;
        for(unsigned offset=0; offset<n_terms; offset+= new_blk_size ) {
            inv_basis_conversion_recursive( poly_coeffs, offset+offset_in , new_blk_size , blk_size );
        }
        inv_basis_conversion_recursive( poly_coeffs, offset_in, n_terms , new_blk_size );
        i_represent_in_si_128( poly_coeffs, offset_in, n_terms , blk_size , si );
    }


    template<typename FieldT>   
    void inv_basis_conversion( std::vector<FieldT> &poly_coeffs, const size_t n_terms){
        inv_basis_conversion_recursive(poly_coeffs, 0, n_terms, 1);
    }


    template<typename FieldT>   
    static
    void i_butterfly_op( std::vector<FieldT> &poly, unsigned offset, unsigned unit , unsigned ska, size_t shift_bit, FieldT* cantor_combinations )
    {     
        size_t module_shifted = (ska) | (shift_bit<<1), i_256 = 0;
        FieldT mult_factor = FieldT::zero();
        while(module_shifted){
            mult_factor += cantor_combinations[i_256 + (module_shifted & 0xff)];
            i_256 += 256;
            module_shifted >>= 8;
        }

        unsigned unit_2= unit/2;
        for(unsigned i=0;i<unit_2;i++) {
            poly[offset + unit_2+i] += poly[offset + i];
            poly[offset + i] += poly[offset + unit_2+i] * mult_factor;
        }
    }

    template<typename FieldT>   
    void inv_butterfly( std::vector<FieldT> &poly_coeffs,
                        const size_t n_terms, const size_t shift_dim, 
                        FieldT* cantor_combinations){
        if( 1 >= n_terms ) return;

        unsigned log_n = __builtin_ctz( n_terms );

        for(unsigned i=1; i <= log_n; i++) {
            unsigned unit = (1<<i);
            unsigned num = n_terms / unit;
            size_t shift_bit = shift_dim==0 ? 0 :(1<<(log_n-i)) << (shift_dim - log_n);

            // butterfly_0( poly_coeffs , unit );
            for(unsigned j=0;j<num;j++) {
                i_butterfly_op( poly_coeffs, j*unit , unit , get_s_k_a_cantor( i-1 , j*unit ), shift_bit, cantor_combinations); 
            }
        }
    }
}
