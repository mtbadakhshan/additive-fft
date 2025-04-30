
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
