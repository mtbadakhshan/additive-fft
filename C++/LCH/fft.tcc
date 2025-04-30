#include <LCH/fft.hpp>
#include <cstddef>
#include <libff/algebra/field_utils/algorithms.hpp>

#include <iostream>
#include "utils/utils.hpp"
#include "Cantor/cantor_basis.hpp"

namespace lch {
        
    template<typename FieldT>   
    std::vector<FieldT> additive_FFT(const std::vector<FieldT> &poly_coeffs, 
                                    const size_t domain_dim, const size_t shift_dim){
        const size_t m = domain_dim;
        std::vector<FieldT> g(poly_coeffs);
        size_t log_poly_terms = libff::log2(poly_coeffs.size());
        g.resize((1ULL << m), FieldT::zero());
        const size_t n = g.size();
        assert(n == (1ull<<m));

        // std::cout<<"log_poly_terms: "<<log_poly_terms<<std::endl;
        // std::cout<<"poly_coeffs.size(): "<<poly_coeffs.size()<<std::endl;
        // std::cout<<"n: "<<n<<std::endl;

        
        FieldT* cantor_combinations;
        if(FieldT::extension_degree() == 128) cantor_combinations = (FieldT*) cantor::cantor_combinations_8R_in_gf2to128;
        else if(FieldT::extension_degree() == 192) cantor_combinations = (FieldT*) cantor::cantor_combinations_8R_in_gf2to192;
        else if(FieldT::extension_degree() == 256) cantor_combinations = (FieldT*) cantor::cantor_combinations_8R_in_gf2to256;
        else throw std::invalid_argument("The field size should be either 128, or 256 for using the cantor basis");
    

        basis_conversion(g, 1<<log_poly_terms);
        // my_print_vector(g);
        unsigned index = 1<<log_poly_terms;
        while(index < n){
            std::copy(g.begin(), g.begin()+poly_coeffs.size(), g.begin() + index);
            index += (1<<log_poly_terms);
            // std::cout<<"copied\n";
        }
        // my_print_vector(g);


        butterfly(g, 1<<log_poly_terms, shift_dim, cantor_combinations);
        return g;
    }

    template<typename FieldT>   
    std::vector<FieldT> additive_IFFT(const std::vector<FieldT> &poly_coeffs, 
                                    const size_t domain_dim, const size_t shift_dim){
        const size_t m = domain_dim;
        std::vector<FieldT> g(poly_coeffs);
        g.resize((1ULL << m), FieldT::zero());
        const size_t n = g.size();
        assert(n == (1ull<<m));
        
        FieldT* cantor_combinations;
        if(FieldT::extension_degree() == 128) cantor_combinations = (FieldT*) cantor::cantor_combinations_8R_in_gf2to128;
        else if(FieldT::extension_degree() == 192) cantor_combinations = (FieldT*) cantor::cantor_combinations_8R_in_gf2to192;
        else if(FieldT::extension_degree() == 256) cantor_combinations = (FieldT*) cantor::cantor_combinations_8R_in_gf2to256;
        else throw std::invalid_argument("The field size should be either 128, or 256 for using the cantor basis");
    
        inv_butterfly(g, n, shift_dim, cantor_combinations);
        inv_basis_conversion(g, n);
        return g;
    }

} //namespace lch