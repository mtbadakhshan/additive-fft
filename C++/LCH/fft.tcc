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

    // Radix-4 variant of additive_FFT. Same basis_conversion + replication
    // setup as the radix-2 path, but the post-conversion butterfly fuses two
    // consecutive radix-2 stages per round, halving the number of cache
    // sweeps over the coefficient buffer (with a single residual radix-2
    // stage when log2(poly terms) is odd).
    template<typename FieldT>
    std::vector<FieldT> additive_FFT_radix4(const std::vector<FieldT> &poly_coeffs,
                                            const size_t domain_dim, const size_t shift_dim){
        const size_t m = domain_dim;
        std::vector<FieldT> g(poly_coeffs);
        size_t log_poly_terms = libff::log2(poly_coeffs.size());
        g.resize((1ULL << m), FieldT::zero());
        const size_t n = g.size();
        assert(n == (1ull<<m));

        FieldT* cantor_combinations;
        if(FieldT::extension_degree() == 128) cantor_combinations = (FieldT*) cantor::cantor_combinations_8R_in_gf2to128;
        else if(FieldT::extension_degree() == 192) cantor_combinations = (FieldT*) cantor::cantor_combinations_8R_in_gf2to192;
        else if(FieldT::extension_degree() == 256) cantor_combinations = (FieldT*) cantor::cantor_combinations_8R_in_gf2to256;
        else throw std::invalid_argument("The field size should be either 128, 192, or 256 for using the cantor basis");

        basis_conversion(g, 1<<log_poly_terms);

        unsigned index = 1<<log_poly_terms;
        while(index < n){
            std::copy(g.begin(), g.begin()+poly_coeffs.size(), g.begin() + index);
            index += (1<<log_poly_terms);
        }

        butterfly_radix4(g, 1<<log_poly_terms, shift_dim, cantor_combinations);
        return g;
    }

    // Templated radix-2^K LCH AFFT. K=1 reduces to the radix-2 baseline,
    // K=2 to additive_FFT_radix4, and K in {3, 4, 5} fuses 3, 4, 5 consecutive
    // radix-2 stages into a single in-register K-stage butterfly per outer
    // round. Same basis_conversion + replication setup as the other variants.
    template<size_t K, typename FieldT>
    std::vector<FieldT> additive_FFT_radix2k(const std::vector<FieldT> &poly_coeffs,
                                             const size_t domain_dim, const size_t shift_dim){
        const size_t m = domain_dim;
        std::vector<FieldT> g(poly_coeffs);
        size_t log_poly_terms = libff::log2(poly_coeffs.size());
        g.resize((1ULL << m), FieldT::zero());
        const size_t n = g.size();
        assert(n == (1ull<<m));

        FieldT* cantor_combinations;
        if(FieldT::extension_degree() == 128) cantor_combinations = (FieldT*) cantor::cantor_combinations_8R_in_gf2to128;
        else if(FieldT::extension_degree() == 192) cantor_combinations = (FieldT*) cantor::cantor_combinations_8R_in_gf2to192;
        else if(FieldT::extension_degree() == 256) cantor_combinations = (FieldT*) cantor::cantor_combinations_8R_in_gf2to256;
        else throw std::invalid_argument("The field size should be either 128, 192, or 256 for using the cantor basis");

        basis_conversion(g, 1<<log_poly_terms);

        unsigned index = 1<<log_poly_terms;
        while(index < n){
            std::copy(g.begin(), g.begin()+poly_coeffs.size(), g.begin() + index);
            index += (1<<log_poly_terms);
        }

        butterfly_radix2k<K, FieldT>(g, 1<<log_poly_terms, shift_dim, cantor_combinations);
        return g;
    }

    template<typename FieldT>
    std::vector<FieldT> additive_FFT_parallel(const std::vector<FieldT> &poly_coeffs,
                                               const size_t domain_dim, const size_t shift_dim)
    {
        const size_t m = domain_dim;
        std::vector<FieldT> g(poly_coeffs);
        size_t log_poly_terms = libff::log2(poly_coeffs.size());
        g.resize((1ULL << m), FieldT::zero());
        const size_t n = g.size();
        assert(n == (1ull << m));

        FieldT *cantor_combinations;
        if (FieldT::extension_degree() == 128)
            cantor_combinations = (FieldT *)cantor::cantor_combinations_8R_in_gf2to128;
        else if (FieldT::extension_degree() == 192)
            cantor_combinations = (FieldT *)cantor::cantor_combinations_8R_in_gf2to192;
        else if (FieldT::extension_degree() == 256)
            cantor_combinations = (FieldT *)cantor::cantor_combinations_8R_in_gf2to256;
        else
            throw std::invalid_argument(
                "The field size should be either 128, or 256 for using the cantor basis");

        basis_conversion(g, 1 << log_poly_terms);

        unsigned index = 1 << log_poly_terms;
        while (index < n) {
            std::copy(g.begin(), g.begin() + poly_coeffs.size(), g.begin() + index);
            index += (1 << log_poly_terms);
        }

        butterfly_parallel(g, 1 << log_poly_terms, shift_dim, cantor_combinations);
        return g;
    }

    template<typename FieldT>
    std::vector<FieldT> additive_FFT_radix4_parallel(const std::vector<FieldT> &poly_coeffs,
                                                     const size_t domain_dim, const size_t shift_dim)
    {
        return additive_FFT_radix2k_parallel<2, FieldT>(poly_coeffs, domain_dim, shift_dim);
    }

    template<size_t K, typename FieldT>
    std::vector<FieldT> additive_FFT_radix2k_parallel(const std::vector<FieldT> &poly_coeffs,
                                                        const size_t domain_dim, const size_t shift_dim)
    {
        const size_t m = domain_dim;
        std::vector<FieldT> g(poly_coeffs);
        size_t log_poly_terms = libff::log2(poly_coeffs.size());
        g.resize((1ULL << m), FieldT::zero());
        const size_t n = g.size();
        assert(n == (1ull << m));

        FieldT *cantor_combinations;
        if (FieldT::extension_degree() == 128)
            cantor_combinations = (FieldT *)cantor::cantor_combinations_8R_in_gf2to128;
        else if (FieldT::extension_degree() == 192)
            cantor_combinations = (FieldT *)cantor::cantor_combinations_8R_in_gf2to192;
        else if (FieldT::extension_degree() == 256)
            cantor_combinations = (FieldT *)cantor::cantor_combinations_8R_in_gf2to256;
        else
            throw std::invalid_argument(
                "The field size should be either 128, 192, or 256 for using the cantor basis");

        basis_conversion(g, 1 << log_poly_terms);

        unsigned index = 1 << log_poly_terms;
        while (index < n) {
            std::copy(g.begin(), g.begin() + poly_coeffs.size(), g.begin() + index);
            index += (1 << log_poly_terms);
        }

        butterfly_radix2k_parallel<K, FieldT>(g, 1 << log_poly_terms, shift_dim, cantor_combinations);
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