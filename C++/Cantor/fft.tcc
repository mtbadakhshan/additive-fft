#include <Cantor/fft.hpp>
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
