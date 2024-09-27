#include <Cantor/fft.hpp>

#include <iostream>
#include "utils/utils.hpp"


template<typename FieldT>
std::vector<FieldT> cantor_additive_FFT(const std::vector<FieldT> &poly_coeffs,
                                        const libiop::affine_subspace<FieldT> &domain)
{
    std::vector<FieldT> g(poly_coeffs);
    S.resize(domain.num_elements(), FieldT::zero());
    
    FieldT affine_shift = domain.shift();
    const size_t n = S.size();
    const size_t m = domain.dimension();
    assert(n == (1ull<<m));


    size_t input_size = n;
    size_t n_modules = 1;
    size_t S_index = m;
    for (size_t r = 0; r < m; ++r)
    {
        // Computing non-zero indices in S_{m-r}(X) for round r in the reverse order. Becuse in division algorithm we start from 
        --S_index; 
        std::vector<size_t> nz_S;
        nz_S.reserve(S_index);
        nz_S.emplace_back((1<<S_index)-1); // C(x,0) = 1, so we assigned that before loop.

        for (size_t i = 1; i < S_index; ++i){
            bool is_odd = (S_index & 1) | (~i & 1);
            size_t ii = i >> 1;
            size_t rr = S_index >> 1;
            while (is_odd && ii>0){
                is_odd = (rr & 1) | (~ii & 1);
                ii >>= 1;
                rr >>= 1;
            }
            if (is_odd){
                nz_S.emplace_back((1<<S_index) - (1<<i));
            }
        }
        nz_S.emplace_back(0); // C(x,x) = 1

        
    }



    return S;
}