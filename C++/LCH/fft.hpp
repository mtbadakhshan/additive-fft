#ifndef ADDITIVE_FFT_LCH_HPP_
#define ADDITIVE_FFT_LCH_HPP_

#include <cstddef>
#include <vector>

#include "libiop/algebra/field_subset.hpp"
#include "libiop/algebra/subspace.hpp"

namespace lch {

    template<typename FieldT>   
    std::vector<FieldT> additive_FFT(const std::vector<FieldT> &poly_coeffs, 
                                    const size_t domain_dim, const size_t shift_dim);
    
    template<typename FieldT>   
    std::vector<FieldT> additive_IFFT(const std::vector<FieldT> &poly_coeffs, 
                                    const size_t domain_dim, const size_t shift_dim);
}

#include "LCH/utils.tcc"
#include "LCH/fft.tcc"
#endif // ADDITIVE_FFT_LCH_HPP_
