#ifndef ADDITIVE_FFT_LCH_HPP_
#define ADDITIVE_FFT_LCH_HPP_

#include <vector>


namespace lch {

    template<typename FieldT>   
    std::vector<FieldT> additive_FFT(const std::vector<FieldT> &poly_coeffs, 
                                    const size_t domain_dim, const size_t shift_dim);
    
    template<typename FieldT>   
    std::vector<FieldT> additive_IFFT(const std::vector<FieldT> &poly_coeffs, 
                                    const size_t domain_dim, const size_t shift_dim);
}

#include "utils.tcc"
#include "fft.tcc"
#endif // ADDITIVE_FFT_LCH_HPP_
