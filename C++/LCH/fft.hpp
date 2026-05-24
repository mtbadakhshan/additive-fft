#ifndef ADDITIVE_FFT_LCH_HPP_
#define ADDITIVE_FFT_LCH_HPP_

#include <array>
#include <cstddef>
#include <vector>

#include "libiop/algebra/field_subset.hpp"
#include "libiop/algebra/subspace.hpp"

namespace lch {

    template<typename FieldT>   
    std::vector<FieldT> additive_FFT(const std::vector<FieldT> &poly_coeffs, 
                                    const size_t domain_dim, const size_t shift_dim);

    template<size_t K, typename FieldT>
    std::vector<FieldT> additive_FFT_radix2k(const std::vector<FieldT> &poly_coeffs,
                                             const size_t domain_dim, const size_t shift_dim);

    /// OpenMP parallel over modules (radix-2 LCH butterfly). Without \c _OPENMP,
    /// forwards to \ref additive_FFT.
    template<typename FieldT>
    std::vector<FieldT> additive_FFT_parallel(const std::vector<FieldT> &poly_coeffs,
                                              const size_t domain_dim, const size_t shift_dim);

    /// OpenMP parallel \ref additive_FFT_radix2k. Without \c _OPENMP,
    /// forwards to the serial implementation.
    template<size_t K, typename FieldT>
    std::vector<FieldT> additive_FFT_radix2k_parallel(const std::vector<FieldT> &poly_coeffs,
                                                      const size_t domain_dim, const size_t shift_dim);

    template<typename FieldT>   
    std::vector<FieldT> additive_IFFT(const std::vector<FieldT> &poly_coeffs, 
                                    const size_t domain_dim, const size_t shift_dim);

    template<size_t K, typename FieldT>
    std::vector<FieldT> additive_IFFT_radix2k(const std::vector<FieldT> &evals,
                                              const size_t domain_dim, const size_t shift_dim);

    /// OpenMP parallel \ref additive_IFFT. Without \c _OPENMP, forwards to serial.
    template<typename FieldT>
    std::vector<FieldT> additive_IFFT_parallel(const std::vector<FieldT> &evals,
                                               const size_t domain_dim, const size_t shift_dim);

    /// OpenMP parallel \ref additive_IFFT_radix2k. Without \c _OPENMP, forwards to serial.
    template<size_t K, typename FieldT>
    std::vector<FieldT> additive_IFFT_radix2k_parallel(const std::vector<FieldT> &evals,
                                                        const size_t domain_dim, const size_t shift_dim);
}

#include "LCH/utils.tcc"
#include "LCH/fft.tcc"
#endif // ADDITIVE_FFT_LCH_HPP_
