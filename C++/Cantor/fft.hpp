#ifndef ADDITIVE_FFT_CANTOR_HPP_
#define ADDITIVE_FFT_CANTOR_HPP_

#include <cstddef>
#include <vector>

#include "libiop/algebra/field_subset.hpp"
#include "libiop/algebra/subspace.hpp"

namespace cantor {

template<typename FieldT>
struct PreComputedValues {
    size_t dimension;
    std::vector<std::vector<size_t>> nz_S_vecs;
    std::vector<FieldT> mult_factor_vec;

    PreComputedValues(const size_t m){
        dimension = m;
        nz_S_vecs.reserve(m);
        mult_factor_vec.reserve((1<<m)-1);
    }
};

template<typename FieldT>
PreComputedValues<FieldT> pre_computation(const libiop::affine_subspace<FieldT> &domain);

/// Cantor-combination table path (hard-coded \c s_i / \c n_terms), same convention as LCH
/// \c (domain_dim, shift_dim) entry points.
template<typename FieldT>
std::vector<FieldT> additive_FFT(const std::vector<FieldT> &poly_coeffs,
                                const size_t domain_dim,
                                const size_t shift_dim);

template<typename FieldT>
std::vector<FieldT> additive_IFFT(const std::vector<FieldT> &evals,
                                  const size_t domain_dim,
                                  const size_t shift_dim);

template<typename FieldT>
std::vector<FieldT> additive_FFT(const std::vector<FieldT> &poly_coeffs,
                                        const PreComputedValues<FieldT> &values);

template<typename FieldT>
std::vector<FieldT> additive_IFFT(const std::vector<FieldT> &evals,
                                        const PreComputedValues<FieldT> &values);

template<typename FieldT>
std::vector<FieldT> additive_FFT(const std::vector<FieldT> &poly_coeffs,
                                        const libiop::affine_subspace<FieldT> &domain);

/// OpenMP Strategy A over modules for \ref additive_FFT(poly, domain_dim, shift_dim).
template<typename FieldT>
std::vector<FieldT> additive_FFT_parallel(const std::vector<FieldT> &poly_coeffs,
                                          const size_t domain_dim,
                                          const size_t shift_dim);

/// OpenMP Strategy A over modules (radix-2 Cantor FFT). Requires OpenMP at
/// link time for speedup; without \c _OPENMP, forwards to \ref additive_FFT.
template<typename FieldT>
std::vector<FieldT> additive_FFT_parallel(const std::vector<FieldT> &poly_coeffs,
                                          const libiop::affine_subspace<FieldT> &domain);

/// Radix fused paths on \c (domain_dim, shift_dim): same Cantor-combination chart as
/// \ref additive_FFT(poly, domain_dim, shift_dim). For the affine-subspace chart, use the
/// \c libiop::affine_subspace overloads.
template<typename FieldT>
std::vector<FieldT> additive_FFT_radix4(const std::vector<FieldT> &poly_coeffs,
                                        const size_t domain_dim,
                                        const size_t shift_dim);

template<typename FieldT>
std::vector<FieldT> additive_FFT_radix4(const std::vector<FieldT> &poly_coeffs,
                                        const libiop::affine_subspace<FieldT> &domain);

template<size_t K, typename FieldT>
std::vector<FieldT> additive_FFT_radix2k(const std::vector<FieldT> &poly_coeffs,
                                         const size_t domain_dim,
                                         const size_t shift_dim);

template<size_t K, typename FieldT>
std::vector<FieldT> additive_FFT_radix2k(const std::vector<FieldT> &poly_coeffs,
                                         const libiop::affine_subspace<FieldT> &domain);

/// OpenMP Strategy A for \ref additive_FFT_radix2k(poly, domain_dim, shift_dim) on the
/// combination table path. Without \c _OPENMP, forwards to the serial implementation.
template<size_t K, typename FieldT>
std::vector<FieldT> additive_FFT_radix2k_parallel(const std::vector<FieldT> &poly_coeffs,
                                                  const size_t domain_dim,
                                                  const size_t shift_dim);

template<size_t K, typename FieldT>
std::vector<FieldT> additive_FFT_radix2k_parallel(const std::vector<FieldT> &poly_coeffs,
                                                  const libiop::affine_subspace<FieldT> &domain);

template<typename FieldT>
std::vector<FieldT> additive_IFFT(const std::vector<FieldT> &evals,
                                        const libiop::affine_subspace<FieldT> &domain);

} // namespace cantor

#include "Cantor/fft.tcc"
#include "Cantor/fft_hc.tcc"
#include "Cantor/fft-radix.tcc"
#endif // ADDITIVE_FFT_CANTOR_HPP_
