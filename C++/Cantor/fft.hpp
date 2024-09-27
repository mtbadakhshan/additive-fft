#ifndef ADDITIVE_FFT_CANTOR_HPP_
#define ADDITIVE_FFT_CANTOR_HPP_

#include <vector>

#include "libiop/algebra/field_subset.hpp"
#include "libiop/algebra/subspace.hpp"

template<typename FieldT>
std::vector<FieldT> cantor_additive_FFT(const std::vector<FieldT> &poly_coeffs,
                                        const libiop::affine_subspace<FieldT> &domain);

#include "Cantor/fft.tcc"
#endif // ADDITIVE_FFT_CANTOR_HPP_
