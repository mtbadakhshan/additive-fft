#ifndef ADDITIVE_FFT_UTILS_HPP_
#define ADDITIVE_FFT_UTILS_HPP_

#include <vector>


template<typename FieldT>
bool field_trace_binary(const FieldT &element);

template<typename FieldT>
std::vector<FieldT> cantor_basis(size_t m);

#include "utils/utils.tcc"

#endif // ADDITIVE_FFT_UTILS_HPP_
