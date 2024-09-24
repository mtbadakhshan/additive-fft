#include "libiop/algebra/subspace.hpp"
#include "libiop/algebra/utils.hpp"
#include "libiop/fft.hpp"
#include <iostream>
#include <libff/algebra/fields/binary/gf64.hpp>
#include <libff/algebra/fields/binary/gf32.hpp>
#include <libff/common/utils.hpp>
#include <vector>

#include "utils/utils.hpp"

template <typename T>
bool check_equal(const std::vector<T> &v1, const std::vector<T> &v2)
{
    for (size_t i = 0; i < v1.size(); ++i)
    {
        if (v1[i] != v2[i])
        {
            return false;
        }
    }
    return true;
}

template <typename T>
void my_print_vector(const std::vector<T> &v)
{
    std::cout <<"{ ";
    for (auto const& elem : v)
    {
        std::cout << elem << " ";
    }
    std::cout << "}" << std::endl;
}

int main()
{
    std::cout << "Start testing!\n";
    typedef libff::gf32 FieldT;

    std::vector<FieldT> result(cantor_basis<FieldT>(5));
    my_print_vector(result);

    // libff::print_vector<FieldT>(static_cast<std::vector<FieldT>>(result));



    // FieldT element = FieldT::random_element();
    // std::cout << "element.extension_degree() " << element.extension_degree() << std::endl;
    // std::cout << "field_trace_binary of " << element << " is :" << field_trace_binary(element) << std::endl;

    if (false)
    {
        for (size_t m = 1; m <= 11; ++m)
        {

            std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
            libiop::field_subset<FieldT> domain = libiop::field_subset<FieldT>(
                libiop::affine_subspace<FieldT>::random_affine_subspace(m));

            std::cout << "domain dimension: " << domain.dimension() << std::endl;
            std::cout << "domain size: " << domain.num_elements() << std::endl;

            /* Additive equals naive */
            const std::vector<FieldT> naive_result =
                libiop::naive_FFT<FieldT>(poly_coeffs, domain);
            const std::vector<FieldT> additive_result =
                libiop::additive_FFT<FieldT>(poly_coeffs, domain.subspace());

            std::cout << "m = " << m << std::endl;
            std::cout << "naive_result[0] = " << naive_result[0] << std::endl;
            std::cout << "additive_result[0] = " << additive_result[0] << std::endl;

            std::cout << "Equality check" << check_equal<FieldT>(naive_result, naive_result) << std::endl;
        }
    }
    return 0;
}
