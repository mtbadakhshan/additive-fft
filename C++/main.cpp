#include "libiop/algebra/subspace.hpp"
#include "libiop/algebra/utils.hpp"
#include "libiop/fft.hpp"
#include <algorithm>
#include <cstddef>
#include <libff/algebra/fields/binary/gf256.hpp>
#include <libff/algebra/fields/binary/gf128.hpp>
#include <libff/algebra/fields/binary/gf64.hpp>
#include <libff/algebra/fields/binary/gf32.hpp>
#include <libff/common/utils.hpp>
#include <ostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cmath>   // for std::sqrt

#include "utils/utils.hpp"
#include "Cantor/fft.hpp"
#include "Gao/fft.hpp"


template<typename FieldT>
void chrono_timing_libiop(const size_t iterations, const size_t warm_up_iter, const size_t m_min, const size_t m_max){
    const std::string filename = "chrono_output.csv";
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file for writing!" << std::endl;
        return;
    }
    size_t num_acc_sections = 5;
    // Write header
    file << "m,Initialization,Scaling,TaylorEx,BasisComputation,Bitreverse,Span,Aggregate\n";
    for (size_t m = m_min; m <= m_max; ++m){
        std::cout << "m = " << m << ", Start testing!\n";
        std::vector<double> mean_timing;
        std::vector<std::vector<double>> all_timings;
        for (size_t i = 0; i < iterations; ++i){
            //Domain and polynomial Creation
            libiop::field_subset<FieldT> domain = libiop::field_subset<FieldT>(libiop::affine_subspace<FieldT>::random_affine_subspace(m));
            std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
            const std::vector<FieldT> gao_result = libiop::additive_FFT<FieldT>(poly_coeffs, domain.subspace(), all_timings);
        }
        
        for (size_t part_cnt = 0; part_cnt < all_timings[0].size(); ++part_cnt){
            double sum = 0;
            for (size_t run = warm_up_iter; run < iterations; ++run){
                sum += all_timings[run][part_cnt];
            }
            mean_timing.push_back(sum/iterations);
        }

        file << m << "," <<  mean_timing[0];
        double scaling=0, taylorex=0, basis_computation=0;
        for (size_t i = 1; i < 3 * m + 1; i+=3 ){
            scaling += mean_timing[i];
            taylorex += mean_timing[i+1];
            basis_computation += mean_timing[i+2];
            // std::cout<< mean_timing[i] << ", " << mean_timing[i+1] << ", " << mean_timing[i+2] << std::endl;
        }
        file << "," <<  scaling;
        file << "," <<  taylorex;
        file << "," <<  basis_computation;
        file << "," <<  mean_timing[3 * m + 1]; // Bitreversal

        double span=0, aggregate=0;
        for (size_t i = 3 * m + 2; i < 5 * m + 2; i+=2 ){
            span += mean_timing[i];
            aggregate += mean_timing[i+1];
            // std::cout<< mean_timing[i] << ", " << mean_timing[i+1] << std::endl;
        }
        file << "," <<  span;
        file << "," <<  aggregate;
        file << std::endl;
    }
}

int main()
{
    typedef libff::gf256 FieldT;
    chrono_timing_libiop<FieldT>(2, 0, 5, 20);

    return 0;
}
