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
#include <vector>
#include <chrono>
#include <cmath>   // for std::sqrt

#include "utils/utils.hpp"
#include "Cantor/fft.hpp"
#include "LCH/fft.hpp"
#include "Gao/fft.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

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

/* This test uses different domains; hence, the timing report might be more accurate due to minimizing CPU's caching*/
void Cantor_FFT_Test(){
    typedef libff::gf256 FieldT;
    size_t m = 15;
    size_t N_test = 100;
    std::cout << "m = " << m << ", N_test = " << N_test << ", F = GF(2^"<<FieldT::extension_degree()<<")" << std::endl;
    std::vector<double> durations_Cantor(N_test);
    std::vector<double> durations_Gao(N_test);

    double Cantor_sum = 0;
    double Gao_sum = 0;

    std::vector<FieldT> basis(cantor_basis<FieldT>(m));

    for (size_t i=0; i<N_test; ++i){
        std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);

        // Cantor test
        libiop::field_subset<FieldT> cantor_domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
        auto start = std::chrono::high_resolution_clock::now();
        const std::vector<FieldT> cantor_result = cantor::additive_FFT<FieldT>(poly_coeffs, cantor_domain.subspace());
        auto stop = std::chrono::high_resolution_clock::now();
        durations_Cantor[i] = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
        Cantor_sum += durations_Cantor[i];

        // Gao test
        // libiop::field_subset<FieldT> gao_domain = libiop::field_subset<FieldT>(libiop::affine_subspace<FieldT>::random_affine_subspace(m));
        libiop::field_subset<FieldT> gao_domain = cantor_domain;
        start = std::chrono::high_resolution_clock::now();
        const std::vector<FieldT> additive_result = libiop::additive_FFT<FieldT>(poly_coeffs, gao_domain.subspace());    // my_print_vector<FieldT>(cantor_result);
        stop = std::chrono::high_resolution_clock::now();
        durations_Gao[i] = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
        Gao_sum += durations_Gao[i];

        std::cout << "Equality check: " << (check_equal<FieldT>(additive_result, cantor_result) ? "\033[1;32mPass\033[0m" : "\033[1;31mFail\033[0m")  << std::endl;

    }
    
    double cantor_variance_sum = 0.0;
    for (double duration : durations_Cantor) {
        cantor_variance_sum += (duration - (Cantor_sum / N_test)) * (duration - (Cantor_sum / N_test));
    }

    double gao_variance_sum = 0.0;
    for (double duration : durations_Gao) {
        gao_variance_sum += (duration - (Gao_sum / N_test)) * (duration - (Gao_sum / N_test));
    }

    std::cout << "Cantor's average duration: " << Cantor_sum / N_test << " ms \t Standard deviations: " << std::sqrt(cantor_variance_sum / N_test) << std::endl;
    std::cout << "Gao's    average duration: " << Gao_sum    / N_test << " ms \t Standard deviations: " << std::sqrt(gao_variance_sum    / N_test) << std::endl;

}

/* This test is primarely for correctness check. We use the same domain; hence, the timing report might not be accurate due to CPU's caching*/
void Gao_CO_FFT_Test(){
    typedef libff::gf128 FieldT;
    size_t m = 20;
    std::cout << "m = " << m << ", Start testing!\n";

    std::vector<FieldT> basis(cantor_basis<FieldT>(m));
    std::reverse(basis.begin(), basis.end());

    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    gao::PreComputedValues_CO<FieldT> values = gao::pre_computation_co<FieldT>(domain.subspace());
    
    std::chrono::time_point<std::chrono::high_resolution_clock> start, stop; 

    std::cout << "Entering gao_co:\n";
    start = std::chrono::high_resolution_clock::now();
    const std::vector<FieldT> gao_co_result = gao::additive_FFT_CO<FieldT>(poly_coeffs, domain.subspace());
    stop = std::chrono::high_resolution_clock::now();
    auto duration_gao_co = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << "Entering gao_co precomputed:\n";
    start = std::chrono::high_resolution_clock::now();
    const std::vector<FieldT> gao_co_result_precmp = gao::additive_FFT_CO<FieldT>(poly_coeffs, values);
    stop = std::chrono::high_resolution_clock::now();
    auto duration_gao_co_precmp = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);


    std::cout << "Entering libop:\n";
    start = std::chrono::high_resolution_clock::now();
    const std::vector<FieldT> additive_result = libiop::additive_FFT<FieldT>(poly_coeffs, domain.subspace()); 
    stop = std::chrono::high_resolution_clock::now();
    auto duration_Gao = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << "Equality check 1: " << (check_equal<FieldT>(additive_result, gao_co_result) ? "\033[1;32mPass\033[0m" : "\033[1;31mFail\033[0m") << std::endl;
    std::cout << "Equality check 2: " << (check_equal<FieldT>(additive_result, gao_co_result_precmp) ? "\033[1;32mPass\033[0m" : "\033[1;31mFail\033[0m") << std::endl;
    
    std::cout << "Gao_CO's Duration (pre-computed): " << duration_gao_co_precmp.count() << " ms" << std::endl;
    std::cout << "Gao_CO's Duration: " << duration_gao_co.count() << " ms" << std::endl;
    std::cout << "Gao's Duration: " << duration_Gao.count() << " ms" << std::endl;

}

/* This test is primarely for correctness check. We use the same domain; hence, the timing report might not be accurate due to CPU's caching*/
void Cantor_FFT_PreComputation_Test(){
    typedef libff::gf128 FieldT;
    size_t m = 15;
    std::cout << "m = " << m << ", Start testing!\n";

    size_t N_test = 10;

    for (size_t i=0; i<N_test; ++i){
    std::cout<<"test "<< i+1 << std::endl;
    std::vector<FieldT> basis(cantor_basis<FieldT>(m));
    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    
    cantor::PreComputedValues<FieldT> values = cantor::pre_computation(domain.subspace());

    auto start = std::chrono::high_resolution_clock::now();
    const std::vector<FieldT> cantor_precmp_result = cantor::additive_FFT<FieldT>(poly_coeffs, values);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration_Cantor_precmp = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    start = std::chrono::high_resolution_clock::now();
    const std::vector<FieldT> cantor_result = cantor::additive_FFT<FieldT>(poly_coeffs, domain.subspace());
    stop = std::chrono::high_resolution_clock::now();
    auto duration_Cantor = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << "Equality check: " << (check_equal<FieldT>(cantor_precmp_result, cantor_result) ? "\033[1;32mPass\033[0m" : "\033[1;31mFail\033[0m")  << std::endl;
    std::cout << "Cantor's Duration (pre computed): " << duration_Cantor_precmp.count()<< " ms" << std::endl;
    std::cout << "Cantor's Duration: " << duration_Cantor.count()<< " ms" << std::endl;

    }
}

/* This test is primarely for correctness check. We use the same domain; hence, the timing report might not be accurate due to CPU's caching*/
void Gao_FFT_PreComputation_Test(){
    typedef libff::gf256 FieldT;
    size_t m = 18;
    std::cout << "m = " << m << ", Start testing!\n";

    libiop::field_subset<FieldT> domain = libiop::field_subset<FieldT>(libiop::affine_subspace<FieldT>::random_affine_subspace(m));
    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);

    auto start = std::chrono::high_resolution_clock::now();
    const std::vector<FieldT> gao_result = libiop::additive_FFT<FieldT>(poly_coeffs, domain.subspace());
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration_Gao = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    gao::PreComputedValues_Level1<FieldT> values_lvl1 = gao::pre_computation_lvl1<FieldT>(domain.subspace());
    start = std::chrono::high_resolution_clock::now();
    const std::vector<FieldT> gao_precmp_result_lvl1 = gao::additive_FFT<FieldT>(poly_coeffs, values_lvl1);
    stop = std::chrono::high_resolution_clock::now();
    auto duration_gao_precmp_lvl1 = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    gao::PreComputedValues_Level2<FieldT> values_lvl2 = gao::pre_computation_lvl2<FieldT>(domain.subspace());
    start = std::chrono::high_resolution_clock::now();
    const std::vector<FieldT> gao_precmp_result_lvl2 = gao::additive_FFT<FieldT>(poly_coeffs, values_lvl2);
    stop = std::chrono::high_resolution_clock::now();
    auto duration_gao_precmp_lvl2 = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);


    std::cout << "Level1 Equality check: " << (check_equal<FieldT>(gao_precmp_result_lvl1, gao_result) ? "\033[1;32mPass\033[0m" : "\033[1;31mFail\033[0m")  << std::endl;
    std::cout << "Level2 Equality check: " << (check_equal<FieldT>(gao_precmp_result_lvl2, gao_result) ? "\033[1;32mPass\033[0m" : "\033[1;31mFail\033[0m")  << std::endl;
    std::cout << "Gao's Duration (pre computed level1): " << duration_gao_precmp_lvl1.count()<< " ms" << std::endl;
    std::cout << "Gao's Duration (pre computed level2): " << duration_gao_precmp_lvl2.count()<< " ms" << std::endl;
    std::cout << "Gao's Duration: " << duration_Gao.count()<< " ms" << std::endl;

}

void Valgrid_libiop_test(){
    typedef libff::gf256 FieldT;
    size_t m = 20;
    std::cout << "m = " << m << ", Start testing!\n";

    //Domain Creation
    libiop::field_subset<FieldT> domain = libiop::field_subset<FieldT>(libiop::affine_subspace<FieldT>::random_affine_subspace(m));
    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);

    const std::vector<FieldT> gao_result = libiop::additive_FFT<FieldT>(poly_coeffs, domain.subspace());
}

void Valgrid_cantor_test(){
    typedef libff::gf256 FieldT;
    size_t m = 10;
    std::cout << "m = " << m << ", Start testing!\n";

    //Cantor Special Basis Domain Creation
    std::vector<FieldT> basis(cantor_basis<FieldT>(m));
    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
    
    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);

    const std::vector<FieldT> cantor_result = cantor::additive_FFT<FieldT>(poly_coeffs, domain.subspace());
}


void Valgrid_cantorPC_test(){
    typedef libff::gf256 FieldT;
    size_t m = 15;
    std::cout << "m = " << m << ", Start testing!\n";

    //Cantor Special Basis Domain Creation
    std::vector<FieldT> basis(cantor_basis<FieldT>(m));
    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};

    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    cantor::PreComputedValues<FieldT> values = cantor::pre_computation(domain.subspace());

    const std::vector<FieldT> cantor_result = cantor::additive_FFT<FieldT>(poly_coeffs, values);
    const std::vector<FieldT> poly_coeffs_computed = cantor::additive_IFFT<FieldT>(cantor_result, values);

    std::cout << "Equality check: " << (check_equal<FieldT>(poly_coeffs_computed, poly_coeffs) ? "\033[1;32mPass\033[0m" : "\033[1;31mFail\033[0m")  << std::endl;


}


void Pre_Compute_Cantor_basis(){
    typedef libff::gf192 FieldT;
    int N = 32;

    std::vector<FieldT> basis = cantor_basis<FieldT>(N);
    std::cout <<"cantor_in_gf"<<FieldT::extension_degree()<<"["<<N<<"]["<<FieldT::extension_degree()/64<<"] = {\n"; 
    for (int i = 0; i < N; i++){
        std::vector<uint64_t> basis_in_words = basis[i].to_words();
    std::cout <<"{ "; 
        for (int j = 0; j < basis_in_words.size() - 1; j++){
            std::cout<< basis_in_words[j] << ", ";
        }
        std::cout<< basis_in_words.back() << " }, \n";
    }
    std::cout <<"};\n"; 


}

#include "Cantor/cantor_basis.hpp"
void check_the_basis_element_order(){
    typedef libff::gf128 FieldT_128;
    int N = 32;
    
    for (int i = 0; i < N; ++i){
        FieldT_128 element = FieldT_128(cantor::cantor_in_gf2to128[i][1], cantor::cantor_in_gf2to128[i][0]);
        std::cout << "e = " << element << " | ";
        std::cout << "e^2 = "  << (element^2) << " | ";
        std::cout << "e^3 = "  << (element^3) << " | ";
        std::cout << "e^4 = "  << (element^4) << " | ";
        std::cout << "e^5 = "  << (element^2) << std::endl;
    }
}

void test_LCH(){
    typedef libff::gf192 FieldT;
    unsigned m = 17;
    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1656);
    const std::vector<FieldT> cantor_result = cantor::additive_FFT_hc<FieldT>(poly_coeffs, m, m);
    const std::vector<FieldT> lch_result = lch::additive_FFT(poly_coeffs, m, m);

    std::cout << "Equality check: " << (check_equal<FieldT>(lch_result, cantor_result) ? "\033[1;32mPass\033[0m" : "\033[1;31mFail\033[0m")  << std::endl;
}

// Correctness check: lch::additive_FFT_radix2k<K> against the radix-2 baseline
// for K in {1, 2, 3, 4, 5}, exercising even/odd log_n cases.
void test_LCH_radix2k_correctness(){
    typedef libff::gf256 FieldT;
    const std::vector<unsigned> ms = {6, 8, 9, 10, 12, 14, 15, 16, 18, 20};
    bool all_passed = true;

    for (unsigned m : ms) {
        std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
        const auto baseline = lch::additive_FFT<FieldT>(poly_coeffs, m, m);

        const auto r2k1 = lch::additive_FFT_radix2k<1, FieldT>(poly_coeffs, m, m);
        const auto r2k2 = lch::additive_FFT_radix2k<2, FieldT>(poly_coeffs, m, m);
        const auto r2k3 = lch::additive_FFT_radix2k<3, FieldT>(poly_coeffs, m, m);
        const auto r2k4 = lch::additive_FFT_radix2k<4, FieldT>(poly_coeffs, m, m);
        const auto r2k5 = lch::additive_FFT_radix2k<5, FieldT>(poly_coeffs, m, m);

        const bool ok1  = check_equal<FieldT>(baseline, r2k1);
        const bool ok2  = check_equal<FieldT>(baseline, r2k2);
        const bool ok3  = check_equal<FieldT>(baseline, r2k3);
        const bool ok4  = check_equal<FieldT>(baseline, r2k4);
        const bool ok5  = check_equal<FieldT>(baseline, r2k5);

        std::cout << "m=" << m
                  << "  K=1: "      << (ok1  ? "PASS" : "FAIL")
                  << "  K=2: "      << (ok2  ? "PASS" : "FAIL")
                  << "  K=3: "      << (ok3  ? "PASS" : "FAIL")
                  << "  K=4: "      << (ok4  ? "PASS" : "FAIL")
                  << "  K=5: "      << (ok5  ? "PASS" : "FAIL")
                  << std::endl;

        all_passed = all_passed && ok1 && ok2 && ok3 && ok4 && ok5;
    }
    std::cout << (all_passed ? "ALL PASS" : "SOME FAILED") << std::endl;
}

// LCH parallel butterflies vs radix-2 baseline.
void test_LCH_parallel_correctness()
{
    typedef libff::gf256 FieldT;
    const std::vector<unsigned> ms = {6, 8, 9, 10, 12, 14, 15, 16};
    bool all_passed = true;

    for (unsigned m : ms) {
        std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
        const auto baseline = lch::additive_FFT<FieldT>(poly_coeffs, m, m);

#ifdef _OPENMP
        const std::vector<int> thread_counts = {1, 2, 4, 8};
        for (int t : thread_counts) {
            omp_set_num_threads(t);
            const auto pr2  = lch::additive_FFT_parallel<FieldT>(poly_coeffs, m, m);
            const auto pk2  = lch::additive_FFT_radix2k_parallel<2, FieldT>(poly_coeffs, m, m);
            const auto pk3  = lch::additive_FFT_radix2k_parallel<3, FieldT>(poly_coeffs, m, m);
            const auto pk4  = lch::additive_FFT_radix2k_parallel<4, FieldT>(poly_coeffs, m, m);
            const bool ok2  = check_equal<FieldT>(baseline, pr2);
            const bool okk2 = check_equal<FieldT>(baseline, pk2);
            const bool okk3 = check_equal<FieldT>(baseline, pk3);
            const bool okk4 = check_equal<FieldT>(baseline, pk4);
            all_passed      = all_passed && ok2 && okk2 && okk3 && okk4;
            std::cout << "m=" << m << "  threads=" << t << "  r2_par:" << (ok2 ? "PASS" : "FAIL")
                      << "  r2k2_par:" << (okk2 ? "PASS" : "FAIL")
                      << "  r2k3_par:" << (okk3 ? "PASS" : "FAIL") << "  r2k4_par:" << (okk4 ? "PASS" : "FAIL")
                      << std::endl;
        }
#else
        const auto pr2 = lch::additive_FFT_parallel<FieldT>(poly_coeffs, m, m);
        const auto pk4 = lch::additive_FFT_radix2k_parallel<4, FieldT>(poly_coeffs, m, m);
        const bool ok  = check_equal<FieldT>(baseline, pr2) && check_equal<FieldT>(baseline, pk4);
        all_passed     = all_passed && ok;
        std::cout << "m=" << m << "  (no OpenMP) parallel aliases: " << (ok ? "PASS" : "FAIL") << std::endl;
#endif
    }

    std::cout << (all_passed ? "LCH parallel ALL PASS" : "SOME FAILED") << std::endl;
}

// Correctness check: lch::additive_FFT_radix2k<2> against the radix-2 baseline
// across several (m, n_poly, shift_dim) configurations, including:
//   * odd log_n (exercises the residual radix-2 stage),
//   * poly_coeffs.size() < 2^m (exercises the post-conversion replication),
//   * shift_dim = 0 (linear subspace, no affine shift),
//   * shift_dim > m (affine shift outside the FFT domain).
void test_LCH_radix2k2_extended_correctness(){
    typedef libff::gf192 FieldT;

    struct Case { unsigned m; size_t n_poly; size_t shift_dim; };
    const std::vector<Case> cases = {
        {6,  1u<<6,  6},
        {7,  1u<<7,  7},
        {8,  1u<<8,  8},
        {9,  1u<<9,  9},
        {10, 1u<<10, 10},
        {12, 1u<<12, 12},
        {14, 1u<<14, 14},
        {15, 1u<<15, 15},
        {16, 1u<<16, 16},
        {17, 1656,   17},                        // log_poly_terms = 11 (odd)
        {18, 1u<<18, 18},
        {12, 1u<<12, 0 },                        // pure linear subspace
        {12, 1u<<12, 20},                        // shift_dim > m
        {13, 4321,   25},                        // poly < domain, odd log_n, large shift_dim
    };

    bool all_passed = true;
    for (const auto& c : cases) {
        std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(c.n_poly);

        const auto baseline = lch::additive_FFT<FieldT>(poly_coeffs, c.m, c.shift_dim);
        const auto r2k2     = lch::additive_FFT_radix2k<2, FieldT>(poly_coeffs, c.m, c.shift_dim);
        const bool ok       = check_equal<FieldT>(baseline, r2k2);
        all_passed = all_passed && ok;

        std::cout << "m=" << c.m
                  << "  n_poly=" << c.n_poly
                  << "  shift_dim=" << c.shift_dim
                  << "  radix2k<2> vs radix-2: "
                  << (ok ? "\033[1;32mPass\033[0m" : "\033[1;31mFail\033[0m")
                  << std::endl;
    }

    std::cout << (all_passed ? "ALL PASS" : "SOME FAILED") << std::endl;
}

// Correctness check: cantor::additive_FFT_radix2k<K> against the radix-2 baseline
// across several K and m values.
void test_radix2k_correctness(){
    typedef libff::gf256 FieldT;
    const std::vector<size_t> ms = {2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18};
    bool all_passed = true;

    for (size_t m : ms) {
        std::vector<FieldT> basis(cantor_basis<FieldT>(m));
        libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
        std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);

        const auto baseline = cantor::additive_FFT<FieldT>(poly_coeffs, domain.subspace());

        const auto r2k1 = cantor::additive_FFT_radix2k<1, FieldT>(poly_coeffs, domain.subspace());
        const bool ok1 = check_equal<FieldT>(baseline, r2k1);

        const auto r2k2 = cantor::additive_FFT_radix2k<2, FieldT>(poly_coeffs, domain.subspace());
        const bool ok2 = check_equal<FieldT>(baseline, r2k2);

        const auto r2k3 = cantor::additive_FFT_radix2k<3, FieldT>(poly_coeffs, domain.subspace());
        const bool ok3 = check_equal<FieldT>(baseline, r2k3);

        const auto r2k4 = cantor::additive_FFT_radix2k<4, FieldT>(poly_coeffs, domain.subspace());
        const bool ok4 = check_equal<FieldT>(baseline, r2k4);

        std::cout << "m=" << m
                  << "  K=1: " << (ok1 ? "PASS" : "FAIL")
                  << "  K=2: " << (ok2 ? "PASS" : "FAIL")
                  << "  K=3: " << (ok3 ? "PASS" : "FAIL")
                  << "  K=4: " << (ok4 ? "PASS" : "FAIL")
                  << std::endl;

        all_passed = all_passed && ok1 && ok2 && ok3 && ok4;
    }

    std::cout << (all_passed ? "ALL PASS" : "SOME FAILED") << std::endl;
}

// cantor::additive_FFT_parallel vs radix-2 baseline.
void test_additive_fft_parallel_correctness()
{
    typedef libff::gf256 FieldT;
    const std::vector<size_t> ms = {4, 5, 6, 8, 10, 12, 14, 16};
    bool all_passed = true;

    for (size_t m : ms) {
        std::vector<FieldT> basis(cantor_basis<FieldT>(m));
        libiop::field_subset<FieldT> domain{
            libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
        std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);

        const auto baseline = cantor::additive_FFT<FieldT>(poly_coeffs, domain.subspace());

#ifdef _OPENMP
        const std::vector<int> thread_counts = {1, 2, 4, 8};
        for (int t : thread_counts) {
            omp_set_num_threads(t);
            const auto par = cantor::additive_FFT_parallel<FieldT>(poly_coeffs, domain.subspace());
            const bool ok = check_equal<FieldT>(baseline, par);
            all_passed = all_passed && ok;
            std::cout << "m=" << m << "  threads=" << t << "  parallel vs radix-2: "
                      << (ok ? "\033[1;32mPass\033[0m" : "\033[1;31mFail\033[0m") << std::endl;
        }
#else
        const auto par = cantor::additive_FFT_parallel<FieldT>(poly_coeffs, domain.subspace());
        const bool ok = check_equal<FieldT>(baseline, par);
        all_passed = all_passed && ok;
        std::cout << "m=" << m << "  (no OpenMP) parallel alias: "
                  << (ok ? "PASS" : "FAIL") << std::endl;
#endif
    }

    std::cout << (all_passed ? "additive_FFT_parallel ALL PASS" : "SOME FAILED") << std::endl;
}

// radix-2^K parallel vs radix-2 baseline (K in {2,3,4}).
void test_radix2k_parallel_correctness()
{
    typedef libff::gf256 FieldT;
    const std::vector<size_t> ms = {4, 5, 6, 8, 10, 11, 12, 13, 14, 16};
    bool all_passed = true;

    for (size_t m : ms) {
        std::vector<FieldT> basis(cantor_basis<FieldT>(m));
        libiop::field_subset<FieldT> domain{
            libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
        std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);

        const auto baseline = cantor::additive_FFT<FieldT>(poly_coeffs, domain.subspace());

#ifdef _OPENMP
        const std::vector<int> thread_counts = {1, 2, 4, 8};
        for (int t : thread_counts) {
            omp_set_num_threads(t);
            const auto p2 = cantor::additive_FFT_radix2k_parallel<2, FieldT>(poly_coeffs, domain.subspace());
            const auto p3 = cantor::additive_FFT_radix2k_parallel<3, FieldT>(poly_coeffs, domain.subspace());
            const auto p4 = cantor::additive_FFT_radix2k_parallel<4, FieldT>(poly_coeffs, domain.subspace());
            const bool ok2 = check_equal<FieldT>(baseline, p2);
            const bool ok3 = check_equal<FieldT>(baseline, p3);
            const bool ok4 = check_equal<FieldT>(baseline, p4);
            all_passed = all_passed && ok2 && ok3 && ok4;
            std::cout << "m=" << m << "  threads=" << t
                      << "  K=2 par: " << (ok2 ? "PASS" : "FAIL")
                      << "  K=3 par: " << (ok3 ? "PASS" : "FAIL")
                      << "  K=4 par: " << (ok4 ? "PASS" : "FAIL") << std::endl;
        }
#else
        const auto p2 = cantor::additive_FFT_radix2k_parallel<2, FieldT>(poly_coeffs, domain.subspace());
        const auto p3 = cantor::additive_FFT_radix2k_parallel<3, FieldT>(poly_coeffs, domain.subspace());
        const auto p4 = cantor::additive_FFT_radix2k_parallel<4, FieldT>(poly_coeffs, domain.subspace());
        const bool ok2 = check_equal<FieldT>(baseline, p2);
        const bool ok3 = check_equal<FieldT>(baseline, p3);
        const bool ok4 = check_equal<FieldT>(baseline, p4);
        all_passed = all_passed && ok2 && ok3 && ok4;
        std::cout << "m=" << m << "  (no OpenMP) radix2k_parallel K=2,3,4: "
                  << (ok2 && ok3 && ok4 ? "PASS" : "FAIL") << std::endl;
#endif
    }

    std::cout << (all_passed ? "radix2k_parallel ALL PASS" : "SOME FAILED") << std::endl;
}

/// Cantor table path: \c additive_FFT_parallel(poly, m, shift_dim) vs serial baseline.
void test_cantor_precmp_parallel_correctness()
{
    typedef libff::gf256 FieldT;
    const std::vector<unsigned> ms = {4, 5, 6, 8, 10, 12, 14, 16};
    const size_t                   shift_dim = 31;
    bool                           all_passed = true;

    for (unsigned m : ms) {
        std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
        const auto          baseline    = cantor::additive_FFT<FieldT>(poly_coeffs, m, shift_dim);

#ifdef _OPENMP
        const std::vector<int> thread_counts = {1, 2, 4, 8};
        for (int t : thread_counts) {
            omp_set_num_threads(t);
            const auto par = cantor::additive_FFT_parallel<FieldT>(poly_coeffs, m, shift_dim);
            const bool ok  = check_equal<FieldT>(baseline, par);
            all_passed     = all_passed && ok;
            std::cout << "cantor precmp m=" << m << "  threads=" << t << "  parallel: "
                      << (ok ? "\033[1;32mPass\033[0m" : "\033[1;31mFail\033[0m") << std::endl;
        }
#else
        const auto par = cantor::additive_FFT_parallel<FieldT>(poly_coeffs, m, shift_dim);
        const bool ok  = check_equal<FieldT>(baseline, par);
        all_passed     = all_passed && ok;
        std::cout << "cantor precmp m=" << m << "  (no OpenMP) parallel alias: "
                  << (ok ? "PASS" : "FAIL") << std::endl;
#endif
    }

    std::cout << (all_passed ? "cantor precmp parallel ALL PASS" : "SOME FAILED") << std::endl;
}

/// Fused radix-2^K on the Cantor combination table path vs radix-2 baseline.
void test_cantor_precmp_radix2k_correctness()
{
    typedef libff::gf256 FieldT;
    const std::vector<unsigned> ms = {2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
    const size_t                 shift_dim = 31;
    bool                         all_passed = true;

    for (unsigned m : ms) {
        std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
        const auto          baseline    = cantor::additive_FFT<FieldT>(poly_coeffs, m, shift_dim);

        const auto r2k1 = cantor::additive_FFT_radix2k<1, FieldT>(poly_coeffs, m, shift_dim);
        const auto r2k2 = cantor::additive_FFT_radix2k<2, FieldT>(poly_coeffs, m, shift_dim);
        const auto r2k3 = cantor::additive_FFT_radix2k<3, FieldT>(poly_coeffs, m, shift_dim);
        const auto r2k4 = cantor::additive_FFT_radix2k<4, FieldT>(poly_coeffs, m, shift_dim);

        const bool ok1  = check_equal<FieldT>(baseline, r2k1);
        const bool ok2  = check_equal<FieldT>(baseline, r2k2);
        const bool ok3  = check_equal<FieldT>(baseline, r2k3);
        const bool ok4  = check_equal<FieldT>(baseline, r2k4);

        std::cout << "cantor precmp radix2k m=" << m << "  K=1:" << (ok1 ? "PASS" : "FAIL")
                  << "  K=2:" << (ok2 ? "PASS" : "FAIL") << "  K=3:" << (ok3 ? "PASS" : "FAIL")
                  << "  K=4:" << (ok4 ? "PASS" : "FAIL")
                  << std::endl;
        all_passed = all_passed && ok1 && ok2 && ok3 && ok4;
    }

    std::cout << (all_passed ? "cantor precmp radix2k ALL PASS" : "SOME FAILED") << std::endl;
}

/// Parallel fused radix-2^K on the Cantor table path vs serial baseline.
void test_cantor_precmp_radix2k_parallel_correctness()
{
    typedef libff::gf256 FieldT;
    const std::vector<unsigned> ms = {4, 5, 6, 8, 10, 12, 14, 16};
    const size_t                 shift_dim = 31;
    bool                         all_passed = true;

    for (unsigned m : ms) {
        std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
        const auto          baseline    = cantor::additive_FFT<FieldT>(poly_coeffs, m, shift_dim);

#ifdef _OPENMP
        const std::vector<int> thread_counts = {1, 2, 4, 8};
        for (int t : thread_counts) {
            omp_set_num_threads(t);
            const auto p2 = cantor::additive_FFT_radix2k_parallel<2, FieldT>(poly_coeffs, m, shift_dim);
            const auto p3 = cantor::additive_FFT_radix2k_parallel<3, FieldT>(poly_coeffs, m, shift_dim);
            const auto p4 = cantor::additive_FFT_radix2k_parallel<4, FieldT>(poly_coeffs, m, shift_dim);
            const bool ok2 = check_equal<FieldT>(baseline, p2);
            const bool ok3 = check_equal<FieldT>(baseline, p3);
            const bool ok4 = check_equal<FieldT>(baseline, p4);
            all_passed     = all_passed && ok2 && ok3 && ok4;
            std::cout << "cantor precmp radix2k par m=" << m << "  threads=" << t << "  K=2:"
                      << (ok2 ? "PASS" : "FAIL") << "  K=3:" << (ok3 ? "PASS" : "FAIL") << "  K=4:"
                      << (ok4 ? "PASS" : "FAIL") << std::endl;
        }
#else
        const auto p2 = cantor::additive_FFT_radix2k_parallel<2, FieldT>(poly_coeffs, m, shift_dim);
        const bool ok2 = check_equal<FieldT>(baseline, p2);
        all_passed     = all_passed && ok2;
        std::cout << "cantor precmp radix2k par m=" << m << "  (no OpenMP) K=2: " << (ok2 ? "PASS" : "FAIL")
                  << std::endl;
#endif
    }

    std::cout << (all_passed ? "cantor precmp radix2k_parallel ALL PASS" : "SOME FAILED") << std::endl;
}

/// FFT then IFFT round-trip on the Cantor combination-table path.
void test_cantor_precmp_ifft_roundtrip()
{
    typedef libff::gf256 FieldT;
    const std::vector<unsigned> ms = {6, 8, 10, 12, 14, 16};
    const size_t                 shift_dim = 31;
    bool                         all_passed = true;

    for (unsigned m : ms) {
        std::vector<FieldT> poly = libiop::random_vector<FieldT>(1ull << m);

        const auto eval_r2 = cantor::additive_FFT<FieldT>(poly, m, shift_dim);
        const auto rec_r2  = cantor::additive_IFFT<FieldT>(eval_r2, m, shift_dim);
        const bool ok_r2   = check_equal<FieldT>(poly, rec_r2);

        const auto eval_k2 = cantor::additive_FFT_radix2k<2, FieldT>(poly, m, shift_dim);
        const auto rec_k2  = cantor::additive_IFFT_radix2k<2, FieldT>(eval_k2, m, shift_dim);
        const bool ok_k2   = check_equal<FieldT>(poly, rec_k2);

        const auto eval_k3 = cantor::additive_FFT_radix2k<3, FieldT>(poly, m, shift_dim);
        const auto rec_k3  = cantor::additive_IFFT_radix2k<3, FieldT>(eval_k3, m, shift_dim);
        const bool ok_k3   = check_equal<FieldT>(poly, rec_k3);

        const auto eval_k4 = cantor::additive_FFT_radix2k<4, FieldT>(poly, m, shift_dim);
        const auto rec_k4  = cantor::additive_IFFT_radix2k<4, FieldT>(eval_k4, m, shift_dim);
        const bool ok_k4   = check_equal<FieldT>(poly, rec_k4);

        all_passed = all_passed && ok_r2 && ok_k2 && ok_k3 && ok_k4;
        std::cout << "cantor ifft roundtrip m=" << m << "  r2:" << (ok_r2 ? "PASS" : "FAIL")
                  << "  k2:" << (ok_k2 ? "PASS" : "FAIL") << "  k3:" << (ok_k3 ? "PASS" : "FAIL")
                  << "  k4:" << (ok_k4 ? "PASS" : "FAIL") << std::endl;

#ifdef _OPENMP
        const std::vector<int> thread_counts = {1, 2, 4, 8};
        for (int t : thread_counts) {
            omp_set_num_threads(t);
            const auto eval_p = cantor::additive_FFT_parallel<FieldT>(poly, m, shift_dim);
            const auto rec_p  = cantor::additive_IFFT_parallel<FieldT>(eval_p, m, shift_dim);
            const auto eval_pk2 =
                cantor::additive_FFT_radix2k_parallel<2, FieldT>(poly, m, shift_dim);
            const auto rec_pk2 =
                cantor::additive_IFFT_radix2k_parallel<2, FieldT>(eval_pk2, m, shift_dim);
            const bool ok_p   = check_equal<FieldT>(poly, rec_p);
            const bool ok_pk2 = check_equal<FieldT>(poly, rec_pk2);
            all_passed        = all_passed && ok_p && ok_pk2;
            std::cout << "  threads=" << t << "  par_r2:" << (ok_p ? "PASS" : "FAIL")
                      << "  par_k2:" << (ok_pk2 ? "PASS" : "FAIL") << std::endl;
        }
#endif
    }

    std::cout << (all_passed ? "cantor precmp ifft roundtrip ALL PASS" : "SOME FAILED") << std::endl;
}

/// FFT then IFFT round-trip for LCH (basis conversion included).
void test_lch_ifft_roundtrip()
{
    typedef libff::gf256 FieldT;
    const std::vector<unsigned> ms = {6, 8, 10, 12, 14, 16};
    bool                         all_passed = true;

    for (unsigned m : ms) {
        std::vector<FieldT> poly = libiop::random_vector<FieldT>(1ull << m);

        const auto eval_r2 = lch::additive_FFT<FieldT>(poly, m, m);
        const auto rec_r2  = lch::additive_IFFT<FieldT>(eval_r2, m, m);
        const bool ok_r2   = check_equal<FieldT>(poly, rec_r2);

        const auto eval_k2 = lch::additive_FFT_radix2k<2, FieldT>(poly, m, m);
        const auto rec_k2  = lch::additive_IFFT_radix2k<2, FieldT>(eval_k2, m, m);
        const bool ok_k2   = check_equal<FieldT>(poly, rec_k2);

        const auto eval_k3 = lch::additive_FFT_radix2k<3, FieldT>(poly, m, m);
        const auto rec_k3  = lch::additive_IFFT_radix2k<3, FieldT>(eval_k3, m, m);
        const bool ok_k3   = check_equal<FieldT>(poly, rec_k3);

        const auto eval_k5 = lch::additive_FFT_radix2k<5, FieldT>(poly, m, m);
        const auto rec_k5  = lch::additive_IFFT_radix2k<5, FieldT>(eval_k5, m, m);
        const bool ok_k5   = check_equal<FieldT>(poly, rec_k5);

        all_passed = all_passed && ok_r2 && ok_k2 && ok_k3 && ok_k5;
        std::cout << "lch ifft roundtrip m=" << m << "  r2:" << (ok_r2 ? "PASS" : "FAIL")
                  << "  k2:" << (ok_k2 ? "PASS" : "FAIL") << "  k3:" << (ok_k3 ? "PASS" : "FAIL")
                  << "  k5:" << (ok_k5 ? "PASS" : "FAIL") << std::endl;

#ifdef _OPENMP
        const std::vector<int> thread_counts = {1, 2, 4, 8};
        for (int t : thread_counts) {
            omp_set_num_threads(t);
            const auto eval_p = lch::additive_FFT_parallel<FieldT>(poly, m, m);
            const auto rec_p  = lch::additive_IFFT_parallel<FieldT>(eval_p, m, m);
            const auto eval_pk3 = lch::additive_FFT_radix2k_parallel<3, FieldT>(poly, m, m);
            const auto rec_pk3  = lch::additive_IFFT_radix2k_parallel<3, FieldT>(eval_pk3, m, m);
            const bool ok_p     = check_equal<FieldT>(poly, rec_p);
            const bool ok_pk3   = check_equal<FieldT>(poly, rec_pk3);
            all_passed          = all_passed && ok_p && ok_pk3;
            std::cout << "  threads=" << t << "  par_r2:" << (ok_p ? "PASS" : "FAIL")
                      << "  par_k3:" << (ok_pk3 ? "PASS" : "FAIL") << std::endl;
        }
#endif
    }

    std::cout << (all_passed ? "lch ifft roundtrip ALL PASS" : "SOME FAILED") << std::endl;
}

int main()
{
    // Pre_Compute_Cantor_basis();
    // Cantor_FFT_Test();
    // Gao_CO_FFT_Test();
    // Cantor_FFT_PreComputation_Test();
    // Gao_FFT_PreComputation_Test();

    // Valgrid_libiop_test();
    // Valgrid_cantorPC_test();
    // check_the_basis_element_order();

    // test_LCH();
    test_LCH_radix2k2_extended_correctness();
    test_LCH_radix2k_correctness();
    test_LCH_parallel_correctness();
    test_additive_fft_parallel_correctness();
    test_radix2k_parallel_correctness();
    test_cantor_precmp_parallel_correctness();
    test_cantor_precmp_radix2k_correctness();
    test_cantor_precmp_radix2k_parallel_correctness();
    test_cantor_precmp_ifft_roundtrip();
    test_lch_ifft_roundtrip();

    // test_radix2k_correctness();

    return 0;
}
