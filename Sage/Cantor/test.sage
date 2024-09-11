from time import time
import copy


load('fft.sage')
load('../utils/utils.sage')

def test_fft(a, FF, ext_degree):
    N_tests = 10
    DIRECT_EVAUATION_TEST = True
    
    direct_eval_time = [0] * N_tests
    cantors_fft_no_precmp_time = [0] * N_tests
    cantors_fft_with_precmp_time = [0] * N_tests

    m = 10
    print("Entering Pre-computation")
    W, nz_hdt_S, table = fft_precmp(a, m, ext_degree)
    print("Finished Pre-computation")

    if DIRECT_EVAUATION_TEST:
        evaluation_set = [0] * 2**m
        for i in range(len(evaluation_set)):
            evaluation_set[i] = an_element_in_basis(W, i)
    
    for iter in range(N_tests):
        g_coeffs = [FF.random_element() for i in range(2**m)]
        g_coeffs_copy = copy.deepcopy(g_coeffs)

        if DIRECT_EVAUATION_TEST:
            print("Entering Direct Evaluation")
            start = time()
            evaluated_polynomial = evaluate_polynomial(g_coeffs, evaluation_set)
            direct_eval_time[iter] = time() - start

        print("Entering FFT excluding pre-computation")
        start = time()
        fft_no_precmp(g_coeffs, m, nz_hdt_S, table)        
        cantors_fft_no_precmp_time[iter] = time() - start

        if(DIRECT_EVAUATION_TEST and g_coeffs != evaluated_polynomial):
            print("Error: test failed for \"exclude pre-computation\"")
            exit()
        
        print("Entering FFT including pre-computation")
        start = time()
        fft(g_coeffs_copy, m, a, ext_degree)
        cantors_fft_with_precmp_time[iter] = time() - start

        if(g_coeffs_copy != g_coeffs):
            print("Error: test failed for \"includes pre-computation\"")
            exit()

    print("All tests passed")
    if DIRECT_EVAUATION_TEST:
        print("Average direct evaluation time:", sum(direct_eval_time)/N_tests, 's')
    print("Average Cantor's FFT time (excludes pre-computation):", sum(cantors_fft_no_precmp_time)/N_tests, 's')
    print("Average Cantor's FFT time (includes pre-computation):", sum(cantors_fft_with_precmp_time)/N_tests, 's')

    return (sum(direct_eval_time)/N_tests, sum(cantors_fft_no_precmp_time)/N_tests, sum(cantors_fft_with_precmp_time)/N_tests)


def generate_a_map(a, dim):
    for i in range(2**dim):
        print(f"a^{i} = {a**i}")

if __name__ == "__main__":
    # F.<xx> = QQ[]
    # FF.<a> = GF(2**4, modulus= xx**4 + xx + 1)

    F.<x> = GF(2)[]
    ext_degree = 256
    irreducible_poly = F.irreducible_element(ext_degree)
    FF.<a> = GF(2**ext_degree, modulus=irreducible_poly)

    # generate_a_map(a, 4)
    test_fft(a, FF, ext_degree)
    # fast_initial_basis_computation(a, 13)
