from time import time
import copy


load('fft.sage')
load('../utils/utils.sage')

def test_fft(a, FF):
    N_tests = 10
    direct_eval_time = [0] * N_tests
    cantors_fft_no_precmp_time = [0] * N_tests
    cantors_fft_with_precmp_time = [0] * N_tests

    m = 3  
    W, S, nz_hdt_S, table = fft_precmp(a, m)
    # print("W:", W)
    evaluation_set = [0] * 2**m
    for i in range(len(evaluation_set)):
        evaluation_set[i] = an_element_in_basis(W, i)
    
    for iter in range(N_tests):
        g_coeffs = [FF.random_element() for i in range(2**m)]
        g_coeffs_copy = copy.deepcopy(g_coeffs)
        start = time()
        evaluated_polynomial = evaluate_polynomial(g_coeffs, evaluation_set)
        direct_eval_time[iter] = time() - start

        start = time()
        fft_no_precmp(g_coeffs, m, W, S, nz_hdt_S, table)        
        cantors_fft_no_precmp_time[iter] = time() - start

        if(g_coeffs != evaluated_polynomial):
            print("Error: test failed for \"no procomputation\"")
            exit()
            
        start = time()
        fft(g_coeffs_copy, m, a)
        cantors_fft_with_precmp_time[iter] = time() - start

        if(g_coeffs_copy != evaluated_polynomial):
            print("Error: test failed for \"with procomputation\"")
            exit()

    print("All tests passed")
    print("Average direct evaluation time:", sum(direct_eval_time)/N_tests, 's')
    print("Average Cantor's FFT time (no procomputation):", sum(cantors_fft_no_precmp_time)/N_tests, 's')
    print("Average Cantor's FFT time (with procomputation):", sum(cantors_fft_with_precmp_time)/N_tests, 's')

    return (sum(direct_eval_time)/N_tests, sum(cantors_fft_no_precmp_time)/N_tests, sum(cantors_fft_with_precmp_time)/N_tests)


def generate_a_map(a, dim):
    for i in range(2**dim):
        print(f"a^{i} = {a**i}")

if __name__ == "__main__":
    # F.<xx> = QQ[]
    # FF.<a> = GF(2**4, modulus= xx**4 + xx + 1)

    F.<x> = GF(2)[]
    irreducible_poly = F.irreducible_element(16)
    FF.<a> = GF(2**16, modulus=irreducible_poly)

    # generate_a_map(a, 4)
    test_fft(a, FF)