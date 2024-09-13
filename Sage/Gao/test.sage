from time import time
import copy
load('taylor.sage')
load('fft.sage')


def test_taylor(a):
    f_coeffs = [1,a,a^2,a^3,a^4,a^5,a^6]
    print("Input: f_coeffs = ", f_coeffs)
    output = taylor_expansion(f_coeffs)
    print(f"output = {output}")

def test_fft(a, FF):
    N_tests = 10
    DIRECT_EVAUATION_TEST = True
    direct_eval_time = [0] * N_tests
    gaos_fft_with_precmp_time = [0] * N_tests
    gaos_fft_no_precmp_time_lvl1 = [0] * N_tests
    gaos_fft_no_precmp_time_lvl2 = [0] * N_tests

    m = 11
    B = [a**i for i in range(m)]
    if(DIRECT_EVAUATION_TEST):
        evaluation_set = span_basis(B)

    B, G_set, D = fft_precmp_l2(m, B)
    _, G, _ = fft_precmp_l1(m, B)
    for iter in range(N_tests):
        g_coeffs = [FF.random_element() for i in range(2**m)]
        # g_coeffs = [1, a, a**2, 1, a^5, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1]
        # print("input f(x):", polynomial_to_string(g_coeffs))
        g_coeffs_copy_1 = copy.deepcopy(g_coeffs)
        g_coeffs_copy_2 = copy.deepcopy(g_coeffs)


        if(DIRECT_EVAUATION_TEST):
            start = time()
            evaluated_polynomial = evaluate_polynomial(g_coeffs, evaluation_set)
            direct_eval_time[iter] = time() - start

        print(iter, "Entering FFT excluding pre-computation (level2)")
        start = time()
        fft_no_precmp_lvl2(g_coeffs_copy_1, m, B, G_set, D)
        gaos_fft_no_precmp_time_lvl2[iter] = time() - start

        print(iter, "Entering FFT excluding pre-computation (level1)")
        start = time()
        fft_no_precmp_lvl1(g_coeffs_copy_2, m, B, G, D)
        gaos_fft_no_precmp_time_lvl1[iter] = time() - start

        print(iter, "Entering Full FFT including pre-computation")
        start = time()
        fft(g_coeffs, m, B)
        gaos_fft_with_precmp_time[iter] = time() - start

        if(DIRECT_EVAUATION_TEST and g_coeffs != evaluated_polynomial):
            print("Error: test failed")
            exit()

    print("All tests passed")
    if(DIRECT_EVAUATION_TEST):
        print("Average direct evaluation time:", sum(direct_eval_time)/N_tests, 's')
    print("Average Gao's FFT time (excludes pre-computation lvl2):", sum(gaos_fft_no_precmp_time_lvl2)/N_tests, 's')
    print("Average Gao's FFT time (excludes pre-computation lvl1):", sum(gaos_fft_no_precmp_time_lvl1)/N_tests, 's')
    print("Average Gao's FFT time (Full: includes pre-computation):", sum(gaos_fft_with_precmp_time)/N_tests, 's')


def generate_a_map(a, dim):
    for i in range(2**dim):
        print(f"a^{i} = {a**i}")

if __name__ == "__main__":
    # F.<xx> = QQ[]
    # FF.<a> = GF(2**8, modulus= xx**8 + xx**4 + xx**3 + xx**2 +1)
    # FF.<a> = GF(2**4, modulus= xx**4 + xx + 1)

    F.<x> = GF(2)[]
    ext_degree = 256
    irreducible_poly = F.irreducible_element(ext_degree)
    FF.<a> = GF(2**ext_degree, modulus=irreducible_poly)

    # generate_a_map(a, dim=4)

    # test_taylor(a)
    test_fft(a, FF)