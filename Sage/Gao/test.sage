from time import time
import copy
load('taylor.sage')
load('fft.sage')


def test_taylor(a):
    f_coeffs = [1,a,a^2,a^3,a^4,a^5,a^6]
    print("Input: f_coeffs = ", f_coeffs)
    output = taylor_expansion(f_coeffs)
    print(f"output = {output}")

def test_fft_cantor_mode(a, FF, ext_degree):
    mode = CANTOR_BASIS
    m = 3
    B = fast_initial_basis_computation(a, m, ext_degree)
    B.reverse()
    evaluation_set = span_basis(B)
    g_coeffs = [FF.random_element() for i in range(2**m)]
    evaluated_polynomial = evaluate_polynomial(g_coeffs, evaluation_set)
    gm_fft(g_coeffs, m, B, mode)

    if(g_coeffs != evaluated_polynomial):
            print("Error: test failed")
            exit()


def test_fft(a, FF, ext_degree):
    N_tests = 10
    DIRECT_EVAUATION_TEST = True
    CANTOR_BASIS_TEST = False
    direct_eval_time = [0] * N_tests
    gm_fft_with_precmp_time = [0] * N_tests
    gm_fft_no_precmp_time_lvl1 = [0] * N_tests
    gm_fft_no_precmp_time_lvl2 = [0] * N_tests

    m = 11
    if CANTOR_BASIS_TEST:
        B = fast_initial_basis_computation(a, m, ext_degree)
        B.reverse()
    else:
        B = [a**i for i in range(m)]
    
    if(DIRECT_EVAUATION_TEST):
        evaluation_set = span_basis(B)

    G_set, D = gm_fft_precmp_l2(m, B)
    G, _ = gm_fft_precmp_l1(m, B)
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
        gm_fft_no_precmp_lvl2(g_coeffs_copy_1, m, B, G_set, D)
        gm_fft_no_precmp_time_lvl2[iter] = time() - start

        print(iter, "Entering FFT excluding pre-computation (level1)")
        start = time()
        gm_fft_no_precmp_lvl1(g_coeffs_copy_2, m, B, G, D)
        gm_fft_no_precmp_time_lvl1[iter] = time() - start

        print(iter, "Entering Full FFT including pre-computation")
        start = time()
        gm_fft(g_coeffs, m, B)
        gm_fft_with_precmp_time[iter] = time() - start

        if(DIRECT_EVAUATION_TEST and g_coeffs != evaluated_polynomial):
            print("Error: test failed")
            exit()

    print("All tests passed")
    if(DIRECT_EVAUATION_TEST):
        print("Average direct evaluation time:", sum(direct_eval_time)/N_tests, 's')
    print("Average Gao's FFT time (excludes pre-computation lvl2):", sum(gm_fft_no_precmp_time_lvl2)/N_tests, 's')
    print("Average Gao's FFT time (excludes pre-computation lvl1):", sum(gm_fft_no_precmp_time_lvl1)/N_tests, 's')
    print("Average Gao's FFT time (Full: includes pre-computation):", sum(gm_fft_with_precmp_time)/N_tests, 's')


def generate_a_map(a, dim):
    for i in range(2**dim):
        print(f"a^{i} = {a**i}")

if __name__ == "__main__":
    # F.<xx> = QQ[]
    # FF.<a> = GF(2**8, modulus= xx**8 + xx**4 + xx**3 + xx**2 +1)
    # FF.<a> = GF(2**4, modulus= xx**4 + xx + 1)

    F.<x> = GF(2)[]
    ext_degree = 264
    irreducible_poly = F.irreducible_element(ext_degree)
    FF.<a> = GF(2**ext_degree, modulus=irreducible_poly)

    # generate_a_map(a, dim=4)

    # test_taylor(a)
    test_fft(a, FF, ext_degree)
    # test_fft_cantor_mode(a, FF, ext_degree)