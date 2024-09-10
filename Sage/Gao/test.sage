from random import randint
from time import time
load('taylor.sage')
load('fft.sage')


def test_taylor(a):
    f_coeffs = [1,a,a^2,a^3,a^4,a^5,a^6]
    print("Input: f_coeffs = ", f_coeffs)
    output = taylor_expansion(f_coeffs)
    print(f"output = {output}")

def test_fft(a, FF):
    N_tests = 10
    DIRECT_EVAUATION_TEST = False
    direct_eval_time = [0] * N_tests
    gaos_fft_time = [0] * N_tests

    m = 14
    B = [a**i for i in range(m)]
    if(DIRECT_EVAUATION_TEST):
        evaluation_set = [0] * 2**m
        for i in range(len(evaluation_set)):
            evaluation_set[i] = an_element_in_basis(B, i)

    for iter in range(N_tests):
        g_coeffs = [FF.random_element() for i in range(2**m)]

        if(DIRECT_EVAUATION_TEST):
            start = time()
            evaluated_polynomial = evaluate_polynomial(g_coeffs, evaluation_set)
            direct_eval_time[iter] = time() - start

        print(iter, "Entering FFT")
        start = time()
        fft(g_coeffs, m, B)
        gaos_fft_time[iter] = time() - start

        if(DIRECT_EVAUATION_TEST and g_coeffs != evaluated_polynomial):
            print("Error: test failed")
            exit()

    print("All tests passed")
    if(DIRECT_EVAUATION_TEST):
        print("Average direct evaluation time:", sum(direct_eval_time)/N_tests, 's')
    print("Average Gao's FFT time:", sum(gaos_fft_time)/N_tests, 's')


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