from time import time
import copy

load('fft.sage')
load('../utils/utils.sage')

def test_fft_correctness(a, FF, ext_degree):
    """Compare fft() against naive evaluation for several sizes, with and without affine shift."""
    all_ok = True
    for m in range(11, 13):
        cantor_basis = fast_initial_basis_computation(a, m, ext_degree)
        for affine_shift in [FF.zero(), FF.random_element()]:
            g_coeffs = [FF.random_element() for i in range(2**m)]
            evals = fft(copy.deepcopy(g_coeffs), m, a, ext_degree, affine_shift=affine_shift)
            naive = evaluate_polynomial(g_coeffs, span_basis(cantor_basis, affine_shift))
            ok = (evals == naive)
            all_ok = all_ok and ok
            print(f"m={m}, shift={affine_shift}: {'OK' if ok else 'MISMATCH'}")
            if not ok:
                print("  fft  :", evals)
                print("  naive:", naive)
    print("==================================")
    print("Correctness tests passed!" if all_ok else "SOME CORRECTNESS TESTS FAILED")
    return all_ok

def test_fft(a, FF, ext_degree):
    N_tests = 1
    DIRECT_EVAUATION_TEST = False

    direct_eval_time = [0] * N_tests
    dyadic_fft_time = [0] * N_tests

    m = 4
    print("Testing FFT with m =", m, "and extension degree =", ext_degree)

    if DIRECT_EVAUATION_TEST:
        evaluation_set = span_basis(W, affine_shift)

    for iter in range(N_tests):
        g_coeffs = [FF.random_element() for i in range(2**m)]
        g_coeffs_copy = copy.deepcopy(g_coeffs)

        if DIRECT_EVAUATION_TEST:
            print(iter, "Entering Direct Evaluation")
            start = time()
            evaluated_polynomial = evaluate_polynomial(g_coeffs, evaluation_set)
            direct_eval_time[iter] = time() - start

        print(iter, "Entering dyadic AFFT")
        start = time()
        fft(g_coeffs_copy, m, a, ext_degree)        
        dyadic_fft_time[iter] = time() - start

        if(DIRECT_EVAUATION_TEST and g_coeffs_copy != evaluated_polynomial):
            print("Error: test failed for \"dyadic AFFT\"")
            exit()

if __name__ == "__main__" or True:
    F.<x> = GF(2)[]
    ext_degree = 32
    irreducible_poly = F.irreducible_element(ext_degree)
    FF.<a> = GF(2**ext_degree, modulus=irreducible_poly)
    test_fft_correctness(a, FF, ext_degree)
    # test_fft(a, FF, ext_degree)

