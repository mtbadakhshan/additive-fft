load('../utils/utils.sage')

def taylor_phase(g_coeffs, mu, s):
    """
    In-place Taylor expansion phase T(mu, s).

    The array g_coeffs is partitioned into independent subproblems of length
    2^mu stored with stride 2^s: the subproblem bases are all indices b whose
    bits [s, s+mu) are zero. Each subproblem polynomial is Taylor-expanded
    w.r.t. S^{mu1}(x) = x^(2^mu1) + x, where mu1 is the largest power of two
    strictly less than mu (valid since the basis is a Cantor basis).

    All 2^s interleaved subproblems inside a block of size 2^(s+mu) share the
    same additions, so each block is processed with plain strided XOR loops.
    """
    mu1 = 1
    while mu1 * 2 < mu:
        mu1 *= 2
    mu2 = mu - mu1

    sigma = 1 << s
    block_size = sigma << mu
    shift = ((1 << mu1) - 1) * sigma  # (2^mu1 - 1) * 2^s

    for b0 in range(0, len(g_coeffs), block_size):
        for j in range(1 << mu2):
            lo = b0 + ((j + 1) << mu1) * sigma
            for p in range(b0 + block_size - 1, lo - 1, -1):
                g_coeffs[p - shift] += g_coeffs[p]

def butterfly_stage(g_coeffs, s, theta, cantor_basis, m):
    """
    In-place butterfly stage B(2^s): for every pair (p, p + 2^s) with
    bit_s(p) = 0, evaluate the length-2 polynomial [g[p], g[p + 2^s]] at
    {lam, lam + 1}, where

        lam = S^s(theta) + sum_{j > s} bit_j(p) * beta[j - s]

    (equivalently S^s(omega_p) with omega_p the p-th evaluation point),
    using S(t) = t^2 + t and the Cantor property S(beta[j]) = beta[j-1].
    """
    sigma = 1 << s
    theta_s = theta
    for _ in range(s):
        theta_s = theta_s**2 + theta_s  # S^s(theta)

    for p in range(len(g_coeffs)):
        if p & sigma:
            continue
        lam = theta_s
        for j in range(s + 1, m):
            if (p >> j) & 1:
                lam += cantor_basis[j - s]
        g_coeffs[p] += lam * g_coeffs[p + sigma]
        g_coeffs[p + sigma] += g_coeffs[p]

def fft(g_coeffs, m, a, ext_degree, affine_shift=0):
    """
    Computes the Fast Fourier Transform (FFT) of the polynomial g(x) using the dyadic AFFT algorithm.

    Parameters:
    g_coeffs (list): A list of coefficients representing the polynomial g(x), ordered from the constant term 
                     to the highest degree term. After the FFT, this list will hold the evaluations of g(x) 
                     at points in the Cantor's basis.    
    m (int): An upper bound on the degree of g(x), ensuring that deg(g) < 2^m. This controls the size of the FFT.
    a (element of GF(2^ext_degree)): The primitive element of the finite field GF(2^ext_degree), which is used 
                                     to define the FFT points. It serves as the root of the irreducible polynomial 
                                     defining the finite field.
    ext_degree (int): The extension degree of GF(2), such that the finite field is GF(2^ext_degree). This controls 
                      the field in which the FFT computations are performed. In Cantor's algorithm, ext_degree 
                      must be a power of two.
    affine_shift (element of GF(2^ext_degree)): Optional shift theta; the polynomial is evaluated over the affine
                                                space theta + W_m instead of W_m.

    * The FFT is computed in place, so g_coeffs is modified to contain the evaluations of g(x) over Cantor's basis (W).
      Initially, g_coeffs represents the coefficients of the polynomial g(x), but by the end of the algorithm, it will 
      hold the values of g(x) evaluated at specific points from theta + W_m, in binary counting order of the basis.

    The recursive dyadic AFFT
        FFT(f, mu) = Taylor(f) -> column FFTs of length 2^mu2 -> row FFTs of length 2^mu1
    is flattened here into a schedule of Taylor phases and butterfly stages:
        Sched(mu, s) = T(mu, s) + Sched(mu2, s + mu1) + Sched(mu1, s),   Sched(1, s) = B(2^s).
    Each butterfly stride 2^s, s = m-1 .. 0, occurs exactly once as a full radix-2 stage.
    The schedule is generated iteratively with a small stack of (mu, s) tasks.
    """
    g_coeffs += [0]*(2**(m)-len(g_coeffs)) # Pad with zeros to ensure the length is 2^m

    theta = affine_shift
    cantor_basis = fast_initial_basis_computation(a, m, ext_degree)

    stack = [(m, 0)]
    print("stack:", stack)
    while stack:
        mu, s = stack.pop()
        # Follow the column-FFT chain; row tasks are stacked for later.
        while mu > 1:
            print("start: stack:", stack, "mu:", mu, "s:", s)
            mu1 = 1 # mu1 must be the largest power of 2 strictly less than mu
            while mu1 * 2 < mu:
                mu1 *= 2
            mu2 = mu - mu1
            print("taylor_phase: mu:", mu, "s:", s, "mu1:", mu1, "mu2:", mu2)
            taylor_phase(g_coeffs, mu, s)
            stack.append((mu1, s))
            mu, s = mu2, s + mu1
            print("----- stack:", stack, "mu:", mu, "s:", s)
        butterfly_stage(g_coeffs, s, theta, cantor_basis, m)
        print("butterfly_stage: mu:", mu, "s:", s, "m:", m)

    return g_coeffs
