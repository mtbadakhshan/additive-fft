load('taylor.sage')
load('../utils/utils.sage')

def fft(g_coeffs, m, B):
    """
    Computes the Fast Fourier Transform (FFT) of the polynomial g(x) using Gao's algorithm.

    Parameters:
    g_coeffs (list): A list of coefficients representing the polynomial g(x), ordered from the constant term up to the highest degree term.
    m (int): The upper bound on the degree of g(x), such that deg(g) < 2^m.
    B (list): The basis of the evaluation set.

    Returns:
    list: The FFT of the polynomial g(x) evaluated at the points in B.
    """
    g_coeffs += [0]*(2**(m)-len(g_coeffs))
    B, G, D = fft_precmp_l1(m, B)

    fft_f(g_coeffs, m, [B] + D)
    fft_r_l1(g_coeffs, m, G)

def fft_no_precmp_lvl1(g_coeffs, m, B, G, D):
    fft_f(g_coeffs, m, [B] + D)
    fft_r_l1(g_coeffs, m, G)

def fft_no_precmp_lvl2(g_coeffs, m, B, G_set, D):
    fft_f(g_coeffs, m, [B] + D)
    fft_r_l2(g_coeffs, m, G_set)


def fft_precmp_l1(m, B):
    G = [[]] * (m-1); D = [[]] * (m-1)
    G[0], D[0] = G_D_computation(B)
    for i in range(1,m-1):
        G[i], D[i] = G_D_computation(D[i-1])
    return B, G, D

def fft_precmp_l2(m, B):
    B, G, D = fft_precmp_l1(m, B)
    for i in range(0, m-1):
        G[i] = span_basis(G[i])
    return B, G, D


def G_D_computation(B):
    """
    Computes the G and D bases from the basis B according to the method described in the referenced paper.

    Parameters:
    B (list): The input basis from which G and D are derived. It is assumed that B is non-empty.

    Returns:
    tuple: A tuple containing two lists:
        - G (list): The computed G basis, where each element is derived from B.
        - D (list): The computed D basis, where each element is calculated as the square of G[i]^2 - G[i].
    """
    G = [0] * (len(B)-1)
    D = [0] * len(G)
    for i in range(len(G)):
        G[i] = B[i] * (B[-1] ** (-1))
        D[i] = (G[i] ** 2) - G[i]
    return G, D

def scale_polynomial_by_beta(coeffs, input_size, beta):
    """
    Scales the polynomial coefficients by powers of beta, effectively computing the coefficients of g(beta * x).

    Parameters:
    coeffs (list): The list of coefficients of the input polynomial g(x), ordered from the constant term up to the highest degree term.
    input_size (int): The size of each segment in the coeffs list to which beta will be applied.
    beta (int): The coefficient by which to scale the x terms.

    Returns:
    coeffs (list): The updated list of coefficients representing g(beta * x).
    """
    for i in range(0,len(coeffs),input_size):
        beta_p = 1
        for j in range(input_size):
            coeffs[i+j] = beta_p * coeffs[i+j]
            beta_p *= beta
    # return coeffs


def fft_f(coeffs, m, B):
    """
    Forward process in the FFT algorithm.

    Parameters:
    coeffs (list): The list of coefficients of the input polynomial, ordered from the constant term up to the highest degree term.
    m (int): The upper bound of the degree of the input polynomial, such that the degree is less than 2^m.
    B (list of lists): The list of precomputed bases for all rounds. The basis for the first round corresponds to the evaluation set,
              and for later rounds, it is equal to the D basis computed in the G_D_computation function.

    Returns:
    list: The result of the forward FFT process applied to the input polynomial coefficients.
    """
    input_size = len(coeffs)
    for r in range(m-1):
        scale_polynomial_by_beta(coeffs, input_size, B[r][-1])  # B[i][-1] := beta_m
        offset = 0
        for b in range(1<<r):
            coeffs[offset:offset+input_size] = taylor_expansion(coeffs[offset:offset+input_size], input_size)
            offset+=input_size
        input_size >>= 1
    assert(input_size == 2)
    assert(len(B[-1]) == 1)
    for i in range(0, len(coeffs), 2):
        coeffs[i:i+2] = evaluate_polynomial(coeffs[i:i+2], [0] + B[-1])

def fft_f_2(coeffs, m, B):
    """
    Forward process in the FFT algorithm.

    Parameters:
    coeffs (list): The list of coefficients of the input polynomial, ordered from the constant term up to the highest degree term.
    m (int): The upper bound of the degree of the input polynomial, such that the degree is less than 2^m.
    B (list of lists): The list of precomputed bases for all rounds. The basis for the first round corresponds to the evaluation set,
              and for later rounds, it is equal to the D basis computed in the G_D_computation function.

    Returns:
    list: The result of the forward FFT process applied to the input polynomial coefficients.
    """
    input_size = len(coeffs)

    for r in range(m-1):
        scale_polynomial_by_beta(coeffs, input_size, B[r][-1])  # B[i][-1] := beta_m
        offset = 0
        for b in range(1<<r):
            taylor_expansion_no_post_cmp(coeffs, input_size, offset)
            offset+=input_size
        g_0_g_1_extraction(coeffs, input_size, 1<<r)


        input_size >>= 1
    assert(input_size == 2)
    assert(len(B[-1]) == 1)
    for i in range(0, len(coeffs), 2):
        coeffs[i:i+2] = evaluate_polynomial(coeffs[i:i+2], [0] + B[-1])




def fft_r_l1(coeffs, m, G):
    """
    Reverse process in the FFT algorithm.

    Parameters:
    coeffs (list): The list of coefficients after the final round of forward computation. 
                   It consists of 2^(m-1) concatenated tuples of (g_0(x), g_1(x)), where each has a degree of 1.
    m (int): The upper bound of the degree of the input polynomial, such that the degree is less than 2^m.
    G (list of lists): The list of precomputed G bases for all rounds, obtained from the G_D_computation function.

    Returns:
    list: The coefficients of the polynomial after applying the reverse FFT process.
    """
    input_size = 2
    for r in reversed(range(1,m)):
        offset=0
        for b in range(0,2**r,2):
            for i in range(input_size):
                coeffs[offset + i] += an_element_in_basis(G[r-1], index=i) * coeffs[offset + input_size + i]
                coeffs[offset + input_size + i] += coeffs[offset + i]
            offset += input_size<<1      
        input_size <<= 1  


def fft_r_l2(coeffs, m, G_set):
    """
    Reverse process in the FFT algorithm.

    Parameters:
    coeffs (list): The list of coefficients after the final round of forward computation. 
                   It consists of 2^(m-1) concatenated tuples of (g_0(x), g_1(x)), where each has a degree of 1.
    m (int): The upper bound of the degree of the input polynomial, such that the degree is less than 2^m.
    G (list of lists): The list of spanned G bases for all rounds.

    Returns:
    list: The coefficients of the polynomial after applying the reverse FFT process.
    """
    input_size = 2
    for r in reversed(range(1,m)):
        offset=0
        for b in range(0,2**r,2):
            for i in range(input_size):
                coeffs[offset + i] += G_set[r-1][i] * coeffs[offset + input_size + i]
                coeffs[offset + input_size + i] += coeffs[offset + i]
            offset += input_size<<1      
        input_size <<= 1  