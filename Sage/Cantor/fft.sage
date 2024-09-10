import math

def fft(g_coeffs, m, a, ext_degree):
    """
    Computes the Fast Fourier Transform (FFT) of the polynomial g(x) using Cantor's algorithm.

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

    * The FFT is computed in place, so g_coeffs is modified to contain the evaluations of g(x) over Cantor's basis (W).
      Initially, g_coeffs represents the coefficients of the polynomial g(x), but by the end of the algorithm, it will 
      hold the values of g(x) evaluated at specific points from    
    """
    g_coeffs += [0]*(2**(m)-len(g_coeffs)) 
    _, nz_hdt_S, table = fft_precmp(a, m, ext_degree)
    fft_no_precmp(g_coeffs, m, nz_hdt_S, table)

def fft_no_precmp(g_coeffs, m, nz_hdt_S, table):
    """
    Performs the core computations of the Fast Fourier Transform (FFT) algorithm, excluding the precomputation of 
    S(x) polynomials and the shift table. This function assumes that these precomputations have already been done.

    Parameters:
    g_coeffs (list): A list of coefficients representing the polynomial g(x). This list will be modified in place 
                     during the FFT and will ultimately contain the evaluations of g(x) over the points in Cantor's basis.
    m (int): deg(g(x)) < 2^m. This corresponds to the number of stages in the FFT.
    nz_hdt_S (list[list]): A list of lists, where each inner list contains the indices of the non-zero coefficients 
                           in the polynomial s_r(x) for each round r. These indices are ordered from the highest 
                           degree term (hdt) to the lowest degree term.
    table (list[list]): A table of precomputed shifts used for tail module computations.

    * The result of the FFT is stored directly in g_coeffs, which is updated in place.
    """
    g_coeffs += [0]*(2**(m)-len(g_coeffs))
    input_size = len(g_coeffs)
    n_modules = 1

    head_module(g_coeffs, nz_hdt_S[0], input_size)
    for r in range(1,m):
        input_size >>= 1
        n_modules  <<= 1
        head_module(g_coeffs, nz_hdt_S[r], input_size)
        offset = 0
        for i in range(0, n_modules-1):
            offset += input_size
            tail_module(g_coeffs, nz_hdt_S[r], input_size, offset, table[r][i])

def fft_precmp(a, m, ext_degree):
    """
    Performs the precomputations required for Cantor's FFT algorithm. These precomputations include
    (1) constructing Cantor's basis, (2) computing the polynomials s_r(x) for each FFT round and determining 
    the corresponding indices of their non-zero coefficients, and (3) generating a table of precomputed shifts 
    for use in the tail module computations.

    Parameters: 
    a (element of GF(2^ext_degree)): The primitive element of the finite field GF(2^ext_degree). It is the root of the 
                                     irreducible polynomial that defines the field.
    m (int): deg(g(x)) < 2^m. This determines the size of the FFT and the number of rounds.
    ext_degree (int): The extension degree of the field GF(2), such that the field is GF(2^ext_degree). In Cantor's algorithm, 
                      ext_degree must be a power of two to compute Canto's basis.

    Returns:
    W (list): Cantor's basis, which defines the set of points where the polynomial will be evaluated. 
    S (list[list]): A list of polynomials s_r(x) for each FFT round r. Each polynomial s_r(x) is represented as a list of 
                    its coefficients. These polynomials are used in the modular reductions at each stage of the FFT.
    nz_hdt_S (list[list]): A list of lists where each inner list contains the indices of the non-zero coefficients 
                           in the polynomial s_r(x) for each round r, ordered from the highest degree term (hdt) to the lowest.
    table (list[list]): A precomputed table of shifts used in tail module computations during each FFT round.

    """
    W = fast_initial_basis_computation(a, m, ext_degree)
    nz_hdt_S = S_function_computation(m)
    table = fast_S_shifts_table_generator(W)
    return W, nz_hdt_S, table

def initial_basis_computation(a, m):
    """
    Computes Cantor's basis non-optimally. The algorithm initializes W[0] = 1 and derives each subsequent W[i], such that 
    W[i]^2 + W[i] = W[i-1] for i >= 1. It performs a search in the field GF(2^ext_degree) to find each W[i]. 
    This is a brute-force approach and not an efficient method.
    
    Parameters: 
    a (element of GF(2^ext_degree)): The primitive element of the finite field GF(2^ext_degree). It is the root of the 
                                     irreducible polynomial that defines the field.
    m (int): deg(g(x)) < 2^m. This determines the size of the FFT and the number of rounds.
    
    Returns: 
    W (list): Cantor's basis, which defines the set of points where the polynomial will be evaluated. 
    """
    W = [1] * m
    for i in range(1, m):
        b = a
        while b**2 + b != W[i-1]:
            b *= a
        W[i] = b
    return W

def fast_initial_basis_computation(a, m, ext_degree):
    """
    Computes Cantor's basis optimally using the trace function over GF(2^ext_degree). The algorithm first searches for W[m-1] 
    such that tr(W[m-1]) = 1, where tr() is the trace function. Once W[m-1] is found, the algorithm recursively determines 
    W[i-1] for each i, such that W[i]^2 + W[i] = W[i-1].
     
    Parameters: 
    a (element of GF(2^ext_degree)): The primitive element of the finite field GF(2^ext_degree). It is the root of the 
                                     irreducible polynomial that defines the field.
    m (int): deg(g(x)) < 2^m. This determines the size of the FFT and the number of rounds.
    ext_degree (int): The extension degree of the field GF(2), such that the field is GF(2^ext_degree). In Cantor's algorithm, 
                      ext_degree must be a power of two to compute Canto's basis.

    Returns: 
    W (list): Cantor's basis, which defines the set of points where the polynomial will be evaluated. 
    """
    W = [1] * ext_degree
    b_m = a
    while(b_m.trace() != 1):
        b_m *= a
    W[-1] = b_m

    for i in reversed(range(1, ext_degree-1)):
        W[i] = W[i+1]**2 + W[i+1]
    return W[:m]

def S_function_computation(m):
    """
    Computes s_r(x) polynomials for each round. 
     
    Parameters: 
    a (element of GF(2^ext_degree)): The primitive element of the finite field GF(2^ext_degree). It is the root of the 
                                     irreducible polynomial that defines the field.
    m (int): deg(g(x)) < 2^m. This determines the size of the FFT and the number of rounds.
    ext_degree (int): The extension degree of the field GF(2), such that the field is GF(2^ext_degree). In Cantor's algorithm, 
                      ext_degree must be a power of two to compute Canto's basis.

    Returns: 
    W (list): Cantor's basis, which defines the set of points where the polynomial will be evaluated. 
    """
    S = [[]] * (m)
    # nz_hdt_S: Non zero coeffs indices in S. Note that index 0 is the coefficient of the highest degree term in s_r(x),
    nz_hdt_S = [[]] * (m)                   
    for r in range(m):
        S[r] = [0] * (2**r + 1) # deg(S_r) = 2^r because 'C(r,r) = 1' 
        nz_hdt_S[r] = []
        for i in range(r+1): # S_r(t) = \sum_{i=0}^{r} C(r,i)*(t^{2^i}). Therefore, the loop should include i = r
            S[r][2**i] = math.comb(r,i) % 2 # Coefficients are ordered from the constant term up to the highest degree term.
            if (S[r][2**i] == 1): nz_hdt_S[r].append(2**r - 2**i) # Non-zero coefficients indices in S ordered from highest degree term (hdt)       
    nz_hdt_S.reverse()
    S.reverse()
    return nz_hdt_S

def fast_S_function_computation(m):
    # Todo: Implement S function computation using the fact that S_{m+1}(x) = S(S_m(x)) = S_m^2(x) + S_m(x) 
    # S = [[]] * (m)
    # nz_hdt_S = [[]] * (m) 
    # S[0] = [0,1] # S_0 = x
    # # S[1] = [0,0,1]  
    # for r in range(1,m):
    #     S[r] = ([0] * (1<<(r) + 1)) # deg(S_r) = 2^r hence it has at most 2^r + 1 coefficients
    #     # print(f"S[{r}] = {S[r]}")
    #     for i in range(1<<(r-1) + 1):
    #         S[r][2*i] = S[r-1][i] 
    #         S[r][i] += S[r-1][i] 
    #         # print(i, f"S[{r}] = {S[r]}")
    #     for i in range((1<<r)+1):
    #         if (S[r][i] == 1): nz_hdt_S[r].append((1<<r) - i)
    # nz_hdt_S.reverse()
    # S.reverse()
    # return nz_hdt_S
    pass
        

def S_shifts_table_generator(S, W):
    table = []
    for i in range(len(W)):
        row = []
        for j in range((2**(i+1))):
            if i == 0:
                shift = 0
            else:
                shift = table[-1][j>>1]
            if j % 2:
                shift += W[~i]         #~i points to the i-th last element in the list:  ~ is the bitwise NOT where ~i = -i-1. bitwise NOT is used because its time complexity is O(1)
            row += [shift]
            # print(f"eval_at_x(S[{i}], {shift}) = ", eval_at_x(S[i], shift))
        table += [row]    
    for i in range(len(table)):
        for j in range(len(table[i])):
            table[i][j] = eval_at_x(S[i], table[i][j])
        table[i] = table[i][::2]
        table[i] = table[i][1:]

    return table

def fast_S_shifts_table_generator(W):
    table = []
    for r in range(0, len(W)): # r is the row number. The first row is skipped as it only has a head module.
        row = [0] * ((1<<(r)) - 1)
        for module in range(1<<(r)): # 1<<(r) = 2^r which is the number of modlues in the row r. The first module is skipped as it is a head module.
            for i in range(r):
                if module & (1<<i): 
                    row[module-1] += W[i+1]
        table += [row]
    return table


def eval_at_x(s, x):
    result = 0
    for i in range(len(s)):
        result += s[i] * (x ** i)
    return result

def divide(coeffs, nz_hdt_S, input_size, offset):
    """
    Computes the quotient q(x) = g(x) / s(x), where the dividend polynomial g(x) is represented 
    by the coefficients in coeffs[offset:offset+input_size], and the divisor polynomial s(x) 
    is represented by the non-zero indices in nz_hdt_S. s(x) is in GF(2).

    Parameters:
    coeffs (list): List of coefficients representing the polynomials. This is modified in place 
                   to store the remainder after division.
    nz_hdt_S (list): List of indices representing the non-zero terms of the divisor polynomial s(x), 
                     where the coefficients of s(x) are in GF(2).
    input_size (int): Number of coefficients in the dividend polynomial g(x).
    offset (int): Starting index in 'coeffs' where the dividend polynomial g(x) begins.

    Returns:
    q (list): Coefficients of the quotient polynomial q(x).
    """
    # Initialize the quotient polynomial q(x), whose degree is at most half the degree of the dividend g(x).
    q = [0] * (input_size>>1)       # The degree of q(x) is derived as deg(g(x)) - deg(s(x)).
                                    # For input size 2^(r+1), we have deg(q) <= 2^r - 1 = input_size // 2 - 1.        
    offset += len(q)

    # Process each term in the quotient by iterating over the coefficients in reverse order.
    for q_ind, i in enumerate(reversed(range(offset,offset+len(q)))):   # The reason to use 'reversed()' is that coefficients are ordered from 
                                                                        # the constant term up to the highest degree term. 
                                                                        # Therefore, the last element is the highest degree coefficient.
        q[~q_ind] = coeffs[i]       # Using ~q_ind accesses elements from the end.

        for nz in nz_hdt_S[:-1]:    # Exclude the highest degree term of s(x) (the last element of nz_hdt_S),
                                    # because it cancels out with coeffs[i].
            coeffs[i - nz] += coeffs[i]
        coeffs[i] = 0 
        q_ind += 1
    return (q)

def head_module(coeffs, nz_hdt_S, input_size):
    q = divide(coeffs, nz_hdt_S, input_size, offset=0)
    for i in range(len(q)):
        coeffs[i+len(q)] = coeffs[i] + q[i]

def tail_module(coeffs, nz_hdt_S, input_size, offset, s_shift1):
    q = divide(coeffs, nz_hdt_S, input_size, offset)
    for i in range(len(q)):
        coeffs[offset+i] = coeffs[offset+i] + q[i] * s_shift1
        coeffs[offset+i+len(q)] = coeffs[offset+i] + q[i] 