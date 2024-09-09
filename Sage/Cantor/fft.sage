import math
# import copy

n_add = 0
n_mult = 0

def initial_basis_computation(a, dim):
    W = [1] * dim
    for i in range(1, dim):
        b = a
        while b**2 + b != W[i-1]:
            b *= a
        W[i] = b
    return W

# def fast_initial_basis_computation(a, dim):
#     for b in range(1, 2^dim):
#         sum = 0
#         for i in range(dim):

    
    W = [1] * dim
    for i in range(1, dim):
        b = a
        while b**2 + b != W[i-1]:
            b *= a
        W[i] = b
    return W


def S_function_computation(dim):
    S = [[]] * (dim)
    # nz_hdt_S: Non zero coeffs indices in S. Note that index 0 is the coefficient of the highest degree term
    nz_hdt_S = [[]] * (dim)                   
    for r in range(dim):
        S[r] = [0] * (2**r + 1) # deg(S_r) = 2^r because 'C(r,r) = 1' 
        nz_hdt_S[r] = []
        for i in range(r+1): # S_r(t) = \sum_{i=0}^{r} C(r,i)*(t^{2^i}). Therefore, the loop should include i = r
            S[r][2**i] = math.comb(r,i) % 2 # Coefficients are ordered from the constant term up to the highest degree term.
            if (math.comb(r,i) % 2 == 1): nz_hdt_S[r].append(2**r - 2**i) # Non-zero coefficients indices in S ordered from highest degree term (hdt)       
    S.reverse()
    nz_hdt_S.reverse()
    return S, nz_hdt_S

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
    return table

def eval_at_x(s, x):
    result = 0
    for i in range(len(s)):
        result += s[i] * (x ** i)
    return result

def divide(coeffs, nz_hdt_S, input_size, offset):
    """
    input_size (int): It is the number of coefficients in coeffs that we see as the dividend.
    offset (int): The index in coeffs where the dividend is started
    """
    global n_add
    # print("coeffs: ", coeffs)
    # print("nz_hdt_S: ", nz_hdt_S)

    q = [0] * (input_size>>1) # In general case,  deg(q) = deg(coeffs) - deg(S_r). Also, deg(coeffs) <= input_size - 1 and the input_size = 2^(r+1).
                                  # Hence, in each round r, deg(coeffs) <= 2^(r+1)-1 and deg(S_r) = 2^r, hence deg(q) <= 2^r -1 = input_size/2 - 1.         
    offset += len(q)
    for q_ind, i in enumerate(reversed(range(offset,offset+len(q)))):
        # Coefficients are ordered from the constant term up to the highest degree term. 
        # Therefore, the last element is the highest degree coefficient. 
        q[~q_ind] = coeffs[i]  
        for nz in nz_hdt_S[:-1]: 
            # print("nz =", nz)
            # We ignor the highest degree of S (last element of nz_hdt_S), because it will be canceled with coeffs[i].
            coeffs[i - nz] += coeffs[i]; n_add += 1
        coeffs[i] = 0 
        q_ind += 1
    return (q)

def head_module(coeffs, nz_hdt_S, input_size):
    q = divide(coeffs, nz_hdt_S, input_size, offset=0)
    for i in range(len(q)):
        coeffs[i+len(q)] = coeffs[i] + q[i]

def tail_module(coeffs, nz_hdt_S, input_size, offset, s_shift1, s_shift2):
    # print("tail module::")
    # print("g_coeffs: ", coeffs)
    q = divide(coeffs, nz_hdt_S, input_size, offset)
    # print("g_coeffs: ", coeffs)
    # print("s_shift1:", s_shift1)
    # print("s_shift2:", s_shift2)
    for i in range(len(q)):
        # print(f"q[{i}]:", q[i])
        # print(f"q[{i}] * s_shift2:", q[i] * s_shift2)
        # print(f"coeffs[{offset+i+len(q)}] + q[{i}] * s_shift2:", q[i] * s_shift2)

        # print(f"q[{i}] * s_shift2:", q[i] * s_shift2)
        # print(f"coeffs[i] + q[i] * s_shift2:", q[i] * s_shift2)
        coeffs[offset+i+len(q)] = coeffs[offset+i] + q[i] * s_shift2
        coeffs[offset+i] = coeffs[offset+i] + q[i] * s_shift1
    
def fft_precmp(a, m):
    W = initial_basis_computation(a, dim=m)
    S, nz_hdt_S = S_function_computation(dim=m)
    table = S_shifts_table_generator(S, W)
    return W, S, nz_hdt_S, table

def fft_no_precmp(g_coeffs, m, W, S, nz_hdt_S, table):
    input_size = len(g_coeffs)
    n_modules = 1
    for r in range(m):
        # print("r: ", r)
        # print("g_coeffs: ", g_coeffs)
        head_module(g_coeffs, nz_hdt_S[r], input_size)
        # print("g_coeffs (after head): ", g_coeffs)
        offset = 0
        for i in range(1, n_modules):
            offset += input_size
            tail_module(g_coeffs, nz_hdt_S[r], input_size, offset, table[r][2*i], table[r][2*i+1])
            # print("g_coeffs (after tail): ", g_coeffs)

        n_modules  <<= 1
        input_size >>= 1

def fft(g_coeffs, m, a):
    # global n_add, n_mult

    g_coeffs += [0]*(2**(m)-len(g_coeffs))
    W, S, nz_hdt_S, table = fft_precmp(a, m)
    fft_no_precmp(g_coeffs, m, W, S, nz_hdt_S, table)

    