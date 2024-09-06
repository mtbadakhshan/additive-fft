import math
import copy

n_add = 0
n_mult = 0

def initial_basis_computation(a, dim):
    B = [1] * dim
    for i in range(1, dim):
        b = a
        while True:
            if b**2 + b == B[i-1]:
                break
            b *= a
        B[i] = b
    return B

def S_function_computation(dim):
    S = [[]] * (dim)
    # ns_hdt_S: Non zero coeffs indices in S. Note that index 0 is the coefficient of the highest degree term
    ns_hdt_S = [[]] * (dim)                   
    for r in range(dim):
        S[r] = [0] * (2**r + 1) ## deg(S_r) = 2^r because 'C(r,r) = 1' 
        ns_hdt_S[r] = []
        for i in range(r+1): # S_r(t) = \sum_{i=0}^{r} C(r,i)*(t^{2^i}). Therefore, the loop should include i = r
            S[r][2**i] = math.comb(r,i) % 2 # Coefficients are ordered from the constant term up to the highest degree term.
            if (math.comb(r,i) % 2 == 1): ns_hdt_S[r].append(2**r - 2**i) #Non-zero coefficients indices in S ordered from highest degree term (hdt) 
            
    return S, ns_hdt_S

def divide(coeffs, ns_hdt_S, input_size, offset):
    """
    input_size (int): It is the number of coefficients in coeffs that we see as the dividend.
    offset (int): The index in coeffs where the dividend is started
    """
    global n_add
    print("coeffs: ", coeffs)
    # print("s:      ", s)
    print("ns_hdt_S: ", ns_hdt_S)
    q = [0] * (len(s)-1) # In general case,  deg(q) = deg(coeffs) - deg(s). 
                         # However, in each round r, deg(coeffs) <= 2^(r+1)-1 and deg(s) = 2^r, hence deg(q) <= 2^r -1 = deg(s) - 1
    for i in range(len(q)):
        if coeffs[~i] == 0: # Coefficients are ordered from the constant term up to the highest degree term. 
                            # Therefore, the last element is the highest degree coefficient. ~i points to the i-th last element in the list
                            # ~ is the bitwise NOT where ~i = -i-1. bitwise NOT is used because its time complexity is O(1)
            continue
        q[~i] = coeffs[~i]  
        for j in range(0,len(ns_hdt_S)-1): # We ignor the highest degree of S (last element of ns_hdt_S) because in this loop we know it will be canceled with coeffs[i]
            coeffs[~i - ns_hdt_S[j]] += coeffs[~i]; n_add += 1
            
        coeffs[~i] = 0 

    return (q, coeffs[:len(q)])

# def b_n_computation(s, s_eval_x, q, b_0):
#     m = len(s)
#     b_1 = copy.deepcopy(b_0)
#     s_eval_y = 0
#     for i in range(len(s)): ## ToDo: Make this 'for-loop' optimized by changing the data structures.
#         if (s[i]==1): 
#             s_eval_y += s_eval_x ** (len(s)-i-1)

#     for i in range(len(q)): 
#         if(q[i]!=0):
#             q[i] *=  s_eval_y
#         b_1[len(b_1) - len(q) + i] += q[i] ## ToDo: changing the coeff data structure makes this operation optimized

#     return(b_1)    
    
    

def fft(f_coeffs, m, a):
    global n_add, n_mult

    dim = m #math.ceil(math.log(len(f_coeffs), 2)) 
    W = initial_basis_computation(a, dim=m)
    print("W: ", W)

    S, ns_hdt_S = S_function_computation(dim=m)
    # R = m - 1 # number of rounds
    print("S:", S)
    print("ns_hdt_S:", ns_hdt_S)
    
    # r = R # r: round indicator
    for r in reversed(range(m)):
        print("r: ", r)

        q, b_0 = divide(f_coeffs, S[r], ns_hdt_S[r])
        print("q: ", q)
        print("b_0: ", b_0)
        # print("b_1: ", b_1)
        f_coeffs = b_0

    quit()
    print("m:", m)
    print("S[0]:", S[0])
    print("W:", W)

    b = [0] * 2**dim
    b[0] = b_0[0]
    b[1] = q[0] + b[0]
    for i in range(2, dim): 
        # print("q[0] * (S[0][0]): ", q[0] * (S[0][0]))
        # print("q[0] * (S[0][1] + W[i-1]): ", q[0] * (S[0][1] + W[i-1]))
        # print("b[0]",  b[0])
        b[i] = q[0] * (S[0][0] * + W[i-1]) + q[0] * (S[0][1]) + b[0]
    print("q:", q)
    
    for i in range(2**dim):
        print(i, b[i])
        

    
