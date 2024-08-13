import math
import copy

n_add = 0
n_mult = 0

def initial_basis_computation(a, dim):
    B = [1] * dim
    for i in range(1, dim):
        b = a
        while true:
            if b**2 + b == B[i-1]:
                break
            b *= a
        B[i] = b
    return B

def S_function_computation(dim):
    S = [[]] * (dim+1)
    # nz_S: Non zero coeffs indices in S. Note that index 0 is the coefficient of the largest degree
    nz_S = [[]] * (dim+1)                   
    for m in range(dim+1):
        S[m] = [0] * (2**m + 1)
        nz_S[m] = []
        for i in range(m+1):
            S[m][2**m - 2**i] = math.comb(m,i) % 2
            if (math.comb(m,i) % 2 == 1): nz_S[m].append(2**m - 2**i) 
            
    return S, nz_S

def divide(coeffs, s, nz_s):
    global n_add

    print("coeffs: ", coeffs)
    print("s:      ", s)
    print("nz_s: ", nz_s)
    q = [0] * (len(s)-1)
    for i in range(len(q)):
        if coeffs[i] == 0:
            continue
        q[i] = coeffs[i]
        for j in range(0,len(nz_s)-1):
                coeffs[i + nz_s[j]] += coeffs[i]; n_add += 1
            
        coeffs[i] = 0

    return (q, coeffs[-len(q):])

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
    
    

def fft(f_coeffs, a):
    global n_add, n_mult

    dim = math.ceil(math.log(len(f_coeffs), 2)) 
    W = initial_basis_computation(a, dim)
    print("W: ", W)
    
    R = dim - 1 # number of rounds
    S, nz_S = S_function_computation(R)
    print("S:", S)
    print("nz_S:", nz_S)
    
    # r = R # r: round indicator
    for r in reversed(range(R+1)):
        print("r: ", r)

        q, b_0 = divide(copy.deepcopy(f_coeffs), S[r], nz_S[r])
        print("q: ", q)
        print("b_0: ", b_0)
        # print("b_1: ", b_1)
        f_coeffs = b_0

    print("R:", R)
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
        

    
