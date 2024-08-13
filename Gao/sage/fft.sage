load('taylor.sage')

def fft(g_coeffs, m, B):
    '''
    g_coeffs: a list of coefficients of the input function g(x)
    m:  deg(g) < 2^m
    B: basis of the evaluation set
    '''
    # Basis precomputation
    G = [[]] * (m-1); D = [[]] * (m-1)
    G[0], D[0] = sets_computation(B)
    for i in range(1,m-1):
        G[i], D[i] = sets_computation(D[i-1])

    # print("G: ", G)
    # print("D: ", D)
    # eval_basis(G[0], 3)

    for i in range(1, m-1):



    

def eval_basis(G, index):
    index_bin = Integer(index).binary()
    result = 0
    for i, bit in enumerate(reversed(index_bin)):
        # print(i, type(bit), G[i])
        result += int(bit) * G[i]
    # print("result: ", result)

def sets_computation(B):
    G = [0] * (len(B)-1)
    D = [0] * len(G)
    for i in range(len(G)):
        G[i] = B[i] * (B[-1] ** (-1))
        D[i] = (G[i] ** 2) - G[i]
    return G, D


def fft_f(g_b):
    pass