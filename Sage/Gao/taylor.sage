
def taylor_expansion(coeffs):
    R = math.ceil(math.log(len(coeffs), 2)) - 1            # R: Number of the rounds
    # coeffs.extend([0]*(2**(R + 1)-len(coeffs)))
    coeffs[:0] = [0]*(2**(R + 1)-len(coeffs))
    # print("> zero padded:  coeffs =",coeffs)
    # print("> Number of rounds:  R =",R)

    r = 0   # r: round counter 
    for r in range(R):
        # print("> r =", r)
        input_size = len(coeffs) / (2**r)
        for b in range(2**r): 
            taylor_module(coeffs, input_size, b) 
        # print(">> coeffs =", coeffs)

    return (coeffs[0::2] + coeffs[1::2]) # Concats g_0 and g_1
    

def taylor_module(coeffs, input_size, b):
    chunk_size = input_size / 4     # Divides the vector of coeffs into four chunks
    for i in range(chunk_size):
        coeffs[ b * input_size + 2 * chunk_size + i ] += coeffs[ b * input_size + 3 * chunk_size + i ]
        coeffs[ b * input_size + 1 * chunk_size + i ] += coeffs[ b * input_size + 2 * chunk_size + i ]
