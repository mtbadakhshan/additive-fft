
def taylor_expansion(coeffs):
    R = math.ceil(math.log(len(coeffs), 2)) - 1            # R: Number of the rounds
    coeffs.extend([0]*(2**(R + 1)-len(coeffs)))
    print("> zero padded:  coeffs =",coeffs)

    r = 0   # r: round number (r = 0 is the input round)
    while r < R:
        r += 1
        print("> r =", r)
        input_size = len(coeffs) / (r)
        for b in range(r): 
            taylor_module(coeffs, input_size, b) 
        print(">> coeffs =", coeffs)

    return (coeffs[0::2], coeffs[1::2]) 
    

def taylor_module(coeffs, input_size, b):
    chunk_size = input_size / 4     # Divides the vector of coeffs into four chunks
    for i in range(chunk_size):
        coeffs[ b * input_size + 2 * chunk_size + i ] += coeffs[ b * input_size + 3 * chunk_size + i ]
        coeffs[ b * input_size + 1 * chunk_size + i ] += coeffs[ b * input_size + 2 * chunk_size + i ]
