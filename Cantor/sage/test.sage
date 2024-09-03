load('fft.sage')

def test_fft(a):
    # f_coeffs = [a^0,0,0,a^0,0,0,0,0,0,a^0,0,a^5,a^0,a^2,a,a^0] # f(x) = x^15 + x^12 + x^6 + a^5x^4 + x^3 + a^2x^2 + ax + 1   
    # for i, c in enumerate(reversed(f)):
    #     print(f"C_{i} = {c}")
    g_coeffs = [a^3,a^2,a,1] # g(x) = a^3 + a^2x + ax^2 + x^3
    m = 2
    fft(g_coeffs, m, a)
    print("a^13 + a: ", a^13 + a)
    # print('a^8 + a^4 + a^2 + a =', a^8 + a^4 + a^2 + a)

def generate_a_map(a, dim):
    for i in range(2**dim):
        print(f"a^{i} = {a**i}")

if __name__ == "__main__":
    F.<xx> = QQ[]
    FF.<a> = GF(2**4, modulus= xx**4 + xx + 1)
    # generate_a_map(a, 4)
    test_fft(a)