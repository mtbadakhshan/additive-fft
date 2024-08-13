load('taylor.sage')
load('fft.sage')


def test_taylor(a):
    f_coeffs = [1,a,a^2,a^3,a^4,a^5,a^6,a^7]
    print("Input: f_coeffs = ", f_coeffs)
    g_0, g_1 = taylor_expansion(f_coeffs)
    print(f"g_0 = {g_0}, g_1 = {g_1}")

def test_fft(a):
    g_coeffs = [1,a,a^2,a^3,a^4,a^5,a^6,a^7]
    m = 8
    B = [a**i for i in range(m)]
    print("Inputs: g_coeffs = ", g_coeffs)
    print("               m = ", m)
    print("               B = ", B)
    fft(g_coeffs, m, B)



if __name__ == "__main__":
    F.<xx> = QQ[]
    FF.<a> = GF(2**8, modulus= xx**8 + xx**4 + xx**3 + xx**2 +1)

    # test_taylor(a)
    test_fft(a)