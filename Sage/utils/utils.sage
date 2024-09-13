def span_basis(B):
    subset = [0] * 2**len(B)
    for i in range(len(subset)):
        subset[i] = an_element_in_basis(B, i)
    return subset

def polynomial_to_string(coeffs):
    string = f"{coeffs[0]}"
    for i, c in enumerate(coeffs[1:]):
        if c:
            if c==1:
                string += f" + x^{i+1}"
            else:
                string += f" + ({c})x^{i+1}"
    return string

def evaluate_polynomial(coeffs, eval_set):
    """
    Evaluates the polynomial represented by coeffs at each point in the evaluation set eval_set.

    Parameters:
    coeffs (list): The list of coefficients of the polynomial, ordered from the constant term up to the highest degree term.
    eval_set (list): The list of points at which to evaluate the polynomial.

    Returns:
    evaluations (list): A list of the polynomial evaluations at each point in eval_set.
    """
    evaluations = [0] * len(eval_set)
    for i in range(len(eval_set)):
        if eval_set[i] == 0:
            evaluations[i] = coeffs[0]
        else:
            x = 1
            for c in coeffs:
                evaluations[i] += c * x
                x *= eval_set[i]
    return evaluations


def an_element_in_basis(B, index):
    """
    Returns the element at the specified index in the basis set B.

    Parameters:
    B (list): The basis set.
    index (int): The index of the element to retrieve from the basis set B.

    Returns:
    result: The element at the specified index in B, s.t., B[index]
    """
    index_bin = Integer(index).binary()
    result = 0
    for i, bit in enumerate(reversed(index_bin)):
        result += int(bit) * B[i]
    return result