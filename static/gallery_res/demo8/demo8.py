"""
`I(X;Y)` = ∑_i ∑_j x_j p_i,j log(p_ij/∑_k x_k p_ik)

where

x: ℝ^n
p: ℝ^(n×m)

"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo8(x, p):
    """
    :param :x : ℝ^n
    :param :p : ℝ^(n×m)
    """
    x = np.asarray(x, dtype=np.float64)
    p = np.asarray(p, dtype=np.float64)

    n = p.shape[0]
    m = p.shape[1]
    assert x.shape == (n,)
    assert p.shape == (n, m)

    _sum_0 = 0
    for i in range(1, len(p)+1):
        _sum_1 = 0
        for j in range(1, len(x)+1):
            _sum_1 += x[j-1] * p[i-1, j-1] * np.log10(p[i-1, j-1] / _sum_2)
        _sum_0 += _sum_1
    I_left_parenthesis_X_semicolon_Y_right_parenthesis = _sum_0

    return I_left_parenthesis_X_semicolon_Y_right_parenthesis


def generateRandomData():
    n = np.random.randint(10)
    m = np.random.randint(10)
    x = np.random.randn(n)
    p = np.random.randn(n, m)
    return x, p


if __name__ == '__main__':
    x, p = generateRandomData()
    print("x:", x)
    print("p:", p)
    func_value = demo8(x, p)
    print("func_value: ", func_value)