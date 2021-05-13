"""
`I(X;Y)` = ∑_i ∑_j x_j p_i,j log₂(p_i,j/∑_k x_k p_i,k)

where

x ∈ ℝ^n
p ∈ ℝ^(m×n)
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class convex_optimization_208ResultType:
    def __init__( self, I_left_parenthesis_X_semicolon_Y_right_parenthesis):
        self.I_left_parenthesis_X_semicolon_Y_right_parenthesis = I_left_parenthesis_X_semicolon_Y_right_parenthesis


def convex_optimization_208(x, p):
    x = np.asarray(x, dtype=np.float64)
    p = np.asarray(p, dtype=np.float64)

    n = x.shape[0]
    m = p.shape[0]
    assert x.shape == (n,)
    assert p.shape == (m, n)

    sum_0 = 0
    for i in range(1, len(p)+1):
        sum_1 = 0
        for j in range(1, len(x)+1):
            sum_2 = 0
            for k in range(1, len(x)+1):
                sum_2 += x[k-1] * p[i-1, k-1]
            sum_1 += x[j-1] * p[i-1, j-1] * np.log2(p[i-1, j-1] / sum_2)
        sum_0 += sum_1
    I_left_parenthesis_X_semicolon_Y_right_parenthesis = sum_0
    return convex_optimization_208ResultType(I_left_parenthesis_X_semicolon_Y_right_parenthesis)


def generateRandomData():
    n = np.random.randint(10)
    m = np.random.randint(10)
    x = np.random.randn(n)
    p = np.random.randn(m, n)
    return x, p


if __name__ == '__main__':
    x, p = generateRandomData()
    print("x:", x)
    print("p:", p)
    func_value = convex_optimization_208(x, p)
    print("return value: ", func_value.I_left_parenthesis_X_semicolon_Y_right_parenthesis)