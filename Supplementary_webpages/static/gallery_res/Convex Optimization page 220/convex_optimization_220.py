"""
`L(x,ν)` = xᵀWx + ∑_i ν_i(x_i²-1)

where

x ∈ ℝ^n
W ∈ ℝ^(n×n)
ν ∈ ℝ^n
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class convex_optimization_220ResultType:
    def __init__( self, L_left_parenthesis_x_comma_ν_right_parenthesis):
        self.L_left_parenthesis_x_comma_ν_right_parenthesis = L_left_parenthesis_x_comma_ν_right_parenthesis


def convex_optimization_220(x, W, ν):
    x = np.asarray(x, dtype=np.float64)
    W = np.asarray(W, dtype=np.float64)
    ν = np.asarray(ν, dtype=np.float64)

    n = x.shape[0]
    assert x.shape == (n,)
    assert W.shape == (n, n)
    assert ν.shape == (n,)

    sum_0 = 0
    for i in range(1, len(ν)+1):
        sum_0 += ν[i-1] * (np.power(x[i-1], 2) - 1)
    L_left_parenthesis_x_comma_ν_right_parenthesis = (x.T.reshape(1, n) @ W @ x).item() + sum_0
    return convex_optimization_220ResultType(L_left_parenthesis_x_comma_ν_right_parenthesis)


def generateRandomData():
    n = np.random.randint(10)
    x = np.random.randn(n)
    W = np.random.randn(n, n)
    ν = np.random.randn(n)
    return x, W, ν


if __name__ == '__main__':
    x, W, ν = generateRandomData()
    print("x:", x)
    print("W:", W)
    print("ν:", ν)
    func_value = convex_optimization_220(x, W, ν)
    print("return value: ", func_value.L_left_parenthesis_x_comma_ν_right_parenthesis)