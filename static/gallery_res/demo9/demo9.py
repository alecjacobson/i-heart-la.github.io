"""
`L(x,v)` = xᵀWx + ∑_i v_i(x_i^2-1)

where

x: ℝ^n
W: ℝ^(n×n)
v: ℝ^n
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo9(x, W, v):
    """
    :param :x : ℝ^n
    :param :W : ℝ^(n×n)
    :param :v : ℝ^n
    """
    x = np.asarray(x, dtype=np.float64)
    W = np.asarray(W, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)

    n = v.shape[0]
    assert x.shape == (n,)
    assert W.shape == (n, n)
    assert v.shape == (n,)

    _sum_0 = 0
    for i in range(1, len(x)+1):
        _sum_0 += v[i-1] * (np.power(x[i-1], 2) - 1)
    L_left_parenthesis_x_comma_v_right_parenthesis = x.T @ W @ x + _sum_0

    return L_left_parenthesis_x_comma_v_right_parenthesis


def generateRandomData():
    n = np.random.randint(10)
    x = np.random.randn(n)
    W = np.random.randn(n, n)
    v = np.random.randn(n)
    return x, W, v


if __name__ == '__main__':
    x, W, v = generateRandomData()
    print("x:", x)
    print("W:", W)
    print("v:", v)
    func_value = demo9(x, W, v)
    print("func_value: ", func_value)