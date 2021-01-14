"""
min_(x ∈ ℝ^n) ∑_i ||A_i x + b_i ||_2 +(1/2)||x-`x₀`||^2_2

where

A_i: ℝ^(m × n)  
`x₀`: ℝ^n  
b_i: ℝ^m  
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo10(A, x0, b):
    """
    :param :A : ℝ^(m × n)
    :param :x0 : ℝ^n
    :param :b : ℝ^m
    """
    A = np.asarray(A, dtype=np.float64)
    x0 = np.asarray(x0, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)

    _dim_0 = A.shape[0]
    m = b.shape[1]
    n = x0.shape[0]
    assert A.shape == (_dim_0, m, n)
    assert x0.shape == (n,)
    assert b.shape == (_dim_0, m, )

    def _target_0(x):
        _sum_0 = 0
        for i in range(1, len(A)+1):
            _sum_0 += np.linalg.norm(A[i-1] @ x + b[i-1], 2)
        return _sum_0 + (1 / 2) * np.power(np.linalg.norm(x - x0, 2), 2)
    ret = minimize(_target_0, np.zeros(n)).fun
    return ret


def generateRandomData():
    _dim_0 = np.random.randint(10)
    m = np.random.randint(10)
    n = np.random.randint(10)
    A = np.random.randn(_dim_0, m, n)
    x0 = np.random.randn(n)
    b = np.random.randn(_dim_0, m, )
    return A, x0, b


if __name__ == '__main__':
    A, x0, b = generateRandomData()
    print("A:", A)
    print("x0:", x0)
    print("b:", b)
    func_value = demo10(A, x0, b)
    print("func_value: ", func_value)