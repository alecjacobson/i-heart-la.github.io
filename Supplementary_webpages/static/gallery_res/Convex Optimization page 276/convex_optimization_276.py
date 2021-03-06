"""
min_(x ∈ ℝ^n) ∑_i ‖A_i x + b_i‖ + (1/2)‖x-`x₀`‖²

where

A_i ∈ ℝ^(m × n)  
`x₀` ∈ ℝ^n  
b_i ∈ ℝ^m  
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class convex_optimization_276ResultType:
    def __init__( self, ret):
        self.ret = ret


def convex_optimization_276(A, x0, b):
    A = np.asarray(A, dtype=np.float64)
    x0 = np.asarray(x0, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)

    dim_0 = A.shape[0]
    m = A.shape[1]
    n = A.shape[2]
    assert A.shape == (dim_0, m, n)
    assert x0.shape == (n,)
    assert b.shape == (dim_0, m, )

    def target_0(x):
        sum_0 = 0
        for i in range(1, len(A)+1):
            sum_0 += np.linalg.norm(A[i-1] @ x + b[i-1], 2)
        return sum_0 + (1 / 2) * np.power(np.linalg.norm(x - x0, 2), 2)
    ret = minimize(target_0, np.zeros(n)).fun
    return convex_optimization_276ResultType(ret)


def generateRandomData():
    dim_0 = np.random.randint(10)
    m = np.random.randint(10)
    n = np.random.randint(10)
    A = np.random.randn(dim_0, m, n)
    x0 = np.random.randn(n)
    b = np.random.randn(dim_0, m, )
    return A, x0, b


if __name__ == '__main__':
    A, x0, b = generateRandomData()
    print("A:", A)
    print("x0:", x0)
    print("b:", b)
    func_value = convex_optimization_276(A, x0, b)
    print("return value: ", func_value.ret)