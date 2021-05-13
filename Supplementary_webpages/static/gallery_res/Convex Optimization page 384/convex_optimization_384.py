"""
y_i = a_iᵀ x + w_i
x̂ = (∑_i a_i a_iᵀ)⁻¹ ∑_i y_i a_i

where

a_i ∈ ℝ^m: the measurement vectors  
w_i ∈ ℝ: measurement noise 
x ∈ ℝ^m: original vector 
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class convex_optimization_384ResultType:
    def __init__( self, y, x̂):
        self.y = y
        self.x̂ = x̂


def convex_optimization_384(a, w, x):
    """
    :param :a : the measurement vectors  
    :param :w : measurement noise 
    :param :x : original vector 
    """
    a = np.asarray(a, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    x = np.asarray(x, dtype=np.float64)

    dim_0 = w.shape[0]
    m = a.shape[1]
    assert a.shape == (dim_0, m, )
    assert w.shape == (dim_0,)
    assert x.shape == (m,)

    y = np.zeros(dim_0)
    for i in range(1, dim_0+1):
        y[i-1] = (a[i-1].T.reshape(1, m) @ x).item() + w[i-1]
    sum_0 = np.zeros((m, m))
    for i in range(1, len(a)+1):
        sum_0 += (a[i-1]).reshape(m, 1) @ a[i-1].T.reshape(1, m)
    sum_1 = np.zeros((m, ))
    for i in range(1, len(y)+1):
        sum_1 += y[i-1] * a[i-1]
    x̂ = np.linalg.inv((sum_0)) @ sum_1
    return convex_optimization_384ResultType(y, x̂)


def generateRandomData():
    dim_0 = np.random.randint(10)
    m = np.random.randint(10)
    a = np.random.randn(dim_0, m, )
    w = np.random.randn(dim_0)
    x = np.random.randn(m)
    return a, w, x


if __name__ == '__main__':
    a, w, x = generateRandomData()
    print("a:", a)
    print("w:", w)
    print("x:", x)
    func_value = convex_optimization_384(a, w, x)
    print("return value: ", func_value.x̂)