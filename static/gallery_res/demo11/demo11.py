"""
y_i = (a_i)ᵀ x + w_i
`x_bar` = (sum_i a_i(a_i)ᵀ)^(-1) sum_i y_i a_i

where

a_i: ℝ^n: the measurement vectors  
w_i: ℝ: measurement noise 
x: ℝ^n: measurement noise 
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo11(a, w, x):
    """
    :param :a : ℝ^n: the measurement vectors  
    :param :w : ℝ: measurement noise 
    :param :x : ℝ^n: measurement noise 
    """
    a = np.asarray(a, dtype=np.float64)
    w = np.asarray(w)
    x = np.asarray(x, dtype=np.float64)

    _dim_0 = w.shape[0]
    n = x.shape[0]
    assert a.shape == (_dim_0, n, )
    assert w.shape == (_dim_0,)
    assert x.shape == (n,)

    y = np.zeros(_dim_0)
    for i in range(1, _dim_0+1):
        y[i-1] = (a[i-1]).T @ x + w[i-1]

    _sum_0 = np.zeros((n, n))
    for i in range(1, len(a)+1):
        _sum_0 += a[i-1] @ (a[i-1]).T
    _sum_1 = np.zeros((n, ))
    for i in range(1, len(a)+1):
        _sum_1 += y[i-1] * a[i-1]
    x_bar = np.linalg.inv((_sum_0)) @ _sum_1

    return x_bar


def generateRandomData():
    _dim_0 = np.random.randint(10)
    n = np.random.randint(10)
    a = np.random.randn(_dim_0, n, )
    w = np.random.randn(_dim_0)
    x = np.random.randn(n)
    return a, w, x


if __name__ == '__main__':
    a, w, x = generateRandomData()
    print("a:", a)
    print("w:", w)
    print("x:", x)
    func_value = demo11(a, w, x)
    print("func_value: ", func_value)