"""
min_(u ∈ ℝ^6) uᵀ(sum_i [x_i×n_i; n_i][(x_i×n_i)ᵀ n_iᵀ])u - 2uᵀ(sum_i [x_i×n_i; n_i]n_iᵀ(p_i-x_i)) + sum_i(p_i-x_i)ᵀn_i n_iᵀ(p_i-x_i)

where

x_i: ℝ^3 
n_i: ℝ^3  
p_i: ℝ^3  
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo25(x, n, p):
    """
    :param :x : ℝ^3
    :param :n : ℝ^3
    :param :p : ℝ^3
    """
    x = np.asarray(x, dtype=np.float64)
    n = np.asarray(n, dtype=np.float64)
    p = np.asarray(p, dtype=np.float64)

    _dim_0 = x.shape[0]
    assert x.shape == (_dim_0, 3, )
    assert n.shape == (_dim_0, 3, )
    assert p.shape == (_dim_0, 3, )

    def _target_0(u):
        _sum_0 = np.zeros((6, 6))
        for i in range(1, len(x)+1):
            _ret_0 = np.vstack((np.cross(x[i-1], n[i-1]), n[i-1]))
            _ret_1 = np.hstack(((np.cross(x[i-1], n[i-1])).T, n[i-1].T))
            _sum_0 += _ret_0 @ _ret_1
        _sum_1 = np.zeros((6, ))
        for i in range(1, len(p)+1):
            _ret_2 = np.vstack((np.cross(x[i-1], n[i-1]), n[i-1]))
            _sum_1 += _ret_2 @ n[i-1].T @ (p[i-1] - x[i-1])
        _sum_2 = 0
        for i in range(1, len(p)+1):
            _sum_2 += (p[i-1] - x[i-1]).T @ n[i-1] * n[i-1].T @ (p[i-1] - x[i-1])
        return u.T @ (_sum_0) @ u - 2 * u.T @ (_sum_1) + _sum_2
    ret = minimize(_target_0, np.zeros(6)).fun
    return ret


def generateRandomData():
    _dim_0 = np.random.randint(10)
    x = np.random.randn(_dim_0, 3, )
    n = np.random.randn(_dim_0, 3, )
    p = np.random.randn(_dim_0, 3, )
    return x, n, p


if __name__ == '__main__':
    x, n, p = generateRandomData()
    print("x:", x)
    print("n:", n)
    print("p:", p)
    func_value = demo25(x, n, p)
    print("func_value: ", func_value)