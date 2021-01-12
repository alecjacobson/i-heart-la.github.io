"""
min_(u ∈ ℝ^6) uᵀ(∑_i [x_i×n̂_i; n̂_i][(x_i×n̂_i)ᵀ n̂_iᵀ])u - 2uᵀ(∑_i [x_i×n̂_i; n̂_i]n̂_iᵀ(p_i-x_i)) + ∑_i(p_i-x_i)ᵀn̂_i n̂_iᵀ(p_i-x_i)

where

x_i: ℝ^3 
n̂_i: ℝ^3  
p_i: ℝ^3  
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo25(x, n̂, p):
    """
    :param :x : ℝ^3
    :param :n̂ : ℝ^3
    :param :p : ℝ^3
    """
    x = np.asarray(x, dtype=np.float64)
    n̂ = np.asarray(n̂, dtype=np.float64)
    p = np.asarray(p, dtype=np.float64)

    _dim_0 = x.shape[0]
    assert x.shape == (_dim_0, 3, )
    assert n̂.shape == (_dim_0, 3, )
    assert p.shape == (_dim_0, 3, )

    def _target_0(u):
        _sum_0 = np.zeros((6, 6))
        for i in range(1, len(x)+1):
            _ret_0 = np.vstack(((np.cross(x[i-1], n̂[i-1])).reshape(3, 1), (n̂[i-1]).reshape(3, 1)))
            _ret_1 = np.hstack(((np.cross(x[i-1], n̂[i-1])).T.reshape(1, 3), n̂[i-1].T.reshape(1, 3)))
            _sum_0 += _ret_0 @ _ret_1
        _sum_1 = np.zeros((6, ))
        for i in range(1, len(n̂)+1):
            _ret_2 = np.vstack(((np.cross(x[i-1], n̂[i-1])).reshape(3, 1), (n̂[i-1]).reshape(3, 1)))
            _sum_1 += _ret_2 @ n̂[i-1].T.reshape(1, 3) @ (p[i-1] - x[i-1])
        _sum_2 = 0
        for i in range(1, len(n̂)+1):
            _sum_2 += (p[i-1] - x[i-1]).T.reshape(1, 3) @ n̂[i-1] * n̂[i-1].T.reshape(1, 3) @ (p[i-1] - x[i-1])
        return u.T.reshape(1, 6) @ (_sum_0) @ u - 2 * u.T.reshape(1, 6) @ (_sum_1) + _sum_2
    ret = minimize(_target_0, np.zeros(6)).fun
    return ret


def generateRandomData():
    _dim_0 = np.random.randint(10)
    x = np.random.randn(_dim_0, 3, )
    n̂ = np.random.randn(_dim_0, 3, )
    p = np.random.randn(_dim_0, 3, )
    return x, n̂, p


if __name__ == '__main__':
    x, n̂, p = generateRandomData()
    print("x:", x)
    print("n̂:", n̂)
    print("p:", p)
    func_value = demo25(x, n̂, p)
    print("func_value: ", func_value)