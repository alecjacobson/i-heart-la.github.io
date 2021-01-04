"""
v_i = sum_j w_i,j M_j u_i

where

w: ℝ^(4×4)
M_j: ℝ^(4×4)
u_i: ℝ^4
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo39(w, M, u):
    """
    :param :w : ℝ^(4×4)
    :param :M : ℝ^(4×4)
    :param :u : ℝ^4
    """
    w = np.asarray(w, dtype=np.float64)
    M = np.asarray(M, dtype=np.float64)
    u = np.asarray(u, dtype=np.float64)

    _dim_0 = M.shape[0]
    _dim_1 = u.shape[0]
    assert w.shape == (4, 4)
    assert M.shape == (_dim_0, 4, 4)
    assert u.shape == (_dim_1, 4, )

    v = np.zeros((_dim_1, 4, ))
    for i in range(1, _dim_1+1):
        _sum_0 = np.zeros((4, ))
        for j in range(1, len(M)+1):
            _sum_0 += w[i-1, j-1] * M[j-1] @ u[i-1]
        v[i-1] = _sum_0

    return v


def generateRandomData():
    _dim_0 = np.random.randint(10)
    _dim_1 = np.random.randint(10)
    w = np.random.randn(4, 4)
    M = np.random.randn(_dim_0, 4, 4)
    u = np.random.randn(_dim_1, 4, )
    return w, M, u


if __name__ == '__main__':
    w, M, u = generateRandomData()
    print("w:", w)
    print("M:", M)
    print("u:", u)
    func_value = demo39(w, M, u)
    print("func_value: ", func_value)