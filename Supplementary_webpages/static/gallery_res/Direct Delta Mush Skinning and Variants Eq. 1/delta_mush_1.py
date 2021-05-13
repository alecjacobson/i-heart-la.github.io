"""
v_i = ∑_j w_i,j M_j u_i

where

w ∈ ℝ^(n×m)
M_j ∈ ℝ^(4×4)
u_i ∈ ℝ^4
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class delta_mush_1ResultType:
    def __init__( self, v):
        self.v = v


def delta_mush_1(w, M, u):
    w = np.asarray(w, dtype=np.float64)
    M = np.asarray(M, dtype=np.float64)
    u = np.asarray(u, dtype=np.float64)

    n = w.shape[0]
    m = w.shape[1]
    dim_0 = M.shape[0]
    dim_1 = u.shape[0]
    assert w.shape == (n, m)
    assert M.shape == (dim_0, 4, 4)
    assert u.shape == (dim_1, 4, )
    assert dim_0 == m 
    assert dim_1 == n 

    v = np.zeros((dim_1, 4, ))
    for i in range(1, dim_1+1):
        sum_0 = np.zeros((4, ))
        for j in range(1, len(w)+1):
            sum_0 += w[i-1, j-1] * M[j-1] @ u[i-1]
        v[i-1] = sum_0
    return delta_mush_1ResultType(v)


def generateRandomData():
    n = np.random.randint(10)
    dim_1 = n
    m = np.random.randint(10)
    dim_0 = m
    w = np.random.randn(n, m)
    M = np.random.randn(dim_0, 4, 4)
    u = np.random.randn(dim_1, 4, )
    return w, M, u


if __name__ == '__main__':
    w, M, u = generateRandomData()
    print("w:", w)
    print("M:", M)
    print("u:", u)
    func_value = delta_mush_1(w, M, u)
    print("return value: ", func_value.v)