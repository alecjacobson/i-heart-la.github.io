"""
n = ∑_T A_T||M_T v_T - [0 -1
                        1  0] M_T u_T||²
where
 
v_i: ℝ^3
u_i: ℝ^3
M_i: ℝ^(2×3)
A_i: ℝ
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo4(v, u, M, A):
    """
    :param :v : ℝ^3
    :param :u : ℝ^3
    :param :M : ℝ^(2×3)
    :param :A : ℝ
    """
    v = np.asarray(v, dtype=np.float64)
    u = np.asarray(u, dtype=np.float64)
    M = np.asarray(M, dtype=np.float64)
    A = np.asarray(A)

    _dim_0 = A.shape[0]
    assert v.shape == (_dim_0, 3, )
    assert u.shape == (_dim_0, 3, )
    assert M.shape == (_dim_0, 2, 3)
    assert A.shape == (_dim_0,)

    _sum_0 = 0
    for T in range(1, len(v)+1):
        _n_0 = np.zeros((2, 2))
        _n_0[0] = [0, -1]
        _n_0[1] = [1, 0]
        _sum_0 += A[T-1] * np.power(np.linalg.norm(M[T-1] @ v[T-1] - _n_0 @ M[T-1] @ u[T-1], 2), 2)
    n = _sum_0

    return n


def generateRandomData():
    _dim_0 = np.random.randint(10)
    v = np.random.randn(_dim_0, 3, )
    u = np.random.randn(_dim_0, 3, )
    M = np.random.randn(_dim_0, 2, 3)
    A = np.random.randn(_dim_0)
    return v, u, M, A


if __name__ == '__main__':
    v, u, M, A = generateRandomData()
    print("v:", v)
    print("u:", u)
    print("M:", M)
    print("A:", A)
    func_value = demo4(v, u, M, A)
    print("func_value: ", func_value)