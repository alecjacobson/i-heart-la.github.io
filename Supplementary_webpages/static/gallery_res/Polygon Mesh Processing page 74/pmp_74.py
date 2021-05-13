"""
`E_LSCM` = ∑_T A_T‖M_T v_T - [0 -1
                              1  0] M_T u_T‖²
where
 
v_i ∈ ℝ^3
u_i ∈ ℝ^3
M_i ∈ ℝ^(2×3)
A_i ∈ ℝ
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class pmp_74ResultType:
    def __init__( self, E_LSCM):
        self.E_LSCM = E_LSCM


def pmp_74(v, u, M, A):
    v = np.asarray(v, dtype=np.float64)
    u = np.asarray(u, dtype=np.float64)
    M = np.asarray(M, dtype=np.float64)
    A = np.asarray(A, dtype=np.float64)

    dim_0 = A.shape[0]
    assert v.shape == (dim_0, 3, )
    assert u.shape == (dim_0, 3, )
    assert M.shape == (dim_0, 2, 3)
    assert A.shape == (dim_0,)

    sum_0 = 0
    for T in range(1, len(u)+1):
        E_LSCM_0 = np.zeros((2, 2))
        E_LSCM_0[0] = [0, -1]
        E_LSCM_0[1] = [1, 0]
        sum_0 += A[T-1] * np.power(np.linalg.norm(M[T-1] @ v[T-1] - E_LSCM_0 @ M[T-1] @ u[T-1], 2), 2)
    E_LSCM = sum_0
    return pmp_74ResultType(E_LSCM)


def generateRandomData():
    dim_0 = np.random.randint(10)
    v = np.random.randn(dim_0, 3, )
    u = np.random.randn(dim_0, 3, )
    M = np.random.randn(dim_0, 2, 3)
    A = np.random.randn(dim_0)
    return v, u, M, A


if __name__ == '__main__':
    v, u, M, A = generateRandomData()
    print("v:", v)
    print("u:", u)
    print("M:", M)
    print("A:", A)
    func_value = pmp_74(v, u, M, A)
    print("return value: ", func_value.E_LSCM)