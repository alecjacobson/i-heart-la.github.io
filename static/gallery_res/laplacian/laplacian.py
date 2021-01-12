"""
from trigonometry: cot

L_i,j = { cot(α_ij) + cot(β_ij) if j ∈ N(i)

L_i,i = -sum_(k for k != i) L_i,k

where

L: ℝ^(n×n)
α: ℝ^(n×n)
β: ℝ^(n×n)
N: ℤ -> {ℤ}
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def laplacian(α, β, N):
    """
    :param :α : ℝ^(n×n)
    :param :β : ℝ^(n×n)
    :param :N : ℤ -> {ℤ}
    """
    α = np.asarray(α, dtype=np.float64)
    β = np.asarray(β, dtype=np.float64)

    n = β.shape[1]
    assert α.shape == (n, n)
    assert β.shape == (n, n)

    _Lij_0 = []
    _Lvals_0 = []
    for i in range(1, n+1):
        for j in range(1, n+1):
            if (j) in N(i):
                _Lij_0.append((i-1, j-1))
                _Lvals_0.append(1/np.tan(α[i-1, j-1]) + 1/np.tan(β[i-1, j-1]))
    _sparse_0 = scipy.sparse.coo_matrix((_Lvals_0, np.asarray(_Lij_0).T), shape=(n, n))
    L = _sparse_0

    for i in range(1, n+1):
        _sum_0 = 0
        for k in range(1, L.shape[0]+1):
            if(k != i):
                _sum_0 += L.tocsr()[i-1, k-1]
        _Lij_0.append((i - 1, i - 1))
        _Lvals_0.append(-_sum_0)
    L = scipy.sparse.coo_matrix((_Lvals_0, np.asarray(_Lij_0).T), shape=(n, n))


    return L


def generateRandomData():
    n = np.random.randint(10)
    α = np.random.randn(n, n)
    β = np.random.randn(n, n)
    def N(p0):
        tmp = []
        tmp_0 = np.random.randint(1, 10)
        for i in range(tmp_0):
            tmp.append((np.random.randint(10)))
        return tmp
    return α, β, N


if __name__ == '__main__':
    α, β, N = generateRandomData()
    print("α:", α)
    print("β:", β)
    print("N:", N)
    func_value = laplacian(α, β, N)
    print("func_value: ", func_value)