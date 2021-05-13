"""
L_i,j = { w_i,j if (i,j) ∈ E
          0 otherwise

L_i,i = -sum_(l for l ≠ i) L_i,l

where

L ∈ ℝ^(n×n)
w ∈ ℝ^(n×n): edge weight matrix
E ∈ {ℤ²} index: edges
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class course_parameterizationResultType:
    def __init__( self, L):
        self.L = L


def course_parameterization(w, E):
    """
    :param :w : edge weight matrix
    :param :E : edges
    """
    w = np.asarray(w, dtype=np.float64)

    n = w.shape[1]
    assert w.shape == (n, n)
    assert isinstance(E, list) and len(E) > 0
    assert len(E[0]) == 2

    Lij_0 = []
    Lvals_0 = []
    for i in range(1, n+1):
        for j in range(1, n+1):
            if (i-1, j-1) in E:
                Lij_0.append((i-1, j-1))
                Lvals_0.append(w[i-1, j-1])
    sparse_0 = scipy.sparse.coo_matrix((Lvals_0, np.asarray(Lij_0).T), shape=(n, n))
    L = sparse_0
    for i in range(1, n+1):
        sum_0 = 0
        for l in range(1, L.shape[1]+1):
            if(l != i):
                sum_0 += L.tocsr()[i-1, l-1]
        Lij_0.append((i - 1, i - 1))
        Lvals_0.append(-sum_0)
    L = scipy.sparse.coo_matrix((Lvals_0, np.asarray(Lij_0).T), shape=(n, n))

    return course_parameterizationResultType(L)


def generateRandomData():
    n = np.random.randint(10)
    w = np.random.randn(n, n)
    E = []
    dim_0 = np.random.randint(1, 10)
    for i in range(dim_0):
        E.append((np.random.randint(10), np.random.randint(10)))
    return w, E


if __name__ == '__main__':
    w, E = generateRandomData()
    print("w:", w)
    print("E:", E)
    func_value = course_parameterization(w, E)
    print("return value: ", func_value.L)