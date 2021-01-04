"""
sum_i α_i + 1/M sum_i sum_j (f(X_i,j)/`p_c`(X_i,j) - (sum_k α_k p_k X_i,j)/`p_c`(X_i,j))

where


α: ℝ^m
p: ℝ^m
X: ℝ^(m×n)
M: ℝ  
f: ℝ -> ℝ 
`p_c`: ℝ -> ℝ 
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo29(α, p, X, M, f, p_c):
    """
    :param :α : ℝ^m
    :param :p : ℝ^m
    :param :X : ℝ^(m×n)
    :param :M : ℝ
    :param :f : ℝ -> ℝ
    :param :p_c : ℝ -> ℝ
    """
    α = np.asarray(α, dtype=np.float64)
    p = np.asarray(p, dtype=np.float64)
    X = np.asarray(X, dtype=np.float64)

    m = X.shape[0]
    n = X.shape[1]
    assert α.shape == (m,)
    assert p.shape == (m,)
    assert X.shape == (m, n)
    assert np.ndim(M) == 0

    _sum_0 = 0
    for i in range(1, len(α)+1):
        _sum_0 += α[i-1]
    _sum_1 = 0
    for i in range(1, len(X)+1):
        _sum_2 = 0
        for j in range(1, len(X)+1):
            _sum_3 = 0
            for k in range(1, len(p)+1):
                _sum_3 += α[k-1] * p[k-1] * X[i-1, j-1]
            _sum_2 += (f(X[i-1, j-1]) / p_c(X[i-1, j-1]) - (_sum_3) / p_c(X[i-1, j-1]))
        _sum_1 += _sum_2
    ret = _sum_0 + 1 / M * _sum_1
    return ret


def generateRandomData():
    M = np.random.randn()
    m = np.random.randint(10)
    n = np.random.randint(10)
    α = np.random.randn(m)
    p = np.random.randn(m)
    X = np.random.randn(m, n)
    def f(p0):
        return np.random.randn()
    def p_c(p0):
        return np.random.randn()
    return α, p, X, M, f, p_c


if __name__ == '__main__':
    α, p, X, M, f, p_c = generateRandomData()
    print("α:", α)
    print("p:", p)
    print("X:", X)
    print("M:", M)
    print("f:", f)
    print("p_c:", p_c)
    func_value = demo29(α, p, X, M, f, p_c)
    print("func_value: ", func_value)