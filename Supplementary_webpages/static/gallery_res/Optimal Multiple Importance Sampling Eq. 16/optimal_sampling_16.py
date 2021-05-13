"""
∑_i α_i + 1/M ∑_i ∑_(j for j ⩽ n_i) (f(X_i,j)/`p_c`(X_i,j) - (∑_k α_k p_k(X_i,j))/`p_c`(X_i,j))

where


α ∈ ℝ^N
p_j ∈ ℝ → ℝ 
X ∈ ℝ^(N×m)
M ∈ ℝ  
n_i ∈ ℝ
f: ℝ → ℝ 
`p_c`: ℝ → ℝ 
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class optimal_sampling_16ResultType:
    def __init__( self, ret):
        self.ret = ret


def optimal_sampling_16(α, p, X, M, n, f, p_c):
    """
    :param :f : ℝ → ℝ
    :param :p_c : ℝ → ℝ
    """
    α = np.asarray(α, dtype=np.float64)
    X = np.asarray(X, dtype=np.float64)
    n = np.asarray(n, dtype=np.float64)

    dim_0 = n.shape[0]
    N = α.shape[0]
    m = X.shape[1]
    dim_1 = p.shape[0]
    assert α.shape == (N,)
    assert X.shape == (N, m)
    assert np.ndim(M) == 0
    assert n.shape == (dim_0,)
    assert N == dim_1 

    sum_0 = 0
    for i in range(1, len(α)+1):
        sum_0 += α[i-1]
    sum_1 = 0
    for i in range(1, len(X)+1):
        sum_2 = 0
        for j in range(1, len(X)+1):
            sum_3 = 0
            for k in range(1, len(α)+1):
                sum_3 += α[k-1] * p[k-1](X[i-1, j-1])
            if(j <= n[i-1]):
                sum_2 += (f(X[i-1, j-1]) / p_c(X[i-1, j-1]) - (sum_3) / p_c(X[i-1, j-1]))
        sum_1 += sum_2
    ret = sum_0 + 1 / M * sum_1
    return optimal_sampling_16ResultType(ret)


def generateRandomData():
    M = np.random.randn()
    dim_0 = np.random.randint(10)
    N = np.random.randint(10)
    dim_1 = N
    m = np.random.randint(10)
    α = np.random.randn(N)
    p = []
    for i in range(dim_1):
        def p_f(p0):
            return np.random.randn()
        p.append(p_f)
    p = np.asarray(p)
    X = np.random.randn(N, m)
    n = np.random.randn(dim_0)
    def f(p0):
        return np.random.randn()
    def p_c(p0):
        return np.random.randn()
    return α, p, X, M, n, f, p_c


if __name__ == '__main__':
    α, p, X, M, n, f, p_c = generateRandomData()
    print("α:", α)
    print("p:", p)
    print("X:", X)
    print("M:", M)
    print("n:", n)
    print("f:", f)
    print("p_c:", p_c)
    func_value = optimal_sampling_16(α, p, X, M, n, f, p_c)
    print("return value: ", func_value.ret)