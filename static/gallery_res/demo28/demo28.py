"""
`G_σ(s_i^k)` = ∑_j l_j exp(-dist(`bᵢ`, b_j)/(2σ²))(s_j)^k

where

l_j: ℝ : the length of bj
dist: ℝ^n, ℝ^n -> ℝ : measures the geodesic distance between the centers of bi and bj along the boundary
σ: ℝ
`bᵢ`: ℝ^n
b_j: ℝ^n
s_j: ℝ : unit direction vector of bi
k: ℝ : iteration number
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo28(l, dist, σ, bᵢ, b, s, k):
    """
    :param :l : ℝ : the length of bj
    :param :dist : ℝ^n, ℝ^n -> ℝ : measures the geodesic distance between the centers of bi and bj along the boundary
    :param :σ : ℝ
    :param :bᵢ : ℝ^n
    :param :b : ℝ^n
    :param :s : ℝ : unit direction vector of bi
    :param :k : ℝ : iteration number
    """
    l = np.asarray(l)
    bᵢ = np.asarray(bᵢ, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    s = np.asarray(s)

    _dim_0 = l.shape[0]
    n = b.shape[1]
    assert l.shape == (_dim_0,)
    assert np.ndim(σ) == 0
    assert bᵢ.shape == (n,)
    assert b.shape == (_dim_0, n, )
    assert s.shape == (_dim_0,)
    assert np.ndim(k) == 0

    _sum_0 = 0
    for j in range(1, len(b)+1):
        _sum_0 += l[j-1] * np.exp(-dist(bᵢ, b[j-1]) / (2 * np.power(σ, 2))) * np.power((s[j-1]), k)
    G_σ_left_parenthesis_s_i_circumflex_accent_k_right_parenthesis = _sum_0

    return G_σ_left_parenthesis_s_i_circumflex_accent_k_right_parenthesis


def generateRandomData():
    σ = np.random.randn()
    k = np.random.randn()
    _dim_0 = np.random.randint(10)
    n = np.random.randint(10)
    l = np.random.randn(_dim_0)
    def dist(p0, p1):
        return np.random.randn()
    bᵢ = np.random.randn(n)
    b = np.random.randn(_dim_0, n, )
    s = np.random.randn(_dim_0)
    return l, dist, σ, bᵢ, b, s, k


if __name__ == '__main__':
    l, dist, σ, bᵢ, b, s, k = generateRandomData()
    print("l:", l)
    print("dist:", dist)
    print("σ:", σ)
    print("bᵢ:", bᵢ)
    print("b:", b)
    print("s:", s)
    print("k:", k)
    func_value = demo28(l, dist, σ, bᵢ, b, s, k)
    print("func_value: ", func_value)