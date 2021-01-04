"""
`G_σ(s_i^k)` = sum_j l_j exp(-`dist`(`b_i`, b_j)/(2`σ`^2))(s_j)^k

where

l_j: ℝ : the length of b_j
`dist`: ℝ, ℝ -> ℝ : measures the geodesic distance between the centers of b_i and b_j along the boundary
`σ`: ℝ
`b_i`: ℝ
b_j: ℝ
s_j: ℝ : unit direction vector of b_i
k: ℝ : iteration number
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo28(l, dist, σ, b_i, b, s, k):
    """
    :param :l : ℝ : the length of b_j
    :param :dist : ℝ, ℝ -> ℝ : measures the geodesic distance between the centers of b_i and b_j along the boundary
    :param :σ : ℝ
    :param :b_i : ℝ
    :param :b : ℝ
    :param :s : ℝ : unit direction vector of b_i
    :param :k : ℝ : iteration number
    """
    l = np.asarray(l)
    b = np.asarray(b)
    s = np.asarray(s)

    _dim_0 = l.shape[0]
    assert l.shape == (_dim_0,)
    assert np.ndim(σ) == 0
    assert np.ndim(b_i) == 0
    assert b.shape == (_dim_0,)
    assert s.shape == (_dim_0,)
    assert np.ndim(k) == 0

    _sum_0 = 0
    for j in range(1, len(b)+1):
        _sum_0 += l[j-1] * np.exp(-dist(b_i, b[j-1]) / (2 * np.power(σ, 2))) * np.power((s[j-1]), k)
    G_σ_left_parenthesis_s_i_circumflex_accent_k_right_parenthesis = _sum_0

    return G_σ_left_parenthesis_s_i_circumflex_accent_k_right_parenthesis


def generateRandomData():
    σ = np.random.randn()
    b_i = np.random.randn()
    k = np.random.randn()
    _dim_0 = np.random.randint(10)
    l = np.random.randn(_dim_0)
    def dist(p0, p1):
        return np.random.randn()
    b = np.random.randn(_dim_0)
    s = np.random.randn(_dim_0)
    return l, dist, σ, b_i, b, s, k


if __name__ == '__main__':
    l, dist, σ, b_i, b, s, k = generateRandomData()
    print("l:", l)
    print("dist:", dist)
    print("σ:", σ)
    print("b_i:", b_i)
    print("b:", b)
    print("s:", s)
    print("k:", k)
    func_value = demo28(l, dist, σ, b_i, b, s, k)
    print("func_value: ", func_value)