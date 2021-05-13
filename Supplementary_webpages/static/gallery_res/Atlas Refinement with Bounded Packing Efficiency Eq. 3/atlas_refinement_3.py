"""
`G_σ(s^k_i)` = ∑_j l_j exp(-dist(`bᵢ`, b_j)²/(2σ²)) `s^k`_j

where

l_j ∈ ℝ : the length of bj
dist: ℝ^n, ℝ^n → ℝ : measures the geodesic distance 
σ ∈ ℝ
`bᵢ` ∈ ℝ^n
b_j ∈ ℝ^n
`s^k`_j ∈ ℝ^n : direction vector
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class atlas_refinement_3ResultType:
    def __init__( self, G_σ_left_parenthesis_s_circumflex_accent_k_i_right_parenthesis):
        self.G_σ_left_parenthesis_s_circumflex_accent_k_i_right_parenthesis = G_σ_left_parenthesis_s_circumflex_accent_k_i_right_parenthesis


def atlas_refinement_3(l, dist, σ, bᵢ, b, s_circumflex_accent_k):
    """
    :param :l : the length of bj
    :param :dist : ℝ^n, ℝ^n → ℝ : measures the geodesic distance 
    :param :s_circumflex_accent_k : direction vector
    """
    l = np.asarray(l, dtype=np.float64)
    bᵢ = np.asarray(bᵢ, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    s_circumflex_accent_k = np.asarray(s_circumflex_accent_k, dtype=np.float64)

    dim_0 = l.shape[0]
    n = bᵢ.shape[0]
    assert l.shape == (dim_0,)
    assert np.ndim(σ) == 0
    assert bᵢ.shape == (n,)
    assert b.shape == (dim_0, n, )
    assert s_circumflex_accent_k.shape == (dim_0, n, )

    sum_0 = np.zeros((n, ))
    for j in range(1, len(s_circumflex_accent_k)+1):
        sum_0 += l[j-1] * np.exp(-np.power(dist(bᵢ, b[j-1]), 2) / (2 * np.power(σ, 2))) * s_circumflex_accent_k[j-1]
    G_σ_left_parenthesis_s_circumflex_accent_k_i_right_parenthesis = sum_0
    return atlas_refinement_3ResultType(G_σ_left_parenthesis_s_circumflex_accent_k_i_right_parenthesis)


def generateRandomData():
    σ = np.random.randn()
    dim_0 = np.random.randint(10)
    n = np.random.randint(10)
    l = np.random.randn(dim_0)
    def dist(p0, p1):
        return np.random.randn()
    bᵢ = np.random.randn(n)
    b = np.random.randn(dim_0, n, )
    s_circumflex_accent_k = np.random.randn(dim_0, n, )
    return l, dist, σ, bᵢ, b, s_circumflex_accent_k


if __name__ == '__main__':
    l, dist, σ, bᵢ, b, s_circumflex_accent_k = generateRandomData()
    print("l:", l)
    print("dist:", dist)
    print("σ:", σ)
    print("bᵢ:", bᵢ)
    print("b:", b)
    print("s_circumflex_accent_k:", s_circumflex_accent_k)
    func_value = atlas_refinement_3(l, dist, σ, bᵢ, b, s_circumflex_accent_k)
    print("return value: ", func_value.G_σ_left_parenthesis_s_circumflex_accent_k_i_right_parenthesis)