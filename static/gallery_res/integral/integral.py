"""
`H(p)` = 1/(2π)int_[0, 2π] `kₙ`(φ, p) ∂φ

where 

p ∈ ℝ^3 : point on the surface
`kₙ`: ℝ,ℝ^3->ℝ : normal curvature
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def integral(p, kₙ):
    """
    :param :p : point on the surface
    :param :kₙ : ℝ,ℝ^3->ℝ : normal curvature
    """
    p = np.asarray(p, dtype=np.float64)

    assert p.shape == (3,)

    H_left_parenthesis_p_right_parenthesis = 1 / (2 * np.pi) * quad(lambda φ: kₙ(φ, p), 0, 2 * np.pi)[0]

    return H_left_parenthesis_p_right_parenthesis


def generateRandomData():
    p = np.random.randn(3)
    def kₙ(p0, p1):
        return np.random.randn()
    return p, kₙ


if __name__ == '__main__':
    p, kₙ = generateRandomData()
    print("p:", p)
    print("kₙ:", kₙ)
    func_value = integral(p, kₙ)
    print("func_value: ", func_value)