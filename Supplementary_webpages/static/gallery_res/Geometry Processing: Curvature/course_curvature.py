"""
`H(p)` = 1/(2π)int_[0, 2π] `kₙ`(φ, p) ∂φ

where 

p ∈ ℝ^3 : point on the surface
`kₙ`: ℝ,ℝ^3 → ℝ : normal curvature
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class course_curvatureResultType:
    def __init__( self, H_left_parenthesis_p_right_parenthesis):
        self.H_left_parenthesis_p_right_parenthesis = H_left_parenthesis_p_right_parenthesis


def course_curvature(p, kₙ):
    """
    :param :p : point on the surface
    :param :kₙ : ℝ,ℝ^3 → ℝ : normal curvature
    """
    p = np.asarray(p, dtype=np.float64)

    assert p.shape == (3,)

    H_left_parenthesis_p_right_parenthesis = 1 / (2 * np.pi) * quad(lambda φ: kₙ(φ, p), 0, 2 * np.pi)[0]
    return course_curvatureResultType(H_left_parenthesis_p_right_parenthesis)


def generateRandomData():
    p = np.random.randn(3)
    def kₙ(p0, p1):
        return np.random.randn()
    return p, kₙ


if __name__ == '__main__':
    p, kₙ = generateRandomData()
    print("p:", p)
    print("kₙ:", kₙ)
    func_value = course_curvature(p, kₙ)
    print("return value: ", func_value.H_left_parenthesis_p_right_parenthesis)