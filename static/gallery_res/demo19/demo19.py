"""
from linearalgebra: tr

`J₃` = [1_3,3]
`k_angle(Dₘ)` = 3(sqrt(2)v)^(2/3)(7/4||Dₘ||_F^2-1/4tr(`J₃` (Dₘ)ᵀDₘ))⁻¹

where

Dₘ: ℝ^(3×3)  
v: ℝ
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo19(Dₘ, v):
    """
    :param :Dₘ : ℝ^(3×3)
    :param :v : ℝ
    """
    Dₘ = np.asarray(Dₘ, dtype=np.float64)

    assert Dₘ.shape == (3, 3)
    assert np.ndim(v) == 0

    _J3_0 = np.vstack((np.ones((3, 3))))
    J3 = _J3_0

    k_angle_left_parenthesis_Dₘ_right_parenthesis = 3 * np.power((np.sqrt(2) * v), (2 / 3)) * 1 / ((7 / 4 * np.power(np.linalg.norm(Dₘ, 'fro'), 2) - 1 / 4 * np.trace(J3 @ (Dₘ).T @ Dₘ)))

    return k_angle_left_parenthesis_Dₘ_right_parenthesis


def generateRandomData():
    v = np.random.randn()
    Dₘ = np.random.randn(3, 3)
    return Dₘ, v


if __name__ == '__main__':
    Dₘ, v = generateRandomData()
    print("Dₘ:", Dₘ)
    print("v:", v)
    func_value = demo19(Dₘ, v)
    print("func_value: ", func_value)