"""
from linearalgebra: tr

`J₃` = 1₃,₃
`κ_angle(Dₘ)` = 3(√2 v)^(2/3)(7/4‖`Dₘ`‖_F^2-1/4tr(`J₃``Dₘ`ᵀ`Dₘ`))⁻¹

where

`Dₘ` ∈ ℝ^(3×3)  
v ∈ ℝ
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class anisotropic_elasticity_47ResultType:
    def __init__( self, J3, κ_angle_left_parenthesis_Dₘ_right_parenthesis):
        self.J3 = J3
        self.κ_angle_left_parenthesis_Dₘ_right_parenthesis = κ_angle_left_parenthesis_Dₘ_right_parenthesis


def anisotropic_elasticity_47(Dₘ, v):
    Dₘ = np.asarray(Dₘ, dtype=np.float64)

    assert Dₘ.shape == (3, 3)
    assert np.ndim(v) == 0

    J3 = np.ones((3, 3))
    κ_angle_left_parenthesis_Dₘ_right_parenthesis = 3 * np.power((np.sqrt(2) * v), (2 / 3)) * 1 / ((7 / 4 * np.power(np.linalg.norm(Dₘ, 'fro'), 2) - 1 / 4 * np.trace(J3 @ Dₘ.T @ Dₘ)))
    return anisotropic_elasticity_47ResultType(J3, κ_angle_left_parenthesis_Dₘ_right_parenthesis)


def generateRandomData():
    v = np.random.randn()
    Dₘ = np.random.randn(3, 3)
    return Dₘ, v


if __name__ == '__main__':
    Dₘ, v = generateRandomData()
    print("Dₘ:", Dₘ)
    print("v:", v)
    func_value = anisotropic_elasticity_47(Dₘ, v)
    print("return value: ", func_value.κ_angle_left_parenthesis_Dₘ_right_parenthesis)