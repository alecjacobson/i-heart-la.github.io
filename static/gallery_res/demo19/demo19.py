"""
from linearalgebra: tr

`k_angle(Dₘ)` = 3(sqrt(2)v)^(2/3)(7/4||`Dₘ`||_F^2-1/4tr(J_3 `Dₘ`ᵀ`Dₘ`))⁻¹

where

`Dₘ`: ℝ^(n×n) 
J_i: ℝ^(n×n) 
v: ℝ
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo19(Dₘ, J, v):
    """
    :param :Dₘ : ℝ^(n×n)
    :param :J : ℝ^(n×n)
    :param :v : ℝ
    """
    Dₘ = np.asarray(Dₘ, dtype=np.float64)
    J = np.asarray(J, dtype=np.float64)

    n = J.shape[2]
    _dim_0 = J.shape[0]
    assert Dₘ.shape == (n, n)
    assert J.shape == (_dim_0, n, n)
    assert np.ndim(v) == 0

    k_angle_left_parenthesis_Dₘ_right_parenthesis = 3 * np.power((np.sqrt(2) * v), (2 / 3)) * 1 / ((7 / 4 * np.power(np.linalg.norm(Dₘ, 'fro'), 2) - 1 / 4 * np.trace(J[3-1] @ Dₘ.T @ Dₘ)))

    return k_angle_left_parenthesis_Dₘ_right_parenthesis


def generateRandomData():
    v = np.random.randn()
    n = np.random.randint(10)
    _dim_0 = np.random.randint(10)
    Dₘ = np.random.randn(n, n)
    J = np.random.randn(_dim_0, n, n)
    return Dₘ, J, v


if __name__ == '__main__':
    Dₘ, J, v = generateRandomData()
    print("Dₘ:", Dₘ)
    print("J:", J)
    print("v:", v)
    func_value = demo19(Dₘ, J, v)
    print("func_value: ", func_value)