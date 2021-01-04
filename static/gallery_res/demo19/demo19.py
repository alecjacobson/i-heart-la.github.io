"""
from linearalgebra: tr

`k_angle(D_m)` = 3(sqrt(2)v)^(2/3)(7/4||`D_m`||_F^2 -1/4tr(J_3 `D_m`^T`D_m`))^(-1)

where

`D_m`: ℝ^(n×n) 
J_i: ℝ^(n×n) 
v: ℝ
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo19(D_m, J, v):
    """
    :param :D_m : ℝ^(n×n)
    :param :J : ℝ^(n×n)
    :param :v : ℝ
    """
    D_m = np.asarray(D_m, dtype=np.float64)
    J = np.asarray(J, dtype=np.float64)

    n = J.shape[2]
    _dim_0 = J.shape[0]
    assert D_m.shape == (n, n)
    assert J.shape == (_dim_0, n, n)
    assert np.ndim(v) == 0

    k_angle_left_parenthesis_D_m_right_parenthesis = 3 * np.power((np.sqrt(2) * v), (2 / 3)) * 1 / ((7 / 4 * np.power(np.linalg.norm(D_m, 'fro'), 2) - 1 / 4 * np.trace(J[3-1] @ D_m.T @ D_m)))

    return k_angle_left_parenthesis_D_m_right_parenthesis


def generateRandomData():
    v = np.random.randn()
    n = np.random.randint(10)
    _dim_0 = np.random.randint(10)
    D_m = np.random.randn(n, n)
    J = np.random.randn(_dim_0, n, n)
    return D_m, J, v


if __name__ == '__main__':
    D_m, J, v = generateRandomData()
    print("D_m:", D_m)
    print("J:", J)
    print("v:", v)
    func_value = demo19(D_m, J, v)
    print("func_value: ", func_value)