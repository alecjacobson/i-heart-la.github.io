"""
`T_1` = 1/sqrt(2)U V^T
c = [0 0 0
                   0 0 -1
                   0 1 0]

where

U: ℝ^(3×3) 
V: ℝ^(3×3)
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo22(U, V):
    """
    :param :U : ℝ^(3×3)
    :param :V : ℝ^(3×3)
    """
    U = np.asarray(U, dtype=np.float64)
    V = np.asarray(V, dtype=np.float64)

    assert U.shape == (3, 3)
    assert V.shape == (3, 3)

    T_1 = 1 / np.sqrt(2) * U @ V.T

    _c_0 = np.zeros((3, 3))
    _c_0[0] = [0, 0, 0]
    _c_0[1] = [0, 0, -1]
    _c_0[2] = [0, 1, 0]
    c = _c_0

    return c


def generateRandomData():
    U = np.random.randn(3, 3)
    V = np.random.randn(3, 3)
    return U, V


if __name__ == '__main__':
    U, V = generateRandomData()
    print("U:", U)
    print("V:", V)
    func_value = demo22(U, V)
    print("func_value: ", func_value)