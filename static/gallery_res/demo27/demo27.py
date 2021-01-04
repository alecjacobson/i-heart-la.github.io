"""
`Ω` = [`e_1` `e_2`][`k_1` 0
		     0    `k_2`] [`e_1`^T
				  `e_2`^T]

where

`k_1`: ℝ  : control the desired kernel variance in either edge or orthogonal direction
`k_2`: ℝ  : control the desired kernel variance in either edge or orthogonal direction
`e_1`: ℝ ^ 3: orthogonal direction vectors
`e_2`: ℝ ^ 3: orthogonal direction vectors
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo27(k_1, k_2, e_1, e_2):
    """
    :param :k_1 : ℝ  : control the desired kernel variance in either edge or orthogonal direction
    :param :k_2 : ℝ  : control the desired kernel variance in either edge or orthogonal direction
    :param :e_1 : ℝ ^ 3: orthogonal direction vectors
    :param :e_2 : ℝ ^ 3: orthogonal direction vectors
    """
    e_1 = np.asarray(e_1, dtype=np.float64)
    e_2 = np.asarray(e_2, dtype=np.float64)

    assert np.ndim(k_1) == 0
    assert np.ndim(k_2) == 0
    assert e_1.shape == (3,)
    assert e_2.shape == (3,)

    _Ω_0 = np.hstack((e_1, e_2))
    _Ω_1 = np.zeros((2, 2))
    _Ω_1[0] = [k_1, 0]
    _Ω_1[1] = [0, k_2]
    _Ω_2 = np.vstack((e_1.T, e_2.T))
    Ω = _Ω_0 @ _Ω_1 @ _Ω_2

    return Ω


def generateRandomData():
    k_1 = np.random.randn()
    k_2 = np.random.randn()
    e_1 = np.random.randn(3)
    e_2 = np.random.randn(3)
    return k_1, k_2, e_1, e_2


if __name__ == '__main__':
    k_1, k_2, e_1, e_2 = generateRandomData()
    print("k_1:", k_1)
    print("k_2:", k_2)
    print("e_1:", e_1)
    print("e_2:", e_2)
    func_value = demo27(k_1, k_2, e_1, e_2)
    print("func_value: ", func_value)