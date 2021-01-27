"""
Ω = [`e₁` `e₂`][`k₁`   0
                 0    `k₂`] [`e₁`ᵀ
                             `e₂`ᵀ]

where

`k₁` ∈ ℝ  : control the desired kernel variance in either edge or orthogonal direction
`k₂` ∈ ℝ  : control the desired kernel variance in either edge or orthogonal direction
`e₁` ∈ ℝ ^ 3: orthogonal direction vectors
`e₂` ∈ ℝ ^ 3: orthogonal direction vectors
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def multi_frame_4(k1, k2, e1, e2):
    """
    :param :k1 : control the desired kernel variance in either edge or orthogonal direction
    :param :k2 : control the desired kernel variance in either edge or orthogonal direction
    :param :e1 : orthogonal direction vectors
    :param :e2 : orthogonal direction vectors
    """
    e1 = np.asarray(e1, dtype=np.float64)
    e2 = np.asarray(e2, dtype=np.float64)

    assert np.ndim(k1) == 0
    assert np.ndim(k2) == 0
    assert e1.shape == (3,)
    assert e2.shape == (3,)

    _Ω_0 = np.hstack(((e1).reshape(3, 1), (e2).reshape(3, 1)))
    _Ω_1 = np.zeros((2, 2))
    _Ω_1[0] = [k1, 0]
    _Ω_1[1] = [0, k2]
    _Ω_2 = np.vstack((e1.T.reshape(1, 3), e2.T.reshape(1, 3)))
    Ω = _Ω_0 @ _Ω_1 @ _Ω_2

    return Ω


def generateRandomData():
    k1 = np.random.randn()
    k2 = np.random.randn()
    e1 = np.random.randn(3)
    e2 = np.random.randn(3)
    return k1, k2, e1, e2


if __name__ == '__main__':
    k1, k2, e1, e2 = generateRandomData()
    print("k1:", k1)
    print("k2:", k2)
    print("e1:", e1)
    print("e2:", e2)
    func_value = multi_frame_4(k1, k2, e1, e2)
    print("func_value: ", func_value)