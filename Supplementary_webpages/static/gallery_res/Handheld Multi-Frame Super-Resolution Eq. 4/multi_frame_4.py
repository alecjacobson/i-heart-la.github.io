"""
Ω = [`e₁` `e₂`][`k₁`   0
                 0    `k₂`] [`e₁`ᵀ
                             `e₂`ᵀ]

where

`k₁` ∈ ℝ  
`k₂` ∈ ℝ 
`e₁` ∈ ℝ ^ 3
`e₂` ∈ ℝ ^ 3
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class multi_frame_4ResultType:
    def __init__( self, Ω):
        self.Ω = Ω


def multi_frame_4(k1, k2, e1, e2):
    e1 = np.asarray(e1, dtype=np.float64)
    e2 = np.asarray(e2, dtype=np.float64)

    assert np.ndim(k1) == 0
    assert np.ndim(k2) == 0
    assert e1.shape == (3,)
    assert e2.shape == (3,)

    Ω_0 = np.hstack(((e1).reshape(3, 1), (e2).reshape(3, 1)))
    Ω_1 = np.zeros((2, 2))
    Ω_1[0] = [k1, 0]
    Ω_1[1] = [0, k2]
    Ω_2 = np.vstack((e1.T.reshape(1, 3), e2.T.reshape(1, 3)))
    Ω = Ω_0 @ Ω_1 @ Ω_2
    return multi_frame_4ResultType(Ω)


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
    print("return value: ", func_value.Ω)