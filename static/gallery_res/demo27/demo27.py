"""
Ω = [`e₁` `e₂`][`k₁`   0
                 0    `k₂`] [`e₁`ᵀ
			                 `e₂`ᵀ]

where

`k₁`: ℝ  : control the desired kernel variance in either edge or orthogonal direction
`k₂`: ℝ  : control the desired kernel variance in either edge or orthogonal direction
`e₁`: ℝ ^ 3: orthogonal direction vectors
`e₂`: ℝ ^ 3: orthogonal direction vectors
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo27(k₁, k₂, e₁, e₂):
    """
    :param :k₁ : ℝ  : control the desired kernel variance in either edge or orthogonal direction
    :param :k₂ : ℝ  : control the desired kernel variance in either edge or orthogonal direction
    :param :e₁ : ℝ ^ 3: orthogonal direction vectors
    :param :e₂ : ℝ ^ 3: orthogonal direction vectors
    """
    e₁ = np.asarray(e₁, dtype=np.float64)
    e₂ = np.asarray(e₂, dtype=np.float64)

    assert np.ndim(k₁) == 0
    assert np.ndim(k₂) == 0
    assert e₁.shape == (3,)
    assert e₂.shape == (3,)

    _Ω_0 = np.hstack(((e₁).reshape(3, 1), (e₂).reshape(3, 1)))
    _Ω_1 = np.zeros((2, 2))
    _Ω_1[0] = [k₁, 0]
    _Ω_1[1] = [0, k₂]
    _Ω_2 = np.vstack((e₁.T.reshape(1, 3), e₂.T.reshape(1, 3)))
    Ω = _Ω_0 @ _Ω_1 @ _Ω_2

    return Ω


def generateRandomData():
    k₁ = np.random.randn()
    k₂ = np.random.randn()
    e₁ = np.random.randn(3)
    e₂ = np.random.randn(3)
    return k₁, k₂, e₁, e₂


if __name__ == '__main__':
    k₁, k₂, e₁, e₂ = generateRandomData()
    print("k₁:", k₁)
    print("k₂:", k₂)
    print("e₁:", e₁)
    print("e₂:", e₂)
    func_value = demo27(k₁, k₂, e₁, e₂)
    print("func_value: ", func_value)