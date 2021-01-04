"""
[A⁻¹+A⁻¹BS⁻¹BᵀA⁻¹   -A⁻¹BS⁻¹
    -S⁻¹BᵀA⁻¹           S⁻¹]

where

A: ℝ^(4×4) 
B: ℝ^(4×4) 
S: ℝ^(4×4) 
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo14(A, B, S):
    """
    :param :A : ℝ^(4×4)
    :param :B : ℝ^(4×4)
    :param :S : ℝ^(4×4)
    """
    A = np.asarray(A, dtype=np.float64)
    B = np.asarray(B, dtype=np.float64)
    S = np.asarray(S, dtype=np.float64)

    assert A.shape == (4, 4)
    assert B.shape == (4, 4)
    assert S.shape == (4, 4)

    _ret_0 = np.block([[np.linalg.inv(A) + np.linalg.inv(A) @ B @ np.linalg.inv(S) @ B.T @ np.linalg.inv(A), -np.linalg.inv(A) @ B @ np.linalg.inv(S)], [-np.linalg.inv(S) @ B.T @ np.linalg.inv(A), np.linalg.inv(S)]])
    ret = _ret_0
    return ret


def generateRandomData():
    A = np.random.randn(4, 4)
    B = np.random.randn(4, 4)
    S = np.random.randn(4, 4)
    return A, B, S


if __name__ == '__main__':
    A, B, S = generateRandomData()
    print("A:", A)
    print("B:", B)
    print("S:", S)
    func_value = demo14(A, B, S)
    print("func_value: ", func_value)