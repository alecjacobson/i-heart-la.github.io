"""
`∂²I₅/∂f²` = 2[A₁,₁I₃  A₁,₂I₃  A₁,₃I₃
               A₂,₁I₃  A₂,₂I₃  A₂,₃I₃
               A₃,₁I₃  A₃,₂I₃  A₃,₃I₃] 

where

A ∈ ℝ^(3×3)
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo17(A):
    A = np.asarray(A, dtype=np.float64)

    assert A.shape == (3, 3)

    __partial_differential_2I5_soliduspartial_differential_f2_0 = np.block([[A[1-1, 1-1] * np.identity(3), A[1-1, 2-1] * np.identity(3), A[1-1, 3-1] * np.identity(3)], [A[2-1, 1-1] * np.identity(3), A[2-1, 2-1] * np.identity(3), A[2-1, 3-1] * np.identity(3)], [A[3-1, 1-1] * np.identity(3), A[3-1, 2-1] * np.identity(3), A[3-1, 3-1] * np.identity(3)]])
    _partial_differential_2I5_soliduspartial_differential_f2 = 2 * __partial_differential_2I5_soliduspartial_differential_f2_0

    return _partial_differential_2I5_soliduspartial_differential_f2


def generateRandomData():
    A = np.random.randn(3, 3)
    return A


if __name__ == '__main__':
    A = generateRandomData()
    print("A:", A)
    func_value = demo17(A)
    print("func_value: ", func_value)