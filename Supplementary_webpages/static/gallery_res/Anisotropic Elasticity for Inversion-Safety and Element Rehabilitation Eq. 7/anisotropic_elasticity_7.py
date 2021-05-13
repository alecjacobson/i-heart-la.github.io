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


class anisotropic_elasticity_7ResultType:
    def __init__( self, _partial_differential_2I5_soliduspartial_differential_f2):
        self._partial_differential_2I5_soliduspartial_differential_f2 = _partial_differential_2I5_soliduspartial_differential_f2


def anisotropic_elasticity_7(A):
    A = np.asarray(A, dtype=np.float64)

    assert A.shape == (3, 3)

    _partial_differential_2I5_soliduspartial_differential_f2_0 = np.block([[A[1-1, 1-1] * np.identity(3), A[1-1, 2-1] * np.identity(3), A[1-1, 3-1] * np.identity(3)], [A[2-1, 1-1] * np.identity(3), A[2-1, 2-1] * np.identity(3), A[2-1, 3-1] * np.identity(3)], [A[3-1, 1-1] * np.identity(3), A[3-1, 2-1] * np.identity(3), A[3-1, 3-1] * np.identity(3)]])
    _partial_differential_2I5_soliduspartial_differential_f2 = 2 * _partial_differential_2I5_soliduspartial_differential_f2_0
    return anisotropic_elasticity_7ResultType(_partial_differential_2I5_soliduspartial_differential_f2)


def generateRandomData():
    A = np.random.randn(3, 3)
    return A


if __name__ == '__main__':
    A = generateRandomData()
    print("A:", A)
    func_value = anisotropic_elasticity_7(A)
    print("return value: ", func_value._partial_differential_2I5_soliduspartial_differential_f2)