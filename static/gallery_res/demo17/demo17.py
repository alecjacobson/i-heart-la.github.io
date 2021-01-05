"""
`∂²I₅/∂f²` = 2[(A_1,1)I_3  (A_1,2)I_3  (A_1,3)I_3
               (A_2,1)I_3  (A_2,2)I_3  (A_2,3)I_3
               (A_3,1)I_3  (A_3,2)I_3  (A_3,3)I_3] 

where

A: ℝ^(3×3) 
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo17(A):
    """
    :param :A : ℝ^(3×3)
    """
    A = np.asarray(A, dtype=np.float64)

    assert A.shape == (3, 3)

    __partial_differential_²I₅_soliduspartial_differential_f²_0 = np.block([[(A[1-1, 1-1]) * np.identity(3), (A[1-1, 2-1]) * np.identity(3), (A[1-1, 3-1]) * np.identity(3)], [(A[2-1, 1-1]) * np.identity(3), (A[2-1, 2-1]) * np.identity(3), (A[2-1, 3-1]) * np.identity(3)], [(A[3-1, 1-1]) * np.identity(3), (A[3-1, 2-1]) * np.identity(3), (A[3-1, 3-1]) * np.identity(3)]])
    _partial_differential_²I₅_soliduspartial_differential_f² = 2 * __partial_differential_²I₅_soliduspartial_differential_f²_0

    return _partial_differential_²I₅_soliduspartial_differential_f²


def generateRandomData():
    A = np.random.randn(3, 3)
    return A


if __name__ == '__main__':
    A = generateRandomData()
    print("A:", A)
    func_value = demo17(A)
    print("func_value: ", func_value)