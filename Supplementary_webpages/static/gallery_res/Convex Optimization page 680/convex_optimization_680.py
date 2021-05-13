"""
[P₁  0
 0   P₃][    L        0
         P₃ᵀCP₂ᵀU⁻¹  -L̃][U  L⁻¹P₁ᵀB
                         0     Ũ   ][P₂   0
                                     0    I₄]

where

P_i ∈ ℝ^(4×4) 
B ∈ ℝ^(4×4) 
C ∈ ℝ^(4×4) 
L ∈ ℝ^(4×4) 
L̃ ∈ ℝ^(4×4) 
U ∈ ℝ^(4×4) 
Ũ ∈ ℝ^(4×4)
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class convex_optimization_680ResultType:
    def __init__( self, ret):
        self.ret = ret


def convex_optimization_680(P, B, C, L, L̃, U, Ũ):
    P = np.asarray(P, dtype=np.float64)
    B = np.asarray(B, dtype=np.float64)
    C = np.asarray(C, dtype=np.float64)
    L = np.asarray(L, dtype=np.float64)
    L̃ = np.asarray(L̃, dtype=np.float64)
    U = np.asarray(U, dtype=np.float64)
    Ũ = np.asarray(Ũ, dtype=np.float64)

    dim_0 = P.shape[0]
    assert P.shape == (dim_0, 4, 4)
    assert B.shape == (4, 4)
    assert C.shape == (4, 4)
    assert L.shape == (4, 4)
    assert L̃.shape == (4, 4)
    assert U.shape == (4, 4)
    assert Ũ.shape == (4, 4)

    ret_0 = np.block([[P[1-1], np.zeros((4, 4))], [np.zeros((4, 4)), P[3-1]]])
    ret_1 = np.block([[L, np.zeros((4, 4))], [P[3-1].T @ C @ P[2-1].T @ np.linalg.inv(U), -L̃]])
    ret_2 = np.block([[U, np.linalg.inv(L) @ P[1-1].T @ B], [np.zeros((4, 4)), Ũ]])
    ret_3 = np.block([[P[2-1], np.zeros((4, 4))], [np.zeros((4, 4)), np.identity(4)]])
    ret = ret_0 @ ret_1 @ ret_2 @ ret_3
    return convex_optimization_680ResultType(ret)


def generateRandomData():
    dim_0 = np.random.randint(10)
    P = np.random.randn(dim_0, 4, 4)
    B = np.random.randn(4, 4)
    C = np.random.randn(4, 4)
    L = np.random.randn(4, 4)
    L̃ = np.random.randn(4, 4)
    U = np.random.randn(4, 4)
    Ũ = np.random.randn(4, 4)
    return P, B, C, L, L̃, U, Ũ


if __name__ == '__main__':
    P, B, C, L, L̃, U, Ũ = generateRandomData()
    print("P:", P)
    print("B:", B)
    print("C:", C)
    print("L:", L)
    print("L̃:", L̃)
    print("U:", U)
    print("Ũ:", Ũ)
    func_value = convex_optimization_680(P, B, C, L, L̃, U, Ũ)
    print("return value: ", func_value.ret)