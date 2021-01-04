"""
[P_1  0
 0  P_3][  L              0
    P_3^TCP_2^TU^(-1)   -Lbar][U  L^(-1)P_1^TB
                               0     `U_bar`][P_2   0
                                               0   I_4]
where

P_i: ℝ^(4×4) 
B: ℝ^(4×4) 
C: ℝ^(4×4) 
L: ℝ^(4×4) 
Lbar: ℝ^(4×4) 
U: ℝ^(4×4) 
`U_bar`: ℝ^(4×4)
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo15(P, B, C, L, Lbar, U, U_bar):
    """
    :param :P : ℝ^(4×4)
    :param :B : ℝ^(4×4)
    :param :C : ℝ^(4×4)
    :param :L : ℝ^(4×4)
    :param :Lbar : ℝ^(4×4)
    :param :U : ℝ^(4×4)
    :param :U_bar : ℝ^(4×4)
    """
    P = np.asarray(P, dtype=np.float64)
    B = np.asarray(B, dtype=np.float64)
    C = np.asarray(C, dtype=np.float64)
    L = np.asarray(L, dtype=np.float64)
    Lbar = np.asarray(Lbar, dtype=np.float64)
    U = np.asarray(U, dtype=np.float64)
    U_bar = np.asarray(U_bar, dtype=np.float64)

    _dim_0 = P.shape[0]
    assert P.shape == (_dim_0, 4, 4)
    assert B.shape == (4, 4)
    assert C.shape == (4, 4)
    assert L.shape == (4, 4)
    assert Lbar.shape == (4, 4)
    assert U.shape == (4, 4)
    assert U_bar.shape == (4, 4)

    _ret_0 = np.block([[P[1-1], np.zeros((4, 4))], [np.zeros((4, 4)), P[3-1]]])
    _ret_1 = np.block([[L, np.zeros((4, 4))], [P[3-1].T @ C @ P[2-1].T @ np.linalg.inv(U), -Lbar]])
    _ret_2 = np.block([[U, np.linalg.inv(L) @ P[1-1].T @ B], [np.zeros((4, 4)), U_bar]])
    _ret_3 = np.block([[P[2-1], np.zeros((4, 4))], [np.zeros((4, 4)), np.identity(4)]])
    ret = _ret_0 @ _ret_1 @ _ret_2 @ _ret_3
    return ret


def generateRandomData():
    _dim_0 = np.random.randint(10)
    P = np.random.randn(_dim_0, 4, 4)
    B = np.random.randn(4, 4)
    C = np.random.randn(4, 4)
    L = np.random.randn(4, 4)
    Lbar = np.random.randn(4, 4)
    U = np.random.randn(4, 4)
    U_bar = np.random.randn(4, 4)
    return P, B, C, L, Lbar, U, U_bar


if __name__ == '__main__':
    P, B, C, L, Lbar, U, U_bar = generateRandomData()
    print("P:", P)
    print("B:", B)
    print("C:", C)
    print("L:", L)
    print("Lbar:", Lbar)
    print("U:", U)
    print("U_bar:", U_bar)
    func_value = demo15(P, B, C, L, Lbar, U, U_bar)
    print("func_value: ", func_value)