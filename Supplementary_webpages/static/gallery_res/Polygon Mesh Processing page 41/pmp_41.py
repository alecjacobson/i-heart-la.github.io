"""
`xᵢ` = T_*,1
`xⱼ` = T_*,2
`xₖ` = T_*,3
`n(T)` = (`xⱼ`-`xᵢ`)×(`xₖ`-`xᵢ`)/‖(`xⱼ`-`xᵢ`)×(`xₖ`-`xᵢ`)‖

where
 
T ∈ ℝ^(3×3)
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class pmp_41ResultType:
    def __init__( self, xᵢ, xⱼ, xₖ, n_left_parenthesis_T_right_parenthesis):
        self.xᵢ = xᵢ
        self.xⱼ = xⱼ
        self.xₖ = xₖ
        self.n_left_parenthesis_T_right_parenthesis = n_left_parenthesis_T_right_parenthesis


def pmp_41(T):
    T = np.asarray(T, dtype=np.float64)

    assert T.shape == (3, 3)

    xᵢ = T[:, 1-1]
    xⱼ = T[:, 2-1]
    xₖ = T[:, 3-1]
    n_left_parenthesis_T_right_parenthesis = np.cross((xⱼ - xᵢ), (xₖ - xᵢ)) / np.linalg.norm(np.cross((xⱼ - xᵢ), (xₖ - xᵢ)), 2)
    return pmp_41ResultType(xᵢ, xⱼ, xₖ, n_left_parenthesis_T_right_parenthesis)


def generateRandomData():
    T = np.random.randn(3, 3)
    return T


if __name__ == '__main__':
    T = generateRandomData()
    print("T:", T)
    func_value = pmp_41(T)
    print("return value: ", func_value.n_left_parenthesis_T_right_parenthesis)