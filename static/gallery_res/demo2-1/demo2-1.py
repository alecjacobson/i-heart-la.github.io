"""
`xᵢ` = T_*,1
`xⱼ` = T_*,2
`xₖ` = T_*,3
`n(T)` = (`xⱼ`-`xᵢ`)×(`xₖ`-`xᵢ`)/||(`xⱼ`-`xᵢ`)×(`xₖ`-`xᵢ`)||

where
 
T: ℝ^(3×3)
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo2-1(T):
    """
    :param :T : ℝ^(3×3)
    """
    T = np.asarray(T, dtype=np.float64)

    assert T.shape == (3, 3)

    xᵢ = T[:, 1-1]

    xⱼ = T[:, 2-1]

    xₖ = T[:, 3-1]

    n_left_parenthesis_T_right_parenthesis = np.cross((xⱼ - xᵢ), (xₖ - xᵢ)) / np.linalg.norm(np.cross((xⱼ - xᵢ), (xₖ - xᵢ)), 2)

    return n_left_parenthesis_T_right_parenthesis


def generateRandomData():
    T = np.random.randn(3, 3)
    return T


if __name__ == '__main__':
    T = generateRandomData()
    print("T:", T)
    func_value = demo2-1(T)
    print("func_value: ", func_value)