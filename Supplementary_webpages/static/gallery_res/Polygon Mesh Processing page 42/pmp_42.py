"""
`n(v)` = (∑_(i for i ∈ `N₁(v)`) α_i n(T_i))/‖∑_(i for i ∈ `N₁(v)`) α_i n(T_i)‖

where
 
T_i ∈ ℝ^(3×3)
α_i ∈ ℝ
`N₁(v)` ∈ {ℤ}
n: ℝ^(3×3) → ℝ^3
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class pmp_42ResultType:
    def __init__( self, n_left_parenthesis_v_right_parenthesis):
        self.n_left_parenthesis_v_right_parenthesis = n_left_parenthesis_v_right_parenthesis


def pmp_42(T, α, N1_left_parenthesis_v_right_parenthesis, n):
    """
    :param :n : ℝ^(3×3) → ℝ^3
    """
    T = np.asarray(T, dtype=np.float64)
    α = np.asarray(α, dtype=np.float64)

    dim_0 = α.shape[0]
    assert T.shape == (dim_0, 3, 3)
    assert α.shape == (dim_0,)
    assert isinstance(N1_left_parenthesis_v_right_parenthesis, list) and len(N1_left_parenthesis_v_right_parenthesis) > 0

    sum_0 = np.zeros((3, ))
    for i in range(1, len(α)+1):
        if((i) in N1_left_parenthesis_v_right_parenthesis):
            sum_0 += α[i-1] * n(T[i-1])
    sum_1 = np.zeros((3, ))
    for i in range(1, len(α)+1):
        if((i) in N1_left_parenthesis_v_right_parenthesis):
            sum_1 += α[i-1] * n(T[i-1])
    n_left_parenthesis_v_right_parenthesis = (sum_0) / np.linalg.norm(sum_1, 2)
    return pmp_42ResultType(n_left_parenthesis_v_right_parenthesis)


def generateRandomData():
    dim_0 = np.random.randint(10)
    T = np.random.randn(dim_0, 3, 3)
    α = np.random.randn(dim_0)
    N1_left_parenthesis_v_right_parenthesis = []
    dim_1 = np.random.randint(1, 10)
    for i in range(dim_1):
        N1_left_parenthesis_v_right_parenthesis.append((np.random.randint(10)))
    def n(p0):
        return np.random.randn(3)
    return T, α, N1_left_parenthesis_v_right_parenthesis, n


if __name__ == '__main__':
    T, α, N1_left_parenthesis_v_right_parenthesis, n = generateRandomData()
    print("T:", T)
    print("α:", α)
    print("N1_left_parenthesis_v_right_parenthesis:", N1_left_parenthesis_v_right_parenthesis)
    print("n:", n)
    func_value = pmp_42(T, α, N1_left_parenthesis_v_right_parenthesis, n)
    print("return value: ", func_value.n_left_parenthesis_v_right_parenthesis)