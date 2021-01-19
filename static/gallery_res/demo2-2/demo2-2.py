"""
`n(v)` = (∑_(i for i ∈ `N₁(v)`) α_i n(T_i))/||∑_(i for i ∈ `N₁(v)`) α_i n(T_i)|| 

where
 
T_i: ℝ^(3×3)
α_i: ℝ
`N₁(v)`: {ℤ}
n: ℝ^(3×3) -> ℝ^3
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo2-2(T, α, N1_left_parenthesis_v_right_parenthesis, n):
    """
    :param :T : ℝ^(3×3)
    :param :α : ℝ
    :param :N1_left_parenthesis_v_right_parenthesis : {ℤ}
    :param :n : ℝ^(3×3) -> ℝ^3
    """
    T = np.asarray(T, dtype=np.float64)
    α = np.asarray(α)

    _dim_0 = α.shape[0]
    assert T.shape == (_dim_0, 3, 3)
    assert α.shape == (_dim_0,)
    assert isinstance(N1_left_parenthesis_v_right_parenthesis, list) and len(N1_left_parenthesis_v_right_parenthesis) > 0

    _sum_0 = np.zeros((3, ))
    for i in range(1, len(α)+1):
        if((i) in N1_left_parenthesis_v_right_parenthesis):
            _sum_0 += α[i-1] * n(T[i-1])
    _sum_1 = np.zeros((3, ))
    for i in range(1, len(α)+1):
        if((i) in N1_left_parenthesis_v_right_parenthesis):
            _sum_1 += α[i-1] * n(T[i-1])
    n_left_parenthesis_v_right_parenthesis = (_sum_0) / np.linalg.norm(_sum_1, 2)

    return n_left_parenthesis_v_right_parenthesis


def generateRandomData():
    _dim_0 = np.random.randint(10)
    T = np.random.randn(_dim_0, 3, 3)
    α = np.random.randn(_dim_0)
    N1_left_parenthesis_v_right_parenthesis = []
    N1_left_parenthesis_v_right_parenthesis_0 = np.random.randint(1, 10)
    for i in range(N1_left_parenthesis_v_right_parenthesis_0):
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
    func_value = demo2-2(T, α, N1_left_parenthesis_v_right_parenthesis, n)
    print("func_value: ", func_value)