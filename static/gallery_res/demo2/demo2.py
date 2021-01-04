"""
n_i = (T_i,*,2-T_i,*,1)×(T_i,*,3-T_i,*,1)/||(T_i,*,2-T_i,*,1)×(T_i,*,3-T_i,*,1)||

`n(v)` = (∑_(i for i ∈ `N1(v)`) α_i n_i)/||∑_(i for i ∈ `N1(v)`) α_i n_i||
where
 
T_i: ℝ^(3×3)
α_i: ℝ
`N1(v)`: {ℤ}
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo2(T, α, N1_left_parenthesis_v_right_parenthesis):
    """
    :param :T : ℝ^(3×3)
    :param :α : ℝ
    :param :N1_left_parenthesis_v_right_parenthesis : {ℤ}
    """
    T = np.asarray(T, dtype=np.float64)
    α = np.asarray(α)

    _dim_0 = α.shape[0]
    assert T.shape == (_dim_0, 3, 3)
    assert α.shape == (_dim_0,)
    assert isinstance(N1_left_parenthesis_v_right_parenthesis, list) and len(N1_left_parenthesis_v_right_parenthesis) > 0

    n = np.zeros((_dim_0, 3, ))
    for i in range(1, _dim_0+1):
        n[i-1] = np.cross((T[i-1][:, 2-1] - T[i-1][:, 1-1]), (T[i-1][:, 3-1] - T[i-1][:, 1-1])) / np.linalg.norm(np.cross((T[i-1][:, 2-1] - T[i-1][:, 1-1]), (T[i-1][:, 3-1] - T[i-1][:, 1-1])), 2)

    _sum_0 = np.zeros((3, ))
    for i in range(1, len(n)+1):
        if((i) in N1_left_parenthesis_v_right_parenthesis):
            _sum_0 += α[i-1] * n[i-1]
    _sum_1 = np.zeros((3, ))
    for i in range(1, len(n)+1):
        if((i) in N1_left_parenthesis_v_right_parenthesis):
            _sum_1 += α[i-1] * n[i-1]
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
    return T, α, N1_left_parenthesis_v_right_parenthesis


if __name__ == '__main__':
    T, α, N1_left_parenthesis_v_right_parenthesis = generateRandomData()
    print("T:", T)
    print("α:", α)
    print("N1_left_parenthesis_v_right_parenthesis:", N1_left_parenthesis_v_right_parenthesis)
    func_value = demo2(T, α, N1_left_parenthesis_v_right_parenthesis)
    print("func_value: ", func_value)