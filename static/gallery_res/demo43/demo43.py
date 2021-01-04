"""
`E(θ,σ²)` = -sum_n log(sum_m 1/||y|| 1/(2σ²)exp(-||x_n - y_m||²/(2σ²)) + 1/||x||)

where


x_i: ℝ^2
y_j: ℝ^2 
σ: ℝ
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo43(x, y, σ):
    """
    :param :x : ℝ^2
    :param :y : ℝ^2
    :param :σ : ℝ
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)

    _dim_0 = x.shape[0]
    _dim_1 = y.shape[0]
    assert x.shape == (_dim_0, 2, )
    assert y.shape == (_dim_1, 2, )
    assert np.ndim(σ) == 0

    _sum_0 = 0
    for n in range(1, len(x)+1):
        _sum_0 += np.log10(_sum_1 + 1 / )
    E_left_parenthesis_θ_comma_σ²_right_parenthesis = -_sum_0

    return E_left_parenthesis_θ_comma_σ²_right_parenthesis


def generateRandomData():
    σ = np.random.randn()
    _dim_0 = np.random.randint(10)
    _dim_1 = np.random.randint(10)
    x = np.random.randn(_dim_0, 2, )
    y = np.random.randn(_dim_1, 2, )
    return x, y, σ


if __name__ == '__main__':
    x, y, σ = generateRandomData()
    print("x:", x)
    print("y:", y)
    print("σ:", σ)
    func_value = demo43(x, y, σ)
    print("func_value: ", func_value)