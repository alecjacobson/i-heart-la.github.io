"""
`C(x,y)` = (sum_n sum_i c_n,i w_n,i R_n) / (sum_n sum_i w_n,i R_n)

where

c: ℝ^(x×y): the value of the Bayer pixel
w: ℝ^(x×y): the local sample weight
R: ℝ^x: the local robustness
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo26(c, w, R):
    """
    :param :c : ℝ^(x×y): the value of the Bayer pixel
    :param :w : ℝ^(x×y): the local sample weight
    :param :R : ℝ^x: the local robustness
    """
    c = np.asarray(c, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    R = np.asarray(R, dtype=np.float64)

    x = R.shape[0]
    y = w.shape[1]
    assert c.shape == (x, y)
    assert w.shape == (x, y)
    assert R.shape == (x,)

    _sum_0 = 0
    for n in range(1, len(R)+1):
        _sum_1 = 0
        for i in range(1, len(w)+1):
            _sum_1 += c[n-1, i-1] * w[n-1, i-1] * R[n-1]
        _sum_0 += _sum_1
    _sum_2 = 0
    for n in range(1, len(R)+1):
        _sum_3 = 0
        for i in range(1, len(w)+1):
            _sum_3 += w[n-1, i-1] * R[n-1]
        _sum_2 += _sum_3
    C_left_parenthesis_x_comma_y_right_parenthesis = (_sum_0) / (_sum_2)

    return C_left_parenthesis_x_comma_y_right_parenthesis


def generateRandomData():
    x = np.random.randint(10)
    y = np.random.randint(10)
    c = np.random.randn(x, y)
    w = np.random.randn(x, y)
    R = np.random.randn(x)
    return c, w, R


if __name__ == '__main__':
    c, w, R = generateRandomData()
    print("c:", c)
    print("w:", w)
    print("R:", R)
    func_value = demo26(c, w, R)
    print("func_value: ", func_value)