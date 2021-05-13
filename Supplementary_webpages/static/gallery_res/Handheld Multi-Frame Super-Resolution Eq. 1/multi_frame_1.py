"""
`C(x,y)` = (∑_n ∑_i c_n,i w_n,i R̂_n) / (∑_n ∑_i w_n,i R̂_n)

where

c ∈ ℝ^(f×s): the value of the Bayer pixel
w ∈ ℝ^(f×s): the local sample weight
R̂ ∈ ℝ^f: the local robustness
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class multi_frame_1ResultType:
    def __init__( self, C_left_parenthesis_x_comma_y_right_parenthesis):
        self.C_left_parenthesis_x_comma_y_right_parenthesis = C_left_parenthesis_x_comma_y_right_parenthesis


def multi_frame_1(c, w, R̂):
    """
    :param :c : the value of the Bayer pixel
    :param :w : the local sample weight
    :param :R̂ : the local robustness
    """
    c = np.asarray(c, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    R̂ = np.asarray(R̂, dtype=np.float64)

    f = c.shape[0]
    s = c.shape[1]
    assert c.shape == (f, s)
    assert w.shape == (f, s)
    assert R̂.shape == (f,)

    sum_0 = 0
    for n in range(1, len(w)+1):
        sum_1 = 0
        for i in range(1, len(w)+1):
            sum_1 += c[n-1, i-1] * w[n-1, i-1] * R̂[n-1]
        sum_0 += sum_1
    sum_2 = 0
    for n in range(1, len(w)+1):
        sum_3 = 0
        for i in range(1, len(w)+1):
            sum_3 += w[n-1, i-1] * R̂[n-1]
        sum_2 += sum_3
    C_left_parenthesis_x_comma_y_right_parenthesis = (sum_0) / (sum_2)
    return multi_frame_1ResultType(C_left_parenthesis_x_comma_y_right_parenthesis)


def generateRandomData():
    f = np.random.randint(10)
    s = np.random.randint(10)
    c = np.random.randn(f, s)
    w = np.random.randn(f, s)
    R̂ = np.random.randn(f)
    return c, w, R̂


if __name__ == '__main__':
    c, w, R̂ = generateRandomData()
    print("c:", c)
    print("w:", w)
    print("R̂:", R̂)
    func_value = multi_frame_1(c, w, R̂)
    print("return value: ", func_value.C_left_parenthesis_x_comma_y_right_parenthesis)