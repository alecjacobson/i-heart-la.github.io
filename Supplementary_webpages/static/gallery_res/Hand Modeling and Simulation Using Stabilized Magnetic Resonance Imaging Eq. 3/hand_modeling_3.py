"""
min_(C ∈ ℝ^3) ∑_i ‖x_i + (R_i - I₃)C‖²

where

x_i ∈ ℝ^3
R_i ∈ ℝ^(3×3)
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class hand_modeling_3ResultType:
    def __init__( self, ret):
        self.ret = ret


def hand_modeling_3(x, R):
    x = np.asarray(x, dtype=np.float64)
    R = np.asarray(R, dtype=np.float64)

    dim_0 = x.shape[0]
    assert x.shape == (dim_0, 3, )
    assert R.shape == (dim_0, 3, 3)

    def target_0(C):
        sum_0 = 0
        for i in range(1, len(R)+1):
            sum_0 += np.power(np.linalg.norm(x[i-1] + (R[i-1] - np.identity(3)) @ C, 2), 2)
        return sum_0
    ret = minimize(target_0, np.zeros(3)).fun
    return hand_modeling_3ResultType(ret)


def generateRandomData():
    dim_0 = np.random.randint(10)
    x = np.random.randn(dim_0, 3, )
    R = np.random.randn(dim_0, 3, 3)
    return x, R


if __name__ == '__main__':
    x, R = generateRandomData()
    print("x:", x)
    print("R:", R)
    func_value = hand_modeling_3(x, R)
    print("return value: ", func_value.ret)