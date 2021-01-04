"""
min_(C ∈ ℝ^3) sum_i ||x_i + (R_i - I_3)C ||²

where

x_i: ℝ^3
R_i: ℝ^(3×3) 
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo38(x, R):
    """
    :param :x : ℝ^3
    :param :R : ℝ^(3×3)
    """
    x = np.asarray(x, dtype=np.float64)
    R = np.asarray(R, dtype=np.float64)

    _dim_0 = x.shape[0]
    assert x.shape == (_dim_0, 3, )
    assert R.shape == (_dim_0, 3, 3)

    def _target_0(C):
        _sum_0 = 0
        for i in range(1, len(x)+1):
            _sum_0 += np.power(np.linalg.norm(x[i-1] + (R[i-1] - np.identity(3)) @ C, 2), 2)
        return _sum_0
    ret = minimize(_target_0, np.zeros(3)).fun
    return ret


def generateRandomData():
    _dim_0 = np.random.randint(10)
    x = np.random.randn(_dim_0, 3, )
    R = np.random.randn(_dim_0, 3, 3)
    return x, R


if __name__ == '__main__':
    x, R = generateRandomData()
    print("x:", x)
    print("R:", R)
    func_value = demo38(x, R)
    print("func_value: ", func_value)