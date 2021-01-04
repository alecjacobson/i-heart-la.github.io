"""
sum_i f_i²p_i - (sum_i f_i p_i)²

where

f_i: ℝ  
p_i: ℝ  
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo6(f, p):
    """
    :param :f : ℝ
    :param :p : ℝ
    """
    f = np.asarray(f)
    p = np.asarray(p)

    _dim_0 = f.shape[0]
    assert f.shape == (_dim_0,)
    assert p.shape == (_dim_0,)

    _sum_0 = 0
    for i in range(1, len(f)+1):
        _sum_0 += np.power(f[i-1], 2) * p[i-1]
    _sum_1 = 0
    for i in range(1, len(f)+1):
        _sum_1 += f[i-1] * p[i-1]
    ret = _sum_0 - np.power((_sum_1), 2)
    return ret


def generateRandomData():
    _dim_0 = np.random.randint(10)
    f = np.random.randn(_dim_0)
    p = np.random.randn(_dim_0)
    return f, p


if __name__ == '__main__':
    f, p = generateRandomData()
    print("f:", f)
    print("p:", p)
    func_value = demo6(f, p)
    print("func_value: ", func_value)