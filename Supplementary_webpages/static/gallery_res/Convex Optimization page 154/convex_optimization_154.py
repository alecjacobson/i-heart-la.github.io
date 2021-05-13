"""
∑_i f_i²p_i - (∑_i f_i p_i)²

where

f_i ∈ ℝ  
p_i ∈ ℝ  
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class convex_optimization_154ResultType:
    def __init__( self, ret):
        self.ret = ret


def convex_optimization_154(f, p):
    f = np.asarray(f, dtype=np.float64)
    p = np.asarray(p, dtype=np.float64)

    dim_0 = f.shape[0]
    assert f.shape == (dim_0,)
    assert p.shape == (dim_0,)

    sum_0 = 0
    for i in range(1, len(p)+1):
        sum_0 += np.power(f[i-1], 2) * p[i-1]
    sum_1 = 0
    for i in range(1, len(p)+1):
        sum_1 += f[i-1] * p[i-1]
    ret = sum_0 - np.power((sum_1), 2)
    return convex_optimization_154ResultType(ret)


def generateRandomData():
    dim_0 = np.random.randint(10)
    f = np.random.randn(dim_0)
    p = np.random.randn(dim_0)
    return f, p


if __name__ == '__main__':
    f, p = generateRandomData()
    print("f:", f)
    print("p:", p)
    func_value = convex_optimization_154(f, p)
    print("return value: ", func_value.ret)