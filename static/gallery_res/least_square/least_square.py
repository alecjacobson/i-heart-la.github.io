"""
given
p_i: ℝ^3: points on lines
d_i: ℝ^3: unit directions along lines

P_i = ( I_3 - d_i d_iᵀ )
q = ( ∑_i P_iᵀP_i )⁻¹ ( ∑_i P_iᵀP_i p_i )
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def least_square(p, d):
    """
    :param :p : ℝ^3: points on lines
    :param :d : ℝ^3: unit directions along lines
    """
    p = np.asarray(p, dtype=np.float64)
    d = np.asarray(d, dtype=np.float64)

    _dim_0 = p.shape[0]
    assert p.shape == (_dim_0, 3, )
    assert d.shape == (_dim_0, 3, )

    P = np.zeros((_dim_0, 3, 3))
    for i in range(1, _dim_0+1):
        P[i-1] = (np.identity(3) - (d[i-1]).reshape(3, 1) @ d[i-1].T.reshape(1, 3))

    _sum_0 = np.zeros((3, 3))
    for i in range(1, len(P)+1):
        _sum_0 += P[i-1].T @ P[i-1]
    _sum_1 = np.zeros((3, ))
    for i in range(1, len(P)+1):
        _sum_1 += P[i-1].T @ P[i-1] @ p[i-1]
    q = np.linalg.inv((_sum_0)) @ (_sum_1)

    return q


def generateRandomData():
    _dim_0 = np.random.randint(10)
    p = np.random.randn(_dim_0, 3, )
    d = np.random.randn(_dim_0, 3, )
    return p, d


if __name__ == '__main__':
    p, d = generateRandomData()
    print("p:", p)
    print("d:", d)
    func_value = least_square(p, d)
    print("func_value: ", func_value)