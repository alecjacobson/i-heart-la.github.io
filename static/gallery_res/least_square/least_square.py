"""
given
p_i: ℝ^3: points on lines
d_i: ℝ^3: unit directions along lines

k_i = (p_i - (p_i⋅d_i)d_i)
a_i = [1,0,0]^T - d_i,0 d_i
b_i = [0,1,0]^T - d_i,1 d_i
c_i = [0,0,1]^T - d_i,2 d_i

 
M = [ (∑_i( a_i,0 - d_i,0 (d_i⋅a_i) ))    (∑_i( a_i,1 - d_i,1 (d_i⋅a_i) ))    (∑_i( a_i,2 - d_i,2 (d_i⋅a_i) ))
      (∑_i( b_i,0 - d_i,0 (d_i⋅b_i) ))    (∑_i( b_i,1 - d_i,1 (d_i⋅b_i) ))    (∑_i( b_i,2 - d_i,2 (d_i⋅b_i) ))
      (∑_i( c_i,0 - d_i,0 (d_i⋅c_i) ))    (∑_i( c_i,1 - d_i,1 (d_i⋅c_i) ))    (∑_i( c_i,2 - d_i,2 (d_i⋅c_i) )) ]

r = [ ∑_i( k_i⋅a_i )
      ∑_i( k_i⋅b_i )
      ∑_i( k_i⋅c_i ) ]

q = M^(-1) r


"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def myExpression(p, d):
    """
    :param :p : ℝ^3: points on lines
    :param :d : ℝ^3: unit directions along lines
    """
    p = np.asarray(p, dtype=np.float64)
    d = np.asarray(d, dtype=np.float64)

    _dim_0 = p.shape[0]
    assert p.shape == (_dim_0, 3, 1)
    assert d.shape == (_dim_0, 3, 1)

    k = np.zeros((_dim_0, 3, 1))
    for i in range(1, _dim_0+1):
        k[i] = (p[i-1] - (np.dot((p[i-1]).ravel(), (d[i-1]).ravel())) * d[i-1])

    _a_i_0 = np.zeros((1, 3))
    _a_i_0[0] = [1, 0, 0]
    a = np.zeros((_dim_0, 3, 1))
    for i in range(1, _dim_0+1):
        a[i] = _a_i_0.T - d[i-1][0-1] * d[i-1]

    _b_i_0 = np.zeros((1, 3))
    _b_i_0[0] = [0, 1, 0]
    b = np.zeros((_dim_0, 3, 1))
    for i in range(1, _dim_0+1):
        b[i] = _b_i_0.T - d[i-1][1-1] * d[i-1]

    _c_i_0 = np.zeros((1, 3))
    _c_i_0[0] = [0, 0, 1]
    c = np.zeros((_dim_0, 3, 1))
    for i in range(1, _dim_0+1):
        c[i] = _c_i_0.T - d[i-1][2-1] * d[i-1]

    _sum_0 = 0
    for i in range(1, len(a)+1):
        _sum_0 += (a[i-1][0-1] - d[i-1][0-1] * (np.dot((d[i-1]).ravel(), (a[i-1]).ravel())))
    _sum_1 = 0
    for i in range(1, len(d)+1):
        _sum_1 += (a[i-1][1-1] - d[i-1][1-1] * (np.dot((d[i-1]).ravel(), (a[i-1]).ravel())))
    _sum_2 = 0
    for i in range(1, len(a)+1):
        _sum_2 += (a[i-1][2-1] - d[i-1][2-1] * (np.dot((d[i-1]).ravel(), (a[i-1]).ravel())))
    _sum_3 = 0
    for i in range(1, len(d)+1):
        _sum_3 += (b[i-1][0-1] - d[i-1][0-1] * (np.dot((d[i-1]).ravel(), (b[i-1]).ravel())))
    _sum_4 = 0
    for i in range(1, len(d)+1):
        _sum_4 += (b[i-1][1-1] - d[i-1][1-1] * (np.dot((d[i-1]).ravel(), (b[i-1]).ravel())))
    _sum_5 = 0
    for i in range(1, len(d)+1):
        _sum_5 += (b[i-1][2-1] - d[i-1][2-1] * (np.dot((d[i-1]).ravel(), (b[i-1]).ravel())))
    _sum_6 = 0
    for i in range(1, len(d)+1):
        _sum_6 += (c[i-1][0-1] - d[i-1][0-1] * (np.dot((d[i-1]).ravel(), (c[i-1]).ravel())))
    _sum_7 = 0
    for i in range(1, len(d)+1):
        _sum_7 += (c[i-1][1-1] - d[i-1][1-1] * (np.dot((d[i-1]).ravel(), (c[i-1]).ravel())))
    _sum_8 = 0
    for i in range(1, len(d)+1):
        _sum_8 += (c[i-1][2-1] - d[i-1][2-1] * (np.dot((d[i-1]).ravel(), (c[i-1]).ravel())))
    _M_0 = np.zeros((3, 3))
    _M_0[0] = [(_sum_0), (_sum_1), (_sum_2)]
    _M_0[1] = [(_sum_3), (_sum_4), (_sum_5)]
    _M_0[2] = [(_sum_6), (_sum_7), (_sum_8)]
    M = _M_0

    _sum_9 = 0
    for i in range(1, len(a)+1):
        _sum_9 += (np.dot((k[i-1]).ravel(), (a[i-1]).ravel()))
    _sum_10 = 0
    for i in range(1, len(b)+1):
        _sum_10 += (np.dot((k[i-1]).ravel(), (b[i-1]).ravel()))
    _sum_11 = 0
    for i in range(1, len(c)+1):
        _sum_11 += (np.dot((k[i-1]).ravel(), (c[i-1]).ravel()))
    _r_0 = np.zeros((3, 1))
    _r_0[0] = [_sum_9]
    _r_0[1] = [_sum_10]
    _r_0[2] = [_sum_11]
    r = _r_0

    q = np.linalg.inv(M) @ r

    return q


def generateRandomData():
    _dim_0 = np.random.randint(10)
    p = np.random.randn(_dim_0, 3, 1)
    d = np.random.randn(_dim_0, 3, 1)
    return p, d


if __name__ == '__main__':
    p, d = generateRandomData()
    print("p:", p)
    print("d:", d)
    func_value = myExpression(p, d)
    print("func_value: ", func_value)