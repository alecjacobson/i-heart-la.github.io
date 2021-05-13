"""
from trigonometry: cos

∑_i cos²(θ)((p_i - q_i)⋅n_i +((p_i+q_i)×n_i)⋅ã+n_i⋅t̃)² 

where

θ ∈ ℝ: angle of rotation
p_i ∈ ℝ^3
q_i ∈ ℝ^3
n_i ∈ ℝ^3
ã ∈ ℝ^3
t̃ ∈ ℝ^3
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class symmetric_objective_function_9ResultType:
    def __init__( self, ret):
        self.ret = ret


def symmetric_objective_function_9(θ, p, q, n, ã, t̃):
    """
    :param :θ : angle of rotation
    """
    p = np.asarray(p, dtype=np.float64)
    q = np.asarray(q, dtype=np.float64)
    n = np.asarray(n, dtype=np.float64)
    ã = np.asarray(ã, dtype=np.float64)
    t̃ = np.asarray(t̃, dtype=np.float64)

    dim_0 = p.shape[0]
    assert np.ndim(θ) == 0
    assert p.shape == (dim_0, 3, )
    assert q.shape == (dim_0, 3, )
    assert n.shape == (dim_0, 3, )
    assert ã.shape == (3,)
    assert t̃.shape == (3,)

    sum_0 = 0
    for i in range(1, len(p)+1):
        sum_0 += np.power(np.cos(θ), 2) * np.power((np.dot(((p[i-1] - q[i-1])).ravel(), (n[i-1]).ravel()) + np.dot(((np.cross((p[i-1] + q[i-1]), n[i-1]))).ravel(), (ã).ravel()) + np.dot((n[i-1]).ravel(), (t̃).ravel())), 2)
    ret = sum_0
    return symmetric_objective_function_9ResultType(ret)


def generateRandomData():
    θ = np.random.randn()
    dim_0 = np.random.randint(10)
    p = np.random.randn(dim_0, 3, )
    q = np.random.randn(dim_0, 3, )
    n = np.random.randn(dim_0, 3, )
    ã = np.random.randn(3)
    t̃ = np.random.randn(3)
    return θ, p, q, n, ã, t̃


if __name__ == '__main__':
    θ, p, q, n, ã, t̃ = generateRandomData()
    print("θ:", θ)
    print("p:", p)
    print("q:", q)
    print("n:", n)
    print("ã:", ã)
    print("t̃:", t̃)
    func_value = symmetric_objective_function_9(θ, p, q, n, ã, t̃)
    print("return value: ", func_value.ret)