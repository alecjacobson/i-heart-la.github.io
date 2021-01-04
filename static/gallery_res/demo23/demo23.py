"""
`x(θ,v)` =  (`r bar`⋅`D_A`(`θ`, v)+`k_r``δ`(`θ`, v))/(`n bar`⋅`D_A`(`θ`, v)+`k_n``δ`(`θ`, v))
`y(θ,v)` =  (`s bar`⋅`D_A`(`θ`, v)+`k_s``δ`(`θ`, v))/(`n bar`⋅`D_A`(`θ`, v)+`k_n``δ`(`θ`, v))

where

`r bar`: ℝ^3
`s bar`: ℝ^3
`n bar`: ℝ^3
`θ`: ℝ 
v: ℝ 
`k_r`: ℝ 
`k_s`: ℝ 
`k_n`: ℝ 
`D_A`: ℝ,ℝ->ℝ^3
`δ`: ℝ,ℝ->ℝ  
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo23(r_bar, s_bar, n_bar, θ, v, k_r, k_s, k_n, D_A, δ):
    """
    :param :r_bar : ℝ^3
    :param :s_bar : ℝ^3
    :param :n_bar : ℝ^3
    :param :θ : ℝ
    :param :v : ℝ
    :param :k_r : ℝ
    :param :k_s : ℝ
    :param :k_n : ℝ
    :param :D_A : ℝ,ℝ->ℝ^3
    :param :δ : ℝ,ℝ->ℝ
    """
    r_bar = np.asarray(r_bar, dtype=np.float64)
    s_bar = np.asarray(s_bar, dtype=np.float64)
    n_bar = np.asarray(n_bar, dtype=np.float64)

    assert r_bar.shape == (3,)
    assert s_bar.shape == (3,)
    assert n_bar.shape == (3,)
    assert np.ndim(θ) == 0
    assert np.ndim(v) == 0
    assert np.ndim(k_r) == 0
    assert np.ndim(k_s) == 0
    assert np.ndim(k_n) == 0

    x_left_parenthesis_θ_comma_v_right_parenthesis = (np.dot((r_bar).ravel(), (D_A(θ, v)).ravel()) + k_r * δ(θ, v)) / (np.dot((n_bar).ravel(), (D_A(θ, v)).ravel()) + k_n * δ(θ, v))

    y_left_parenthesis_θ_comma_v_right_parenthesis = (np.dot((s_bar).ravel(), (D_A(θ, v)).ravel()) + k_s * δ(θ, v)) / (np.dot((n_bar).ravel(), (D_A(θ, v)).ravel()) + k_n * δ(θ, v))

    return y_left_parenthesis_θ_comma_v_right_parenthesis


def generateRandomData():
    θ = np.random.randn()
    v = np.random.randn()
    k_r = np.random.randn()
    k_s = np.random.randn()
    k_n = np.random.randn()
    r_bar = np.random.randn(3)
    s_bar = np.random.randn(3)
    n_bar = np.random.randn(3)
    def D_A(p0, p1):
        return np.random.randn(3)
    def δ(p0, p1):
        return np.random.randn()
    return r_bar, s_bar, n_bar, θ, v, k_r, k_s, k_n, D_A, δ


if __name__ == '__main__':
    r_bar, s_bar, n_bar, θ, v, k_r, k_s, k_n, D_A, δ = generateRandomData()
    print("r_bar:", r_bar)
    print("s_bar:", s_bar)
    print("n_bar:", n_bar)
    print("θ:", θ)
    print("v:", v)
    print("k_r:", k_r)
    print("k_s:", k_s)
    print("k_n:", k_n)
    print("D_A:", D_A)
    print("δ:", δ)
    func_value = demo23(r_bar, s_bar, n_bar, θ, v, k_r, k_s, k_n, D_A, δ)
    print("func_value: ", func_value)