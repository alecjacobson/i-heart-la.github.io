"""

E = 1/`σ_N`^2`E_I` + sum_(j for j>1) α_j²/`σ_S`_j^2  + sum_(j for j>1) β_j²/`σ_T`_j^2   + sum_j (ρ_j-`ρ_bar`_j)²/`σ_ρ`_j^2 

where

`σ_N`: ℝ 
`E_I`: ℝ
α_i : ℝ
β_i : ℝ
`σ_S`_i: ℝ 
`σ_T`_i: ℝ 
ρ_i: ℝ 
`ρ_bar`_i: ℝ 
`σ_ρ`_i: ℝ 
ā_i: ℝ 
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo24(σ_N, E_I, α, β, σ_S, σ_T, ρ, ρ_bar, σ_ρ, ā):
    """
    :param :σ_N : ℝ
    :param :E_I : ℝ
    :param :α : ℝ
    :param :β : ℝ
    :param :σ_S : ℝ
    :param :σ_T : ℝ
    :param :ρ : ℝ
    :param :ρ_bar : ℝ
    :param :σ_ρ : ℝ
    :param :ā : ℝ
    """
    α = np.asarray(α)
    β = np.asarray(β)
    σ_S = np.asarray(σ_S)
    σ_T = np.asarray(σ_T)
    ρ = np.asarray(ρ)
    ρ_bar = np.asarray(ρ_bar)
    σ_ρ = np.asarray(σ_ρ)
    ā = np.asarray(ā)

    _dim_0 = α.shape[0]
    assert np.ndim(σ_N) == 0
    assert np.ndim(E_I) == 0
    assert α.shape == (_dim_0,)
    assert β.shape == (_dim_0,)
    assert σ_S.shape == (_dim_0,)
    assert σ_T.shape == (_dim_0,)
    assert ρ.shape == (_dim_0,)
    assert ρ_bar.shape == (_dim_0,)
    assert σ_ρ.shape == (_dim_0,)
    assert ā.shape == (_dim_0,)

    _sum_0 = 0
    for j in range(1, len(α)+1):
        if(j > 1):
            _sum_0 += np.power(α[j-1], 2) / np.power(σ_S[j-1], 2)
    _sum_1 = 0
    for j in range(1, len(β)+1):
        if(j > 1):
            _sum_1 += np.power(β[j-1], 2) / np.power(σ_T[j-1], 2)
    _sum_2 = 0
    for j in range(1, len(ρ)+1):
        _sum_2 += np.power((ρ[j-1] - ρ_bar[j-1]), 2) / np.power(σ_ρ[j-1], 2)
    E = 1 / np.power(σ_N, 2) * E_I + _sum_0 + _sum_1 + _sum_2

    return E


def generateRandomData():
    σ_N = np.random.randn()
    E_I = np.random.randn()
    _dim_0 = np.random.randint(10)
    α = np.random.randn(_dim_0)
    β = np.random.randn(_dim_0)
    σ_S = np.random.randn(_dim_0)
    σ_T = np.random.randn(_dim_0)
    ρ = np.random.randn(_dim_0)
    ρ_bar = np.random.randn(_dim_0)
    σ_ρ = np.random.randn(_dim_0)
    ā = np.random.randn(_dim_0)
    return σ_N, E_I, α, β, σ_S, σ_T, ρ, ρ_bar, σ_ρ, ā


if __name__ == '__main__':
    σ_N, E_I, α, β, σ_S, σ_T, ρ, ρ_bar, σ_ρ, ā = generateRandomData()
    print("σ_N:", σ_N)
    print("E_I:", E_I)
    print("α:", α)
    print("β:", β)
    print("σ_S:", σ_S)
    print("σ_T:", σ_T)
    print("ρ:", ρ)
    print("ρ_bar:", ρ_bar)
    print("σ_ρ:", σ_ρ)
    print("ā:", ā)
    func_value = demo24(σ_N, E_I, α, β, σ_S, σ_T, ρ, ρ_bar, σ_ρ, ā)
    print("func_value: ", func_value)