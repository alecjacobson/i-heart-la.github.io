"""
from trigonometry: sin, cos
`x(θ, ϕ)` = [Rcos(θ)cos(ϕ)
             Rsin(θ)cos(ϕ)
             Rsin(ϕ)]

where

ϕ: ℝ : angle between 0 and 2π
θ: ℝ : angle between -π/2 and π/2
R: ℝ : the radius of the sphere
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


def demo1(ϕ, θ, R):
    """
    :param :ϕ : ℝ : angle between 0 and 2π
    :param :θ : ℝ : angle between -π/2 and π/2
    :param :R : ℝ : the radius of the sphere
    """
    assert np.ndim(ϕ) == 0
    assert np.ndim(θ) == 0
    assert np.ndim(R) == 0

    _x_left_parenthesis_θ_comma__ϕ_right_parenthesis__0 = np.zeros((3, 1))
    _x_left_parenthesis_θ_comma__ϕ_right_parenthesis__0[0] = [R * np.cos(θ) * np.cos(ϕ)]
    _x_left_parenthesis_θ_comma__ϕ_right_parenthesis__0[1] = [R * np.sin(θ) * np.cos(ϕ)]
    _x_left_parenthesis_θ_comma__ϕ_right_parenthesis__0[2] = [R * np.sin(ϕ)]
    x_left_parenthesis_θ_comma__ϕ_right_parenthesis = _x_left_parenthesis_θ_comma__ϕ_right_parenthesis__0

    return x_left_parenthesis_θ_comma__ϕ_right_parenthesis


def generateRandomData():
    ϕ = np.random.randn()
    θ = np.random.randn()
    R = np.random.randn()
    return ϕ, θ, R


if __name__ == '__main__':
    ϕ, θ, R = generateRandomData()
    print("ϕ:", ϕ)
    print("θ:", θ)
    print("R:", R)
    func_value = demo1(ϕ, θ, R)
    print("func_value: ", func_value)