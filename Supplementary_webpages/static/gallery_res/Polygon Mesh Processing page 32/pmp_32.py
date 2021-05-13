"""
from trigonometry: sin, cos
`x(θ, ϕ)` = [Rcos(θ)cos(ϕ)
             Rsin(θ)cos(ϕ)
             Rsin(ϕ)]

where

ϕ ∈ ℝ : angle between 0 and 2π
θ ∈ ℝ : angle between -π/2 and π/2
R ∈ ℝ : the radius of the sphere
"""
import numpy as np
import scipy
import scipy.linalg
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize


class pmp_32ResultType:
    def __init__( self, x_left_parenthesis_θ_comma__ϕ_right_parenthesis):
        self.x_left_parenthesis_θ_comma__ϕ_right_parenthesis = x_left_parenthesis_θ_comma__ϕ_right_parenthesis


def pmp_32(ϕ, θ, R):
    """
    :param :ϕ : angle between 0 and 2π
    :param :θ : angle between -π/2 and π/2
    :param :R : the radius of the sphere
    """
    assert np.ndim(ϕ) == 0
    assert np.ndim(θ) == 0
    assert np.ndim(R) == 0

    x_left_parenthesis_θ_comma__ϕ_right_parenthesis__0 = np.zeros((3, 1))
    x_left_parenthesis_θ_comma__ϕ_right_parenthesis__0[0] = [R * np.cos(θ) * np.cos(ϕ)]
    x_left_parenthesis_θ_comma__ϕ_right_parenthesis__0[1] = [R * np.sin(θ) * np.cos(ϕ)]
    x_left_parenthesis_θ_comma__ϕ_right_parenthesis__0[2] = [R * np.sin(ϕ)]
    x_left_parenthesis_θ_comma__ϕ_right_parenthesis = x_left_parenthesis_θ_comma__ϕ_right_parenthesis__0
    return pmp_32ResultType(x_left_parenthesis_θ_comma__ϕ_right_parenthesis)


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
    func_value = pmp_32(ϕ, θ, R)
    print("return value: ", func_value.x_left_parenthesis_θ_comma__ϕ_right_parenthesis)