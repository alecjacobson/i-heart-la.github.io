/*
from trigonometry: sin, cos
`x(θ, ϕ)` = [Rcos(θ)cos(ϕ)
             Rsin(θ)cos(ϕ)
             Rsin(ϕ)]

where

ϕ: ℝ : angle between 0 and 2π
θ: ℝ : angle between -π/2 and π/2
R: ℝ : the radius of the sphere
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo1
 *
 * @param ϕ  ℝ : angle between 0 and 2π
 * @param θ  ℝ : angle between -π/2 and π/2
 * @param R  ℝ : the radius of the sphere
 * @return x_left_parenthesis_θ_comma__ϕ_right_parenthesis
 */
Eigen::Matrix<double, 3, 1> demo1(
    const double & ϕ,
    const double & θ,
    const double & R)
{
    Eigen::Matrix<double, 3, 1> _x_left_parenthesis_θ_comma__ϕ_right_parenthesis__0;
    _x_left_parenthesis_θ_comma__ϕ_right_parenthesis__0 << R * cos(θ) * cos(ϕ),
    R * sin(θ) * cos(ϕ),
    R * sin(ϕ);
    Eigen::Matrix<double, 3, 1> x_left_parenthesis_θ_comma__ϕ_right_parenthesis = _x_left_parenthesis_θ_comma__ϕ_right_parenthesis__0;

    return x_left_parenthesis_θ_comma__ϕ_right_parenthesis;
}


void generateRandomData(double & ϕ,
    double & θ,
    double & R)
{
    ϕ = rand() % 10;
    θ = rand() % 10;
    R = rand() % 10;
}


int main(int argc, char *argv[])
{
    double ϕ;
    double θ;
    double R;
    generateRandomData(ϕ, θ, R);
    Eigen::Matrix<double, 3, 1> func_value = demo1(ϕ, θ, R);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}