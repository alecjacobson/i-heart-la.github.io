/*
from trigonometry: sin, cos
`x(θ, ϕ)` = [Rcos(θ)cos(ϕ)
             Rsin(θ)cos(ϕ)
             Rsin(ϕ)]

where

ϕ ∈ ℝ : angle between 0 and 2π
θ ∈ ℝ : angle between -π/2 and π/2
R ∈ ℝ : the radius of the sphere
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct pmp_32ResultType {
    Eigen::Matrix<double, 3, 1> x_left_parenthesis_θ_comma__ϕ_right_parenthesis;
    pmp_32ResultType(const Eigen::Matrix<double, 3, 1> & x_left_parenthesis_θ_comma__ϕ_right_parenthesis)
    : x_left_parenthesis_θ_comma__ϕ_right_parenthesis(x_left_parenthesis_θ_comma__ϕ_right_parenthesis)
    {}
};

/**
 * pmp_32
 *
 * @param ϕ  angle between 0 and 2π
 * @param θ  angle between -π/2 and π/2
 * @param R  the radius of the sphere
 * @return x_left_parenthesis_θ_comma__ϕ_right_parenthesis
 */
pmp_32ResultType pmp_32(
    const double & ϕ,
    const double & θ,
    const double & R)
{
    Eigen::Matrix<double, 3, 1> x_left_parenthesis_θ_comma__ϕ_right_parenthesis__0;
    x_left_parenthesis_θ_comma__ϕ_right_parenthesis__0 << R * cos(θ) * cos(ϕ),
    R * sin(θ) * cos(ϕ),
    R * sin(ϕ);
    Eigen::Matrix<double, 3, 1> x_left_parenthesis_θ_comma__ϕ_right_parenthesis = x_left_parenthesis_θ_comma__ϕ_right_parenthesis__0;

    return pmp_32ResultType(x_left_parenthesis_θ_comma__ϕ_right_parenthesis);
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
    srand((int)time(NULL));
    double ϕ;
    double θ;
    double R;
    generateRandomData(ϕ, θ, R);
    pmp_32ResultType func_value = pmp_32(ϕ, θ, R);
    std::cout<<"return value:\n"<<func_value.x_left_parenthesis_θ_comma__ϕ_right_parenthesis<<std::endl;
    return 0;
}