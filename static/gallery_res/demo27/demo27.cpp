/*
Ω = [`e₁` `e₂`][`k₁`   0
                 0    `k₂`] [`e₁`ᵀ
			                 `e₂`ᵀ]

where

`k₁`: ℝ  : control the desired kernel variance in either edge or orthogonal direction
`k₂`: ℝ  : control the desired kernel variance in either edge or orthogonal direction
`e₁`: ℝ ^ 3: orthogonal direction vectors
`e₂`: ℝ ^ 3: orthogonal direction vectors
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo27
 *
 * @param k₁  ℝ  : control the desired kernel variance in either edge or orthogonal direction
 * @param k₂  ℝ  : control the desired kernel variance in either edge or orthogonal direction
 * @param e₁  ℝ ^ 3: orthogonal direction vectors
 * @param e₂  ℝ ^ 3: orthogonal direction vectors
 * @return Ω
 */
Eigen::Matrix<double, 3, 3> demo27(
    const double & k₁,
    const double & k₂,
    const Eigen::Matrix<double, 3, 1> & e₁,
    const Eigen::Matrix<double, 3, 1> & e₂)
{
    Eigen::Matrix<double, 3, 2> _Ω_0;
    _Ω_0 << e₁, e₂;
    Eigen::Matrix<double, 2, 2> _Ω_1;
    _Ω_1 << k₁, 0,
    0, k₂;
    Eigen::Matrix<double, 2, 3> _Ω_2;
    _Ω_2 << e₁.transpose(),
    e₂.transpose();
    Eigen::Matrix<double, 3, 3> Ω = _Ω_0 * _Ω_1 * _Ω_2;

    return Ω;
}


void generateRandomData(double & k₁,
    double & k₂,
    Eigen::Matrix<double, 3, 1> & e₁,
    Eigen::Matrix<double, 3, 1> & e₂)
{
    k₁ = rand() % 10;
    k₂ = rand() % 10;
    e₁ = Eigen::VectorXd::Random(3);
    e₂ = Eigen::VectorXd::Random(3);
}


int main(int argc, char *argv[])
{
    double k₁;
    double k₂;
    Eigen::Matrix<double, 3, 1> e₁;
    Eigen::Matrix<double, 3, 1> e₂;
    generateRandomData(k₁, k₂, e₁, e₂);
    Eigen::Matrix<double, 3, 3> func_value = demo27(k₁, k₂, e₁, e₂);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}