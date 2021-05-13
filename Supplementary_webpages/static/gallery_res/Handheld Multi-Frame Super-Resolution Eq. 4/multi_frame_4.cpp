/*
Ω = [`e₁` `e₂`][`k₁`   0
                 0    `k₂`] [`e₁`ᵀ
                             `e₂`ᵀ]

where

`k₁` ∈ ℝ  
`k₂` ∈ ℝ 
`e₁` ∈ ℝ ^ 3
`e₂` ∈ ℝ ^ 3
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct multi_frame_4ResultType {
    Eigen::Matrix<double, 3, 3> Ω;
    multi_frame_4ResultType(const Eigen::Matrix<double, 3, 3> & Ω)
    : Ω(Ω)
    {}
};

multi_frame_4ResultType multi_frame_4(
    const double & k₁,
    const double & k₂,
    const Eigen::Matrix<double, 3, 1> & e₁,
    const Eigen::Matrix<double, 3, 1> & e₂)
{
    Eigen::Matrix<double, 3, 2> Ω_0;
    Ω_0 << e₁, e₂;
    Eigen::Matrix<double, 2, 2> Ω_1;
    Ω_1 << k₁, 0,
    0, k₂;
    Eigen::Matrix<double, 2, 3> Ω_2;
    Ω_2 << e₁.transpose(),
    e₂.transpose();
    Eigen::Matrix<double, 3, 3> Ω = Ω_0 * Ω_1 * Ω_2;

    return multi_frame_4ResultType(Ω);
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
    srand((int)time(NULL));
    double k₁;
    double k₂;
    Eigen::Matrix<double, 3, 1> e₁;
    Eigen::Matrix<double, 3, 1> e₂;
    generateRandomData(k₁, k₂, e₁, e₂);
    multi_frame_4ResultType func_value = multi_frame_4(k₁, k₂, e₁, e₂);
    std::cout<<"return value:\n"<<func_value.Ω<<std::endl;
    return 0;
}