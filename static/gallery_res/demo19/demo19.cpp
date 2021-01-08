/*
from linearalgebra: tr

`J₃` = [1_3,3]
`k_angle(Dₘ)` = 3(sqrt(2)v)^(2/3)(7/4||Dₘ||_F^2-1/4tr(`J₃` (Dₘ)ᵀDₘ))⁻¹

where

Dₘ: ℝ^(3×3)  
v: ℝ
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo19
 *
 * @param Dₘ  ℝ^(3×3)
 * @param v  ℝ
 * @return k_angle_left_parenthesis_Dₘ_right_parenthesis
 */
double demo19(
    const Eigen::Matrix<double, 3, 3> & Dₘ,
    const double & v)
{
    Eigen::Matrix<double, 3, 3> _J₃_0;
    _J₃_0 << Eigen::MatrixXd::Ones(3, 3);
    Eigen::Matrix<double, 3, 3> J₃ = _J₃_0;

    double k_angle_left_parenthesis_Dₘ_right_parenthesis = 3 * pow((sqrt(2) * v), (2 / 3)) * 1 / ((7 / 4 * pow((Dₘ).norm(), 2) - 1 / 4 * (J₃ * (Dₘ).transpose() * Dₘ).trace()));

    return k_angle_left_parenthesis_Dₘ_right_parenthesis;
}


void generateRandomData(Eigen::Matrix<double, 3, 3> & Dₘ,
    double & v)
{
    v = rand() % 10;
    Dₘ = Eigen::MatrixXd::Random(3, 3);
}


int main(int argc, char *argv[])
{
    Eigen::Matrix<double, 3, 3> Dₘ;
    double v;
    generateRandomData(Dₘ, v);
    double func_value = demo19(Dₘ, v);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}