/*
from linearalgebra: tr

`J₃` = 1₃,₃
`κ_angle(Dₘ)` = 3(√2 v)^(2/3)(7/4‖`Dₘ`‖_F^2-1/4tr(`J₃``Dₘ`ᵀ`Dₘ`))⁻¹

where

`Dₘ` ∈ ℝ^(3×3)  
v ∈ ℝ
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct anisotropic_elasticity_47ResultType {
    Eigen::Matrix<double, 3, 3> J₃;
    double κ_angle_left_parenthesis_Dₘ_right_parenthesis;
    anisotropic_elasticity_47ResultType(const Eigen::Matrix<double, 3, 3> & J₃,
               const double & κ_angle_left_parenthesis_Dₘ_right_parenthesis)
    : J₃(J₃),
    κ_angle_left_parenthesis_Dₘ_right_parenthesis(κ_angle_left_parenthesis_Dₘ_right_parenthesis)
    {}
};

anisotropic_elasticity_47ResultType anisotropic_elasticity_47(
    const Eigen::Matrix<double, 3, 3> & Dₘ,
    const double & v)
{
    Eigen::Matrix<double, 3, 3> J₃ = Eigen::MatrixXd::Ones(3, 3);

    double κ_angle_left_parenthesis_Dₘ_right_parenthesis = 3 * pow((sqrt(2) * v), (2 / double(3))) * 1 / ((7 / double(4) * pow((Dₘ).norm(), 2) - 1 / double(4) * (J₃ * Dₘ.transpose() * Dₘ).trace()));

    return anisotropic_elasticity_47ResultType(J₃, κ_angle_left_parenthesis_Dₘ_right_parenthesis);
}


void generateRandomData(Eigen::Matrix<double, 3, 3> & Dₘ,
    double & v)
{
    v = rand() % 10;
    Dₘ = Eigen::MatrixXd::Random(3, 3);
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    Eigen::Matrix<double, 3, 3> Dₘ;
    double v;
    generateRandomData(Dₘ, v);
    anisotropic_elasticity_47ResultType func_value = anisotropic_elasticity_47(Dₘ, v);
    std::cout<<"return value:\n"<<func_value.κ_angle_left_parenthesis_Dₘ_right_parenthesis<<std::endl;
    return 0;
}