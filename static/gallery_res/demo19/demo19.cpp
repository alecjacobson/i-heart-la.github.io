/*
from linearalgebra: tr

`k_angle(Dₘ)` = 3(sqrt(2)v)^(2/3)(7/4||`Dₘ`||_F^2-1/4tr(J_3 `Dₘ`ᵀ`Dₘ`))⁻¹

where

`Dₘ`: ℝ^(n×n) 
J_i: ℝ^(n×n) 
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
 * @param Dₘ  ℝ^(n×n)
 * @param J  ℝ^(n×n)
 * @param v  ℝ
 * @return k_angle_left_parenthesis_Dₘ_right_parenthesis
 */
double demo19(
    const Eigen::MatrixXd & Dₘ,
    const std::vector<Eigen::MatrixXd> & J,
    const double & v)
{
    const long n = J[0].cols();
    const long _dim_0 = J.size();
    assert( Dₘ.rows() == n );
    assert( Dₘ.cols() == n );
    assert( J.size() == _dim_0 );
    for( const auto& el : J ) {
        assert( el.rows() == n );
        assert( el.cols() == n );
    }

    double k_angle_left_parenthesis_Dₘ_right_parenthesis = 3 * pow((sqrt(2) * v), (2 / 3)) * 1 / ((7 / 4 * pow((Dₘ).norm(), 2) - 1 / 4 * (J.at(3-1) * Dₘ.transpose() * Dₘ).trace()));

    return k_angle_left_parenthesis_Dₘ_right_parenthesis;
}


void generateRandomData(Eigen::MatrixXd & Dₘ,
    std::vector<Eigen::MatrixXd> & J,
    double & v)
{
    v = rand() % 10;
    const int n = rand()%10;
    const int _dim_0 = rand()%10;
    Dₘ = Eigen::MatrixXd::Random(n, n);
    J.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        J[i] = Eigen::MatrixXd::Random(n, n);
    }
}


int main(int argc, char *argv[])
{
    Eigen::MatrixXd Dₘ;
    std::vector<Eigen::MatrixXd> J;
    double v;
    generateRandomData(Dₘ, J, v);
    double func_value = demo19(Dₘ, J, v);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}