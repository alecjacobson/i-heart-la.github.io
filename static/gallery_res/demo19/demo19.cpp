/*
from linearalgebra: tr

`k_angle(D_m)` = 3(sqrt(2)v)^(2/3)(7/4||`D_m`||_F^2 -1/4tr(J_3 `D_m`^T`D_m`))^(-1)

where

`D_m`: ℝ^(n×n) 
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
 * @param D_m  ℝ^(n×n)
 * @param J  ℝ^(n×n)
 * @param v  ℝ
 * @return k_angle_left_parenthesis_D_m_right_parenthesis
 */
double demo19(
    const Eigen::MatrixXd & D_m,
    const std::vector<Eigen::MatrixXd> & J,
    const double & v)
{
    const long n = J[0].cols();
    const long _dim_0 = J.size();
    assert( D_m.rows() == n );
    assert( D_m.cols() == n );
    assert( J.size() == _dim_0 );
    for( const auto& el : J ) {
        assert( el.rows() == n );
        assert( el.cols() == n );
    }

    double k_angle_left_parenthesis_D_m_right_parenthesis = 3 * pow((sqrt(2) * v), (2 / 3)) * 1 / ((7 / 4 * pow((D_m).norm(), 2) - 1 / 4 * (J.at(3-1) * D_m.transpose() * D_m).trace()));

    return k_angle_left_parenthesis_D_m_right_parenthesis;
}


void generateRandomData(Eigen::MatrixXd & D_m,
    std::vector<Eigen::MatrixXd> & J,
    double & v)
{
    v = rand() % 10;
    const int n = rand()%10;
    const int _dim_0 = rand()%10;
    D_m = Eigen::MatrixXd::Random(n, n);
    J.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        J[i] = Eigen::MatrixXd::Random(n, n);
    }
}


int main(int argc, char *argv[])
{
    Eigen::MatrixXd D_m;
    std::vector<Eigen::MatrixXd> J;
    double v;
    generateRandomData(D_m, J, v);
    double func_value = demo19(D_m, J, v);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}