/*
[P_1  0
 0  P_3][  L              0
    P_3^TCP_2^TU^(-1)   -Lbar][U  L^(-1)P_1^TB
                               0     `U_bar`][P_2   0
                                               0   I_4]
where

P_i: ℝ^(4×4) 
B: ℝ^(4×4) 
C: ℝ^(4×4) 
L: ℝ^(4×4) 
Lbar: ℝ^(4×4) 
U: ℝ^(4×4) 
`U_bar`: ℝ^(4×4)
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo15
 *
 * @param P  ℝ^(4×4)
 * @param B  ℝ^(4×4)
 * @param C  ℝ^(4×4)
 * @param L  ℝ^(4×4)
 * @param Lbar  ℝ^(4×4)
 * @param U  ℝ^(4×4)
 * @param U_bar  ℝ^(4×4)
 * @return ret
 */
Eigen::Matrix<double, 8, 8> demo15(
    const std::vector<Eigen::Matrix<double, 4, 4>> & P,
    const Eigen::Matrix<double, 4, 4> & B,
    const Eigen::Matrix<double, 4, 4> & C,
    const Eigen::Matrix<double, 4, 4> & L,
    const Eigen::Matrix<double, 4, 4> & Lbar,
    const Eigen::Matrix<double, 4, 4> & U,
    const Eigen::Matrix<double, 4, 4> & U_bar)
{
    const long _dim_0 = P.size();
    assert( P.size() == _dim_0 );

    Eigen::Matrix<double, 8, 8> _ret_0;
    _ret_0 << P.at(1-1), Eigen::MatrixXd::Zero(4, 4),
    Eigen::MatrixXd::Zero(4, 4), P.at(3-1);
    Eigen::Matrix<double, 8, 8> _ret_1;
    _ret_1 << L, Eigen::MatrixXd::Zero(4, 4),
    P.at(3-1).transpose() * C * P.at(2-1).transpose() * U.inverse(), -Lbar;
    Eigen::Matrix<double, 8, 8> _ret_2;
    _ret_2 << U, L.inverse() * P.at(1-1).transpose() * B,
    Eigen::MatrixXd::Zero(4, 4), U_bar;
    Eigen::Matrix<double, 8, 8> _ret_3;
    _ret_3 << P.at(2-1), Eigen::MatrixXd::Zero(4, 4),
    Eigen::MatrixXd::Zero(4, 4), Eigen::MatrixXd::Identity(4, 4);
    Eigen::Matrix<double, 8, 8> ret = _ret_0 * _ret_1 * _ret_2 * _ret_3;
    return ret;
}


void generateRandomData(std::vector<Eigen::Matrix<double, 4, 4>> & P,
    Eigen::Matrix<double, 4, 4> & B,
    Eigen::Matrix<double, 4, 4> & C,
    Eigen::Matrix<double, 4, 4> & L,
    Eigen::Matrix<double, 4, 4> & Lbar,
    Eigen::Matrix<double, 4, 4> & U,
    Eigen::Matrix<double, 4, 4> & U_bar)
{
    const int _dim_0 = rand()%10;
    P.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        P[i] = Eigen::MatrixXd::Random(4, 4);
    }
    B = Eigen::MatrixXd::Random(4, 4);
    C = Eigen::MatrixXd::Random(4, 4);
    L = Eigen::MatrixXd::Random(4, 4);
    Lbar = Eigen::MatrixXd::Random(4, 4);
    U = Eigen::MatrixXd::Random(4, 4);
    U_bar = Eigen::MatrixXd::Random(4, 4);
}


int main(int argc, char *argv[])
{
    std::vector<Eigen::Matrix<double, 4, 4>> P;
    Eigen::Matrix<double, 4, 4> B;
    Eigen::Matrix<double, 4, 4> C;
    Eigen::Matrix<double, 4, 4> L;
    Eigen::Matrix<double, 4, 4> Lbar;
    Eigen::Matrix<double, 4, 4> U;
    Eigen::Matrix<double, 4, 4> U_bar;
    generateRandomData(P, B, C, L, Lbar, U, U_bar);
    Eigen::Matrix<double, 8, 8> func_value = demo15(P, B, C, L, Lbar, U, U_bar);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}