/*
`T_1` = 1/sqrt(2)U V^T
c = [0 0 0
                   0 0 -1
                   0 1 0]

where

U: ℝ^(3×3) 
V: ℝ^(3×3)
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo22
 *
 * @param U  ℝ^(3×3)
 * @param V  ℝ^(3×3)
 * @return c
 */
Eigen::Matrix<double, 3, 3> demo22(
    const Eigen::Matrix<double, 3, 3> & U,
    const Eigen::Matrix<double, 3, 3> & V)
{
    Eigen::Matrix<double, 3, 3> T_1 = 1 / sqrt(2) * U * V.transpose();

    Eigen::Matrix<double, 3, 3> _c_0;
    _c_0 << 0, 0, 0,
    0, 0, -1,
    0, 1, 0;
    Eigen::Matrix<double, 3, 3> c = _c_0;

    return c;
}


void generateRandomData(Eigen::Matrix<double, 3, 3> & U,
    Eigen::Matrix<double, 3, 3> & V)
{
    U = Eigen::MatrixXd::Random(3, 3);
    V = Eigen::MatrixXd::Random(3, 3);
}


int main(int argc, char *argv[])
{
    Eigen::Matrix<double, 3, 3> U;
    Eigen::Matrix<double, 3, 3> V;
    generateRandomData(U, V);
    Eigen::Matrix<double, 3, 3> func_value = demo22(U, V);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}