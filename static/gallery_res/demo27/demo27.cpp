/*
`Ω` = [`e_1` `e_2`][`k_1` 0
		     0    `k_2`] [`e_1`^T
				  `e_2`^T]

where

`k_1`: ℝ  : control the desired kernel variance in either edge or orthogonal direction
`k_2`: ℝ  : control the desired kernel variance in either edge or orthogonal direction
`e_1`: ℝ ^ 3: orthogonal direction vectors
`e_2`: ℝ ^ 3: orthogonal direction vectors
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo27
 *
 * @param k_1  ℝ  : control the desired kernel variance in either edge or orthogonal direction
 * @param k_2  ℝ  : control the desired kernel variance in either edge or orthogonal direction
 * @param e_1  ℝ ^ 3: orthogonal direction vectors
 * @param e_2  ℝ ^ 3: orthogonal direction vectors
 * @return Ω
 */
Eigen::Matrix<double, 3, 3> demo27(
    const double & k_1,
    const double & k_2,
    const Eigen::Matrix<double, 3, 1> & e_1,
    const Eigen::Matrix<double, 3, 1> & e_2)
{
    Eigen::Matrix<double, 3, 2> _Ω_0;
    _Ω_0 << e_1, e_2;
    Eigen::Matrix<double, 2, 2> _Ω_1;
    _Ω_1 << k_1, 0,
    0, k_2;
    Eigen::Matrix<double, 2, 3> _Ω_2;
    _Ω_2 << e_1.transpose(),
    e_2.transpose();
    Eigen::Matrix<double, 3, 3> Ω = _Ω_0 * _Ω_1 * _Ω_2;

    return Ω;
}


void generateRandomData(double & k_1,
    double & k_2,
    Eigen::Matrix<double, 3, 1> & e_1,
    Eigen::Matrix<double, 3, 1> & e_2)
{
    k_1 = rand() % 10;
    k_2 = rand() % 10;
    e_1 = Eigen::VectorXd::Random(3);
    e_2 = Eigen::VectorXd::Random(3);
}


int main(int argc, char *argv[])
{
    double k_1;
    double k_2;
    Eigen::Matrix<double, 3, 1> e_1;
    Eigen::Matrix<double, 3, 1> e_2;
    generateRandomData(k_1, k_2, e_1, e_2);
    Eigen::Matrix<double, 3, 3> func_value = demo27(k_1, k_2, e_1, e_2);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}