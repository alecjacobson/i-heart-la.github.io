/*
[A⁻¹+A⁻¹BS⁻¹BᵀA⁻¹   -A⁻¹BS⁻¹
    -S⁻¹BᵀA⁻¹           S⁻¹]

where

A ∈ ℝ^(4×4) 
B ∈ ℝ^(4×4) 
S ∈ ℝ^(4×4) 
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

Eigen::Matrix<double, 8, 8> demo14(
    const Eigen::Matrix<double, 4, 4> & A,
    const Eigen::Matrix<double, 4, 4> & B,
    const Eigen::Matrix<double, 4, 4> & S)
{
    Eigen::Matrix<double, 8, 8> _ret_0;
    _ret_0 << A.inverse() + A.inverse() * B * S.inverse() * B.transpose() * A.inverse(), -A.inverse() * B * S.inverse(),
    -S.inverse() * B.transpose() * A.inverse(), S.inverse();
    Eigen::Matrix<double, 8, 8> ret = _ret_0;
    return ret;
}


void generateRandomData(Eigen::Matrix<double, 4, 4> & A,
    Eigen::Matrix<double, 4, 4> & B,
    Eigen::Matrix<double, 4, 4> & S)
{
    A = Eigen::MatrixXd::Random(4, 4);
    B = Eigen::MatrixXd::Random(4, 4);
    S = Eigen::MatrixXd::Random(4, 4);
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    Eigen::Matrix<double, 4, 4> A;
    Eigen::Matrix<double, 4, 4> B;
    Eigen::Matrix<double, 4, 4> S;
    generateRandomData(A, B, S);
    Eigen::Matrix<double, 8, 8> func_value = demo14(A, B, S);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}