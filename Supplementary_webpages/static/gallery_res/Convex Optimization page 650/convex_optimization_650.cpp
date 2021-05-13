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

struct convex_optimization_650ResultType {
    Eigen::Matrix<double, 8, 8> ret;
    convex_optimization_650ResultType(const Eigen::Matrix<double, 8, 8> & ret)
    : ret(ret)
    {}
};

convex_optimization_650ResultType convex_optimization_650(
    const Eigen::Matrix<double, 4, 4> & A,
    const Eigen::Matrix<double, 4, 4> & B,
    const Eigen::Matrix<double, 4, 4> & S)
{
    Eigen::Matrix<double, 8, 8> ret_0;
    ret_0 << A.inverse() + A.inverse() * B * S.inverse() * B.transpose() * A.inverse(), -A.inverse() * B * S.inverse(),
    -S.inverse() * B.transpose() * A.inverse(), S.inverse();
    Eigen::Matrix<double, 8, 8> ret = ret_0;
    return convex_optimization_650ResultType(ret);
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
    convex_optimization_650ResultType func_value = convex_optimization_650(A, B, S);
    std::cout<<"return value:\n"<<func_value.ret<<std::endl;
    return 0;
}