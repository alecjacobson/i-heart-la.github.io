/*
`T₁` = 1/√2 U[0 0 0
              0 0 -1
              0 1 0]Vᵀ

where 

U ∈ ℝ^(3×3) 
V ∈ ℝ^(3×3)
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct analytic_eigensystems_13ResultType {
    Eigen::Matrix<double, 3, 3> T₁;
    analytic_eigensystems_13ResultType(const Eigen::Matrix<double, 3, 3> & T₁)
    : T₁(T₁)
    {}
};

analytic_eigensystems_13ResultType analytic_eigensystems_13(
    const Eigen::Matrix<double, 3, 3> & U,
    const Eigen::Matrix<double, 3, 3> & V)
{
    Eigen::Matrix<double, 3, 3> T₁_0;
    T₁_0 << 0, 0, 0,
    0, 0, -1,
    0, 1, 0;
    Eigen::Matrix<double, 3, 3> T₁ = 1 / double(sqrt(2)) * U * T₁_0 * V.transpose();

    return analytic_eigensystems_13ResultType(T₁);
}


void generateRandomData(Eigen::Matrix<double, 3, 3> & U,
    Eigen::Matrix<double, 3, 3> & V)
{
    U = Eigen::MatrixXd::Random(3, 3);
    V = Eigen::MatrixXd::Random(3, 3);
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    Eigen::Matrix<double, 3, 3> U;
    Eigen::Matrix<double, 3, 3> V;
    generateRandomData(U, V);
    analytic_eigensystems_13ResultType func_value = analytic_eigensystems_13(U, V);
    std::cout<<"return value:\n"<<func_value.T₁<<std::endl;
    return 0;
}