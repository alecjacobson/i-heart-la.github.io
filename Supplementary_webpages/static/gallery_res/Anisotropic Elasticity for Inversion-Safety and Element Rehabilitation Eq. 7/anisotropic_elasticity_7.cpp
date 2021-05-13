/*
`∂²I₅/∂f²` = 2[A₁,₁I₃  A₁,₂I₃  A₁,₃I₃
               A₂,₁I₃  A₂,₂I₃  A₂,₃I₃
               A₃,₁I₃  A₃,₂I₃  A₃,₃I₃] 

where

A ∈ ℝ^(3×3)
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct anisotropic_elasticity_7ResultType {
    Eigen::Matrix<double, 9, 9> _partial_differential_²I₅_soliduspartial_differential_f²;
    anisotropic_elasticity_7ResultType(const Eigen::Matrix<double, 9, 9> & _partial_differential_²I₅_soliduspartial_differential_f²)
    : _partial_differential_²I₅_soliduspartial_differential_f²(_partial_differential_²I₅_soliduspartial_differential_f²)
    {}
};

anisotropic_elasticity_7ResultType anisotropic_elasticity_7(const Eigen::Matrix<double, 3, 3> & A)
{
    Eigen::Matrix<double, 9, 9> _partial_differential_²I₅_soliduspartial_differential_f²_0;
    _partial_differential_²I₅_soliduspartial_differential_f²_0 << A(1-1, 1-1) * Eigen::MatrixXd::Identity(3, 3), A(1-1, 2-1) * Eigen::MatrixXd::Identity(3, 3), A(1-1, 3-1) * Eigen::MatrixXd::Identity(3, 3),
    A(2-1, 1-1) * Eigen::MatrixXd::Identity(3, 3), A(2-1, 2-1) * Eigen::MatrixXd::Identity(3, 3), A(2-1, 3-1) * Eigen::MatrixXd::Identity(3, 3),
    A(3-1, 1-1) * Eigen::MatrixXd::Identity(3, 3), A(3-1, 2-1) * Eigen::MatrixXd::Identity(3, 3), A(3-1, 3-1) * Eigen::MatrixXd::Identity(3, 3);
    Eigen::Matrix<double, 9, 9> _partial_differential_²I₅_soliduspartial_differential_f² = 2 * _partial_differential_²I₅_soliduspartial_differential_f²_0;

    return anisotropic_elasticity_7ResultType(_partial_differential_²I₅_soliduspartial_differential_f²);
}


void generateRandomData(Eigen::Matrix<double, 3, 3> & A)
{
    A = Eigen::MatrixXd::Random(3, 3);
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    Eigen::Matrix<double, 3, 3> A;
    generateRandomData(A);
    anisotropic_elasticity_7ResultType func_value = anisotropic_elasticity_7(A);
    std::cout<<"return value:\n"<<func_value._partial_differential_²I₅_soliduspartial_differential_f²<<std::endl;
    return 0;
}