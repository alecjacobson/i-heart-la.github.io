/*
[P₁  0
 0   P₃][    L        0
         P₃ᵀCP₂ᵀU⁻¹  -L̃][U  L⁻¹P₁ᵀB
                         0     Ũ   ][P₂   0
                                     0    I₄]

where

P_i ∈ ℝ^(4×4) 
B ∈ ℝ^(4×4) 
C ∈ ℝ^(4×4) 
L ∈ ℝ^(4×4) 
L̃ ∈ ℝ^(4×4) 
U ∈ ℝ^(4×4) 
Ũ ∈ ℝ^(4×4)
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct convex_optimization_680ResultType {
    Eigen::Matrix<double, 8, 8> ret;
    convex_optimization_680ResultType(const Eigen::Matrix<double, 8, 8> & ret)
    : ret(ret)
    {}
};

convex_optimization_680ResultType convex_optimization_680(
    const std::vector<Eigen::Matrix<double, 4, 4>> & P,
    const Eigen::Matrix<double, 4, 4> & B,
    const Eigen::Matrix<double, 4, 4> & C,
    const Eigen::Matrix<double, 4, 4> & L,
    const Eigen::Matrix<double, 4, 4> & L̃,
    const Eigen::Matrix<double, 4, 4> & U,
    const Eigen::Matrix<double, 4, 4> & Ũ)
{
    const long dim_0 = P.size();
    Eigen::Matrix<double, 8, 8> ret_0;
    ret_0 << P.at(1-1), Eigen::MatrixXd::Zero(4, 4),
    Eigen::MatrixXd::Zero(4, 4), P.at(3-1);
    Eigen::Matrix<double, 8, 8> ret_1;
    ret_1 << L, Eigen::MatrixXd::Zero(4, 4),
    P.at(3-1).transpose() * C * P.at(2-1).transpose() * U.inverse(), -L̃;
    Eigen::Matrix<double, 8, 8> ret_2;
    ret_2 << U, L.inverse() * P.at(1-1).transpose() * B,
    Eigen::MatrixXd::Zero(4, 4), Ũ;
    Eigen::Matrix<double, 8, 8> ret_3;
    ret_3 << P.at(2-1), Eigen::MatrixXd::Zero(4, 4),
    Eigen::MatrixXd::Zero(4, 4), Eigen::MatrixXd::Identity(4, 4);
    Eigen::Matrix<double, 8, 8> ret = ret_0 * ret_1 * ret_2 * ret_3;
    return convex_optimization_680ResultType(ret);
}


void generateRandomData(std::vector<Eigen::Matrix<double, 4, 4>> & P,
    Eigen::Matrix<double, 4, 4> & B,
    Eigen::Matrix<double, 4, 4> & C,
    Eigen::Matrix<double, 4, 4> & L,
    Eigen::Matrix<double, 4, 4> & L̃,
    Eigen::Matrix<double, 4, 4> & U,
    Eigen::Matrix<double, 4, 4> & Ũ)
{
    const int dim_0 = rand()%10;
    P.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        P[i] = Eigen::MatrixXd::Random(4, 4);
    }
    B = Eigen::MatrixXd::Random(4, 4);
    C = Eigen::MatrixXd::Random(4, 4);
    L = Eigen::MatrixXd::Random(4, 4);
    L̃ = Eigen::MatrixXd::Random(4, 4);
    U = Eigen::MatrixXd::Random(4, 4);
    Ũ = Eigen::MatrixXd::Random(4, 4);
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    std::vector<Eigen::Matrix<double, 4, 4>> P;
    Eigen::Matrix<double, 4, 4> B;
    Eigen::Matrix<double, 4, 4> C;
    Eigen::Matrix<double, 4, 4> L;
    Eigen::Matrix<double, 4, 4> L̃;
    Eigen::Matrix<double, 4, 4> U;
    Eigen::Matrix<double, 4, 4> Ũ;
    generateRandomData(P, B, C, L, L̃, U, Ũ);
    convex_optimization_680ResultType func_value = convex_optimization_680(P, B, C, L, L̃, U, Ũ);
    std::cout<<"return value:\n"<<func_value.ret<<std::endl;
    return 0;
}