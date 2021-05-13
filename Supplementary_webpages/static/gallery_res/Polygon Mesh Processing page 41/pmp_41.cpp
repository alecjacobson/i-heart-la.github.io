/*
`xᵢ` = T_*,1
`xⱼ` = T_*,2
`xₖ` = T_*,3
`n(T)` = (`xⱼ`-`xᵢ`)×(`xₖ`-`xᵢ`)/‖(`xⱼ`-`xᵢ`)×(`xₖ`-`xᵢ`)‖

where
 
T ∈ ℝ^(3×3)
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct pmp_41ResultType {
    Eigen::Matrix<double, 3, 1> xᵢ;
    Eigen::Matrix<double, 3, 1> xⱼ;
    Eigen::Matrix<double, 3, 1> xₖ;
    Eigen::Matrix<double, 3, 1> n_left_parenthesis_T_right_parenthesis;
    pmp_41ResultType(const Eigen::Matrix<double, 3, 1> & xᵢ,
               const Eigen::Matrix<double, 3, 1> & xⱼ,
               const Eigen::Matrix<double, 3, 1> & xₖ,
               const Eigen::Matrix<double, 3, 1> & n_left_parenthesis_T_right_parenthesis)
    : xᵢ(xᵢ),
    xⱼ(xⱼ),
    xₖ(xₖ),
    n_left_parenthesis_T_right_parenthesis(n_left_parenthesis_T_right_parenthesis)
    {}
};

pmp_41ResultType pmp_41(const Eigen::Matrix<double, 3, 3> & T)
{
    Eigen::Matrix<double, 3, 1> xᵢ = T.col(1-1);

    Eigen::Matrix<double, 3, 1> xⱼ = T.col(2-1);

    Eigen::Matrix<double, 3, 1> xₖ = T.col(3-1);

    Eigen::Matrix<double, 3, 1> n_left_parenthesis_T_right_parenthesis = ((xⱼ - xᵢ)).cross((xₖ - xᵢ)) / double((((xⱼ - xᵢ)).cross((xₖ - xᵢ))).lpNorm<2>());

    return pmp_41ResultType(xᵢ, xⱼ, xₖ, n_left_parenthesis_T_right_parenthesis);
}


void generateRandomData(Eigen::Matrix<double, 3, 3> & T)
{
    T = Eigen::MatrixXd::Random(3, 3);
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    Eigen::Matrix<double, 3, 3> T;
    generateRandomData(T);
    pmp_41ResultType func_value = pmp_41(T);
    std::cout<<"return value:\n"<<func_value.n_left_parenthesis_T_right_parenthesis<<std::endl;
    return 0;
}