/*
`L(x,ν)` = xᵀWx + ∑_i ν_i(x_i²-1)

where

x ∈ ℝ^n
W ∈ ℝ^(n×n)
ν ∈ ℝ^n
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct convex_optimization_220ResultType {
    double L_left_parenthesis_x_comma_ν_right_parenthesis;
    convex_optimization_220ResultType(const double & L_left_parenthesis_x_comma_ν_right_parenthesis)
    : L_left_parenthesis_x_comma_ν_right_parenthesis(L_left_parenthesis_x_comma_ν_right_parenthesis)
    {}
};

convex_optimization_220ResultType convex_optimization_220(
    const Eigen::VectorXd & x,
    const Eigen::MatrixXd & W,
    const Eigen::VectorXd & ν)
{
    const long n = x.size();
    assert( W.rows() == n );
    assert( W.cols() == n );
    assert( ν.size() == n );

    double sum_0 = 0;
    for(int i=1; i<=ν.size(); i++){
        sum_0 += ν[i-1] * (pow(x[i-1], 2) - 1);
    }
    double L_left_parenthesis_x_comma_ν_right_parenthesis = (double)(x.transpose() * W * x) + sum_0;

    return convex_optimization_220ResultType(L_left_parenthesis_x_comma_ν_right_parenthesis);
}


void generateRandomData(Eigen::VectorXd & x,
    Eigen::MatrixXd & W,
    Eigen::VectorXd & ν)
{
    const int n = rand()%10;
    x = Eigen::VectorXd::Random(n);
    W = Eigen::MatrixXd::Random(n, n);
    ν = Eigen::VectorXd::Random(n);
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    Eigen::VectorXd x;
    Eigen::MatrixXd W;
    Eigen::VectorXd ν;
    generateRandomData(x, W, ν);
    convex_optimization_220ResultType func_value = convex_optimization_220(x, W, ν);
    std::cout<<"return value:\n"<<func_value.L_left_parenthesis_x_comma_ν_right_parenthesis<<std::endl;
    return 0;
}