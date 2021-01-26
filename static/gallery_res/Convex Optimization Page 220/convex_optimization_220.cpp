/*
`L(x,v)` = xᵀWx + ∑_i v_i(x_i²-1)

where

x ∈ ℝ^n
W ∈ ℝ^(n×n)
v ∈ ℝ^n
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

double convex_optimization_220(
    const Eigen::VectorXd & x,
    const Eigen::MatrixXd & W,
    const Eigen::VectorXd & v)
{
    const long n = x.size();
    assert( x.size() == n );
    assert( W.rows() == n );
    assert( W.cols() == n );
    assert( v.size() == n );

    double _sum_0 = 0;
    for(int i=1; i<=v.size(); i++){
        _sum_0 += v[i-1] * (pow(x[i-1], 2) - 1);
    }
    double L_left_parenthesis_x_comma_v_right_parenthesis = (double)(x.transpose() * W * x) + _sum_0;

    return L_left_parenthesis_x_comma_v_right_parenthesis;
}


void generateRandomData(Eigen::VectorXd & x,
    Eigen::MatrixXd & W,
    Eigen::VectorXd & v)
{
    const int n = rand()%10;
    x = Eigen::VectorXd::Random(n);
    W = Eigen::MatrixXd::Random(n, n);
    v = Eigen::VectorXd::Random(n);
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    Eigen::VectorXd x;
    Eigen::MatrixXd W;
    Eigen::VectorXd v;
    generateRandomData(x, W, v);
    double func_value = convex_optimization_220(x, W, v);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}