/*
`I(X;Y)` = ∑_i ∑_j x_j p_i,j log₂(p_i,j/∑_k x_k p_i,k)

where

x ∈ ℝ^n
p ∈ ℝ^(m×n)
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct convex_optimization_208ResultType {
    double I_left_parenthesis_X_semicolon_Y_right_parenthesis;
    convex_optimization_208ResultType(const double & I_left_parenthesis_X_semicolon_Y_right_parenthesis)
    : I_left_parenthesis_X_semicolon_Y_right_parenthesis(I_left_parenthesis_X_semicolon_Y_right_parenthesis)
    {}
};

convex_optimization_208ResultType convex_optimization_208(
    const Eigen::VectorXd & x,
    const Eigen::MatrixXd & p)
{
    const long n = x.size();
    const long m = p.rows();
    assert( p.cols() == n );

    double sum_0 = 0;
    for(int i=1; i<=p.rows(); i++){
        double sum_1 = 0;
        for(int j=1; j<=x.size(); j++){
            double sum_2 = 0;
            for(int k=1; k<=x.size(); k++){
                sum_2 += x[k-1] * p(i-1, k-1);
            }
            sum_1 += x[j-1] * p(i-1, j-1) * (log10(p(i-1, j-1) / double(sum_2)) / log10(2));
        }
        sum_0 += sum_1;
    }
    double I_left_parenthesis_X_semicolon_Y_right_parenthesis = sum_0;

    return convex_optimization_208ResultType(I_left_parenthesis_X_semicolon_Y_right_parenthesis);
}


void generateRandomData(Eigen::VectorXd & x,
    Eigen::MatrixXd & p)
{
    const int n = rand()%10;
    const int m = rand()%10;
    x = Eigen::VectorXd::Random(n);
    p = Eigen::MatrixXd::Random(m, n);
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    Eigen::VectorXd x;
    Eigen::MatrixXd p;
    generateRandomData(x, p);
    convex_optimization_208ResultType func_value = convex_optimization_208(x, p);
    std::cout<<"return value:\n"<<func_value.I_left_parenthesis_X_semicolon_Y_right_parenthesis<<std::endl;
    return 0;
}