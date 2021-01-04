/*
`I(X;Y)` = ∑_i ∑_j x_j p_i,j log(p_ij/∑_k x_k p_ik)

where

x: ℝ^n
p: ℝ^(n×m)

*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo8
 *
 * @param x  ℝ^n
 * @param p  ℝ^(n×m)
 * @return I_left_parenthesis_X_semicolon_Y_right_parenthesis
 */
double demo8(
    const Eigen::VectorXd & x,
    const Eigen::MatrixXd & p)
{
    const long n = p.rows();
    const long m = p.cols();
    assert( x.size() == n );
    assert( p.rows() == n );
    assert( p.cols() == m );

    double _sum_0 = 0;
    for(int i=1; i<=p.rows(); i++){
        double _sum_1 = 0;
        for(int j=1; j<=x.size(); j++){
            double _sum_2 = 0;
            for(int k=1; k<=p.cols(); k++){
                _sum_2 += x(k-1) * p(i-1, k-1);
            }
            _sum_1 += x(j-1) * p(i-1, j-1) * log10(p(i-1, j-1) / _sum_2);
        }
        _sum_0 += _sum_1;
    }
    double I_left_parenthesis_X_semicolon_Y_right_parenthesis = _sum_0;

    return I_left_parenthesis_X_semicolon_Y_right_parenthesis;
}


void generateRandomData(Eigen::VectorXd & x,
    Eigen::MatrixXd & p)
{
    const int n = rand()%10;
    const int m = rand()%10;
    x = Eigen::VectorXd::Random(n);
    p = Eigen::MatrixXd::Random(n, m);
}


int main(int argc, char *argv[])
{
    Eigen::VectorXd x;
    Eigen::MatrixXd p;
    generateRandomData(x, p);
    double func_value = demo8(x, p);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}