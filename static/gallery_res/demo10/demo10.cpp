/*
min_(x ∈ ℝ^n) ∑_i ||A_i x + b_i ||_2 +(1/2)||x-`x₀`||^2_2

where

A_i: ℝ^(m × n)  
`x₀`: ℝ^n  
b_i: ℝ^m  
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo10
 *
 * @param A  ℝ^(m × n)
 * @param x₀  ℝ^n
 * @param b  ℝ^m
 * @return ret
 */
double demo10(
    const std::vector<Eigen::MatrixXd> & A,
    const Eigen::VectorXd & x₀,
    const std::vector<Eigen::VectorXd> & b)
{
    const long _dim_0 = A.size();
    const long m = b[0].rows();
    const long n = x₀.size();
    assert( A.size() == _dim_0 );
    for( const auto& el : A ) {
        assert( el.rows() == m );
        assert( el.cols() == n );
    }
    assert( x₀.size() == n );
    assert( b.size() == _dim_0 );
    for( const auto& el : b ) {
        assert( el.size() == m );
    }

    double ret = ;
    return ret;
}


void generateRandomData(std::vector<Eigen::MatrixXd> & A,
    Eigen::VectorXd & x₀,
    std::vector<Eigen::VectorXd> & b)
{
    const int _dim_0 = rand()%10;
    const int m = rand()%10;
    const int n = rand()%10;
    A.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        A[i] = Eigen::MatrixXd::Random(m, n);
    }
    x₀ = Eigen::VectorXd::Random(n);
    b.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        b[i] = Eigen::VectorXd::Random(m);
    }
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    std::vector<Eigen::MatrixXd> A;
    Eigen::VectorXd x₀;
    std::vector<Eigen::VectorXd> b;
    generateRandomData(A, x₀, b);
    double func_value = demo10(A, x₀, b);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}