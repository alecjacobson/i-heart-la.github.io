/*
min_(x ∈ ℝ^n) ∑_i ‖A_i x + b_i‖ + (1/2)‖x-`x₀`‖²

where

A_i ∈ ℝ^(m × n)  
`x₀` ∈ ℝ^n  
b_i ∈ ℝ^m  
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct convex_optimization_276ResultType {
    double ret;
    convex_optimization_276ResultType(const double & ret)
    : ret(ret)
    {}
};

convex_optimization_276ResultType convex_optimization_276(
    const std::vector<Eigen::MatrixXd> & A,
    const Eigen::VectorXd & x₀,
    const std::vector<Eigen::VectorXd> & b)
{
    const long dim_0 = A.size();
    const long m = A[0].rows();
    const long n = A[0].cols();
    assert( A.size() == dim_0 );
    for( const auto& el : A ) {
        assert( el.rows() == m );
        assert( el.cols() == n );
    }
    assert( x₀.size() == n );
    assert( b.size() == dim_0 );
    for( const auto& el : b ) {
        assert( el.size() == m );
    }

    double ret = ;
    return convex_optimization_276ResultType(ret);
}


void generateRandomData(std::vector<Eigen::MatrixXd> & A,
    Eigen::VectorXd & x₀,
    std::vector<Eigen::VectorXd> & b)
{
    const int dim_0 = rand()%10;
    const int m = rand()%10;
    const int n = rand()%10;
    A.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        A[i] = Eigen::MatrixXd::Random(m, n);
    }
    x₀ = Eigen::VectorXd::Random(n);
    b.resize(dim_0);
    for(int i=0; i<dim_0; i++){
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
    convex_optimization_276ResultType func_value = convex_optimization_276(A, x₀, b);
    std::cout<<"return value:\n"<<func_value.ret<<std::endl;
    return 0;
}