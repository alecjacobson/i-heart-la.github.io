/*
y_i = a_iᵀ x + w_i
x̂ = (∑_i a_i a_iᵀ)⁻¹ ∑_i y_i a_i

where

a_i ∈ ℝ^n: the measurement vectors  
w_i ∈ ℝ: measurement noise 
x ∈ ℝ^n: original vector 
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * convex_optimization_384
 *
 * @param a  the measurement vectors  
 * @param w  measurement noise 
 * @param x  original vector 
 * @return x̂
 */
Eigen::VectorXd convex_optimization_384(
    const std::vector<Eigen::VectorXd> & a,
    const std::vector<double> & w,
    const Eigen::VectorXd & x)
{
    const long _dim_0 = w.size();
    const long n = a[0].rows();
    assert( a.size() == _dim_0 );
    for( const auto& el : a ) {
        assert( el.size() == n );
    }
    assert( w.size() == _dim_0 );
    assert( x.size() == n );

    std::vector<double> y(_dim_0);
    for( int i=1; i<=_dim_0; i++){
        y.at(i-1) = (double)(a.at(i-1).transpose() * x) + w.at(i-1);
    }


    Eigen::MatrixXd _sum_0 = Eigen::MatrixXd::Zero(n, n);
    for(int i=1; i<=a.size(); i++){
        _sum_0 += a.at(i-1) * a.at(i-1).transpose();
    }
    Eigen::MatrixXd _sum_1 = Eigen::MatrixXd::Zero(n, 1);
    for(int i=1; i<=a.size(); i++){
        _sum_1 += y.at(i-1) * a.at(i-1);
    }
    Eigen::VectorXd x̂ = (_sum_0).inverse() * _sum_1;

    return x̂;
}


void generateRandomData(std::vector<Eigen::VectorXd> & a,
    std::vector<double> & w,
    Eigen::VectorXd & x)
{
    const int _dim_0 = rand()%10;
    const int n = rand()%10;
    a.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        a[i] = Eigen::VectorXd::Random(n);
    }
    w.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        w[i] = rand() % 10;
    }
    x = Eigen::VectorXd::Random(n);
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    std::vector<Eigen::VectorXd> a;
    std::vector<double> w;
    Eigen::VectorXd x;
    generateRandomData(a, w, x);
    Eigen::VectorXd func_value = convex_optimization_384(a, w, x);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}