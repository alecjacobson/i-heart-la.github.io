/*
y_i = (a_i)ᵀ x + w_i
`x_bar` = (sum_i a_i(a_i)ᵀ)^(-1) sum_i y_i a_i

where

a_i: ℝ^n: the measurement vectors  
w_i: ℝ: measurement noise 
x: ℝ^n: measurement noise 
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo11
 *
 * @param a  ℝ^n: the measurement vectors  
 * @param w  ℝ: measurement noise 
 * @param x  ℝ^n: measurement noise 
 * @return x_bar
 */
Eigen::VectorXd demo11(
    const std::vector<Eigen::VectorXd> & a,
    const std::vector<double> & w,
    const Eigen::VectorXd & x)
{
    const long _dim_0 = w.size();
    const long n = x.size();
    assert( a.size() == _dim_0 );
    for( const auto& el : a ) {
        assert( el.size() == n );
    }
    assert( w.size() == _dim_0 );
    assert( x.size() == n );

    std::vector<double> y(_dim_0);
    for( int i=1; i<=_dim_0; i++){
        y.at(i-1) = (a.at(i-1)).transpose() * x + w.at(i-1);
    }


    Eigen::MatrixXd _sum_0 = Eigen::MatrixXd::Zero(n, n);
    for(int i=1; i<=a.size(); i++){
        _sum_0 += a.at(i-1) * (a.at(i-1)).transpose();
    }
    Eigen::MatrixXd _sum_1 = Eigen::MatrixXd::Zero(n, 1);
    for(int i=1; i<=a.size(); i++){
        _sum_1 += y.at(i-1) * a.at(i-1);
    }
    Eigen::VectorXd x_bar = (_sum_0).inverse() * _sum_1;

    return x_bar;
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
    std::vector<Eigen::VectorXd> a;
    std::vector<double> w;
    Eigen::VectorXd x;
    generateRandomData(a, w, x);
    Eigen::VectorXd func_value = demo11(a, w, x);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}