/*
y_i = a_iᵀ x + w_i
x̂ = (∑_i a_i a_iᵀ)⁻¹ ∑_i y_i a_i

where

a_i ∈ ℝ^m: the measurement vectors  
w_i ∈ ℝ: measurement noise 
x ∈ ℝ^m: original vector 
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct convex_optimization_384ResultType {
    Eigen::VectorXd y;
    Eigen::VectorXd x̂;
    convex_optimization_384ResultType(const Eigen::VectorXd & y,
               const Eigen::VectorXd & x̂)
    : y(y),
    x̂(x̂)
    {}
};

/**
 * convex_optimization_384
 *
 * @param a  the measurement vectors  
 * @param w  measurement noise 
 * @param x  original vector 
 * @return x̂
 */
convex_optimization_384ResultType convex_optimization_384(
    const std::vector<Eigen::VectorXd> & a,
    const std::vector<double> & w,
    const Eigen::VectorXd & x)
{
    const long dim_0 = w.size();
    const long m = a[0].rows();
    assert( a.size() == dim_0 );
    for( const auto& el : a ) {
        assert( el.size() == m );
    }
    assert( x.size() == m );

    Eigen::VectorXd y(dim_0);
    for( int i=1; i<=dim_0; i++){
        y[i-1] = (double)(a.at(i-1).transpose() * x) + w.at(i-1);
    }

    Eigen::MatrixXd sum_0 = Eigen::MatrixXd::Zero(m, m);
    for(int i=1; i<=a.size(); i++){
        sum_0 += a.at(i-1) * a.at(i-1).transpose();
    }
    Eigen::MatrixXd sum_1 = Eigen::MatrixXd::Zero(m, 1);
    for(int i=1; i<=y.size(); i++){
        sum_1 += y[i-1] * a.at(i-1);
    }
    Eigen::VectorXd x̂ = (sum_0).inverse() * sum_1;

    return convex_optimization_384ResultType(y, x̂);
}


void generateRandomData(std::vector<Eigen::VectorXd> & a,
    std::vector<double> & w,
    Eigen::VectorXd & x)
{
    const int dim_0 = rand()%10;
    const int m = rand()%10;
    a.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        a[i] = Eigen::VectorXd::Random(m);
    }
    w.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        w[i] = rand() % 10;
    }
    x = Eigen::VectorXd::Random(m);
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    std::vector<Eigen::VectorXd> a;
    std::vector<double> w;
    Eigen::VectorXd x;
    generateRandomData(a, w, x);
    convex_optimization_384ResultType func_value = convex_optimization_384(a, w, x);
    std::cout<<"return value:\n"<<func_value.x̂<<std::endl;
    return 0;
}