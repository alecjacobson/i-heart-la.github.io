/*
min_(u ∈ ℝ^6) uᵀ(sum_i [x_i×n_i; n_i][(x_i×n_i)ᵀ n_iᵀ])u - 2uᵀ(sum_i [x_i×n_i; n_i]n_iᵀ(p_i-x_i)) + sum_i(p_i-x_i)ᵀn_i n_iᵀ(p_i-x_i)

where

x_i: ℝ^3 
n_i: ℝ^3  
p_i: ℝ^3  
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo25
 *
 * @param x  ℝ^3
 * @param n  ℝ^3
 * @param p  ℝ^3
 * @return ret
 */
double demo25(
    const std::vector<Eigen::Matrix<double, 3, 1>> & x,
    const std::vector<Eigen::Matrix<double, 3, 1>> & n,
    const std::vector<Eigen::Matrix<double, 3, 1>> & p)
{
    const long _dim_0 = x.size();
    double ret = ;
    return ret;
}


void generateRandomData(std::vector<Eigen::Matrix<double, 3, 1>> & x,
    std::vector<Eigen::Matrix<double, 3, 1>> & n,
    std::vector<Eigen::Matrix<double, 3, 1>> & p)
{
    const int _dim_0 = rand()%10;
    x.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        x[i] = Eigen::VectorXd::Random(3);
    }
    n.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        n[i] = Eigen::VectorXd::Random(3);
    }
    p.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        p[i] = Eigen::VectorXd::Random(3);
    }
}


int main(int argc, char *argv[])
{
    std::vector<Eigen::Matrix<double, 3, 1>> x;
    std::vector<Eigen::Matrix<double, 3, 1>> n;
    std::vector<Eigen::Matrix<double, 3, 1>> p;
    generateRandomData(x, n, p);
    double func_value = demo25(x, n, p);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}