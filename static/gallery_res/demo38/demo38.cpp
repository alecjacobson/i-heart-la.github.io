/*
min_(C ∈ ℝ^3) sum_i ||x_i + (R_i - I_3)C ||²

where

x_i: ℝ^3
R_i: ℝ^(3×3) 
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo38
 *
 * @param x  ℝ^3
 * @param R  ℝ^(3×3)
 * @return ret
 */
double demo38(
    const std::vector<Eigen::Matrix<double, 3, 1>> & x,
    const std::vector<Eigen::Matrix<double, 3, 3>> & R)
{
    const long _dim_0 = x.size();
    assert( R.size() == _dim_0 );

    double ret = ;
    return ret;
}


void generateRandomData(std::vector<Eigen::Matrix<double, 3, 1>> & x,
    std::vector<Eigen::Matrix<double, 3, 3>> & R)
{
    const int _dim_0 = rand()%10;
    x.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        x[i] = Eigen::VectorXd::Random(3);
    }
    R.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        R[i] = Eigen::MatrixXd::Random(3, 3);
    }
}


int main(int argc, char *argv[])
{
    std::vector<Eigen::Matrix<double, 3, 1>> x;
    std::vector<Eigen::Matrix<double, 3, 3>> R;
    generateRandomData(x, R);
    double func_value = demo38(x, R);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}