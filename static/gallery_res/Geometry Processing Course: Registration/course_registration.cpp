/*
min_(u ∈ ℝ^6) uᵀ(∑_i [x_i×n̂_i
                        n̂_i  ][(x_i×n̂_i)ᵀ n̂_iᵀ])u - 2uᵀ(∑_i [x_i×n̂_i
                                                               n̂_i  ]n̂_iᵀ(p_i-x_i)) + ∑_i(p_i-x_i)ᵀn̂_i n̂_iᵀ(p_i-x_i)

where

x_i: ℝ^3 
n̂_i: ℝ^3  
p_i: ℝ^3  
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * course_registration
 *
 * @param x  ℝ^3
 * @param n̂  ℝ^3
 * @param p  ℝ^3
 * @return ret
 */
double course_registration(
    const std::vector<Eigen::Matrix<double, 3, 1>> & x,
    const std::vector<Eigen::Matrix<double, 3, 1>> & n̂,
    const std::vector<Eigen::Matrix<double, 3, 1>> & p)
{
    const long _dim_0 = x.size();
    double ret = ;
    return ret;
}


void generateRandomData(std::vector<Eigen::Matrix<double, 3, 1>> & x,
    std::vector<Eigen::Matrix<double, 3, 1>> & n̂,
    std::vector<Eigen::Matrix<double, 3, 1>> & p)
{
    const int _dim_0 = rand()%10;
    x.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        x[i] = Eigen::VectorXd::Random(3);
    }
    n̂.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        n̂[i] = Eigen::VectorXd::Random(3);
    }
    p.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        p[i] = Eigen::VectorXd::Random(3);
    }
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    std::vector<Eigen::Matrix<double, 3, 1>> x;
    std::vector<Eigen::Matrix<double, 3, 1>> n̂;
    std::vector<Eigen::Matrix<double, 3, 1>> p;
    generateRandomData(x, n̂, p);
    double func_value = course_registration(x, n̂, p);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}