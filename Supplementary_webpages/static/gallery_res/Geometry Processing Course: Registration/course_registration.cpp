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

struct course_registrationResultType {
    double ret;
    course_registrationResultType(const double & ret)
    : ret(ret)
    {}
};

/**
 * course_registration
 *
 * @param x  ℝ^3
 * @param n̂  ℝ^3
 * @param p  ℝ^3
 * @return ret
 */
course_registrationResultType course_registration(
    const std::vector<Eigen::Matrix<double, 3, 1>> & x,
    const std::vector<Eigen::Matrix<double, 3, 1>> & n̂,
    const std::vector<Eigen::Matrix<double, 3, 1>> & p)
{
    const long dim_0 = x.size();
    assert( n̂.size() == dim_0 );
    assert( p.size() == dim_0 );

    double ret = ;
    return course_registrationResultType(ret);
}


void generateRandomData(std::vector<Eigen::Matrix<double, 3, 1>> & x,
    std::vector<Eigen::Matrix<double, 3, 1>> & n̂,
    std::vector<Eigen::Matrix<double, 3, 1>> & p)
{
    const int dim_0 = rand()%10;
    x.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        x[i] = Eigen::VectorXd::Random(3);
    }
    n̂.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        n̂[i] = Eigen::VectorXd::Random(3);
    }
    p.resize(dim_0);
    for(int i=0; i<dim_0; i++){
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
    course_registrationResultType func_value = course_registration(x, n̂, p);
    std::cout<<"return value:\n"<<func_value.ret<<std::endl;
    return 0;
}