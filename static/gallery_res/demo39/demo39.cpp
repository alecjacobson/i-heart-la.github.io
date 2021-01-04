/*
v_i = sum_j w_i,j M_j u_i

where

w: ℝ^(4×4)
M_j: ℝ^(4×4)
u_i: ℝ^4
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo39
 *
 * @param w  ℝ^(4×4)
 * @param M  ℝ^(4×4)
 * @param u  ℝ^4
 * @return v
 */
std::vector<Eigen::Matrix<double, 4, 1>> demo39(
    const Eigen::Matrix<double, 4, 4> & w,
    const std::vector<Eigen::Matrix<double, 4, 4>> & M,
    const std::vector<Eigen::Matrix<double, 4, 1>> & u)
{
    const long _dim_0 = M.size();
    const long _dim_1 = u.size();
    assert( M.size() == _dim_0 );

    std::vector<Eigen::Matrix<double, 4, 1>> v(_dim_1);
    for( int i=1; i<=_dim_1; i++){
        Eigen::MatrixXd _sum_0 = Eigen::MatrixXd::Zero(4, 1);
        for(int j=1; j<=M.size(); j++){
            _sum_0 += w(i-1, j-1) * M.at(j-1) * u.at(i-1);
        }
        v.at(i-1) = _sum_0;
    }


    return v;
}


void generateRandomData(Eigen::Matrix<double, 4, 4> & w,
    std::vector<Eigen::Matrix<double, 4, 4>> & M,
    std::vector<Eigen::Matrix<double, 4, 1>> & u)
{
    const int _dim_0 = rand()%10;
    const int _dim_1 = rand()%10;
    w = Eigen::MatrixXd::Random(4, 4);
    M.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        M[i] = Eigen::MatrixXd::Random(4, 4);
    }
    u.resize(_dim_1);
    for(int i=0; i<_dim_1; i++){
        u[i] = Eigen::VectorXd::Random(4);
    }
}


int main(int argc, char *argv[])
{
    Eigen::Matrix<double, 4, 4> w;
    std::vector<Eigen::Matrix<double, 4, 4>> M;
    std::vector<Eigen::Matrix<double, 4, 1>> u;
    generateRandomData(w, M, u);
    std::vector<Eigen::Matrix<double, 4, 1>> func_value = demo39(w, M, u);
    std::cout<<"vector func_value:"<<std::endl;
    for(int i=0; i<func_value.size(); i++){
        std::cout<<"i:"<<i<<", value:\n"<<func_value.at(i)<<std::endl;
    }
    return 0;
}