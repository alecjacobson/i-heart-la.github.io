/*
given
p_i: ℝ^3: points on lines
d_i: ℝ^3: unit directions along lines

P_i = ( I_3 - d_i d_iᵀ )
q = ( ∑_i P_iᵀP_i )⁻¹ ( ∑_i P_iᵀP_i p_i )
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * least_square
 *
 * @param p  ℝ^3: points on lines
 * @param d  ℝ^3: unit directions along lines
 * @return q
 */
Eigen::Matrix<double, 3, 1> least_square(
    const std::vector<Eigen::Matrix<double, 3, 1>> & p,
    const std::vector<Eigen::Matrix<double, 3, 1>> & d)
{
    const long _dim_0 = p.size();
    std::vector<Eigen::Matrix<double, 3, 3>> P(_dim_0);
    for( int i=1; i<=_dim_0; i++){
        P.at(i-1) = (Eigen::MatrixXd::Identity(3, 3) - d.at(i-1) * d.at(i-1).transpose());
    }


    Eigen::MatrixXd _sum_0 = Eigen::MatrixXd::Zero(3, 3);
    for(int i=1; i<=P.size(); i++){
        _sum_0 += P.at(i-1).transpose() * P.at(i-1);
    }
    Eigen::MatrixXd _sum_1 = Eigen::MatrixXd::Zero(3, 1);
    for(int i=1; i<=p.size(); i++){
        _sum_1 += P.at(i-1).transpose() * P.at(i-1) * p.at(i-1);
    }
    Eigen::Matrix<double, 3, 1> q = (_sum_0).inverse() * (_sum_1);

    return q;
}


void generateRandomData(std::vector<Eigen::Matrix<double, 3, 1>> & p,
    std::vector<Eigen::Matrix<double, 3, 1>> & d)
{
    const int _dim_0 = rand()%10;
    p.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        p[i] = Eigen::VectorXd::Random(3);
    }
    d.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        d[i] = Eigen::VectorXd::Random(3);
    }
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    std::vector<Eigen::Matrix<double, 3, 1>> p;
    std::vector<Eigen::Matrix<double, 3, 1>> d;
    generateRandomData(p, d);
    Eigen::Matrix<double, 3, 1> func_value = least_square(p, d);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}