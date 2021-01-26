/*
from trigonometry: cos

∑_i cos(θ)²((p_i - q_i)⋅n_i +((p_i+q_i)×n_i)⋅ã+n_i⋅t̃)² 

where

θ ∈ ℝ: angle of rotation
p_i ∈ ℝ^3
q_i ∈ ℝ^3
n_i ∈ ℝ^3
ã ∈ ℝ^3
t̃ ∈ ℝ^3
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * symmetric_objective_function_9
 *
 * @param θ  angle of rotation
 * @return ret
 */
double symmetric_objective_function_9(
    const double & θ,
    const std::vector<Eigen::Matrix<double, 3, 1>> & p,
    const std::vector<Eigen::Matrix<double, 3, 1>> & q,
    const std::vector<Eigen::Matrix<double, 3, 1>> & n,
    const Eigen::Matrix<double, 3, 1> & ã,
    const Eigen::Matrix<double, 3, 1> & t̃)
{
    const long _dim_0 = p.size();
    double _sum_0 = 0;
    for(int i=1; i<=n.size(); i++){
        _sum_0 += pow(cos(θ), 2) * pow((((p.at(i-1) - q.at(i-1))).dot(n.at(i-1)) + ((((p.at(i-1) + q.at(i-1))).cross(n.at(i-1)))).dot(ã) + (n.at(i-1)).dot(t̃)), 2);
    }
    double ret = _sum_0;
    return ret;
}


void generateRandomData(double & θ,
    std::vector<Eigen::Matrix<double, 3, 1>> & p,
    std::vector<Eigen::Matrix<double, 3, 1>> & q,
    std::vector<Eigen::Matrix<double, 3, 1>> & n,
    Eigen::Matrix<double, 3, 1> & ã,
    Eigen::Matrix<double, 3, 1> & t̃)
{
    θ = rand() % 10;
    const int _dim_0 = rand()%10;
    p.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        p[i] = Eigen::VectorXd::Random(3);
    }
    q.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        q[i] = Eigen::VectorXd::Random(3);
    }
    n.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        n[i] = Eigen::VectorXd::Random(3);
    }
    ã = Eigen::VectorXd::Random(3);
    t̃ = Eigen::VectorXd::Random(3);
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    double θ;
    std::vector<Eigen::Matrix<double, 3, 1>> p;
    std::vector<Eigen::Matrix<double, 3, 1>> q;
    std::vector<Eigen::Matrix<double, 3, 1>> n;
    Eigen::Matrix<double, 3, 1> ã;
    Eigen::Matrix<double, 3, 1> t̃;
    generateRandomData(θ, p, q, n, ã, t̃);
    double func_value = symmetric_objective_function_9(θ, p, q, n, ã, t̃);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}