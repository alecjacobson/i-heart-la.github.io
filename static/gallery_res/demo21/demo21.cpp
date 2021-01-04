/*
from trigonometry: cos

b = ∑_i cos(θ)²((p_i - q_i)⋅n_i +((p_i+q_i)×n_i)⋅a+n_i⋅t)² 

where

θ: ℝ: angle of rotation
p_i: ℝ^3
q_i: ℝ^3
n_i: ℝ^3
a: ℝ^3
t: ℝ^3



*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo21
 *
 * @param θ  ℝ: angle of rotation
 * @param p  ℝ^3
 * @param q  ℝ^3
 * @param n  ℝ^3
 * @param a  ℝ^3
 * @param t  ℝ^3
 * @return b
 */
double demo21(
    const double & θ,
    const std::vector<Eigen::Matrix<double, 3, 1>> & p,
    const std::vector<Eigen::Matrix<double, 3, 1>> & q,
    const std::vector<Eigen::Matrix<double, 3, 1>> & n,
    const Eigen::Matrix<double, 3, 1> & a,
    const Eigen::Matrix<double, 3, 1> & t)
{
    const long _dim_0 = p.size();
    double _sum_0 = 0;
    for(int i=1; i<=q.size(); i++){
        _sum_0 += pow(cos(θ), 2) * pow((((p.at(i-1) - q.at(i-1))).dot(n.at(i-1)) + ((((p.at(i-1) + q.at(i-1))).cross(n.at(i-1)))).dot(a) + (n.at(i-1)).dot(t)), 2);
    }
    double b = _sum_0;

    return b;
}


void generateRandomData(double & θ,
    std::vector<Eigen::Matrix<double, 3, 1>> & p,
    std::vector<Eigen::Matrix<double, 3, 1>> & q,
    std::vector<Eigen::Matrix<double, 3, 1>> & n,
    Eigen::Matrix<double, 3, 1> & a,
    Eigen::Matrix<double, 3, 1> & t)
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
    a = Eigen::VectorXd::Random(3);
    t = Eigen::VectorXd::Random(3);
}


int main(int argc, char *argv[])
{
    double θ;
    std::vector<Eigen::Matrix<double, 3, 1>> p;
    std::vector<Eigen::Matrix<double, 3, 1>> q;
    std::vector<Eigen::Matrix<double, 3, 1>> n;
    Eigen::Matrix<double, 3, 1> a;
    Eigen::Matrix<double, 3, 1> t;
    generateRandomData(θ, p, q, n, a, t);
    double func_value = demo21(θ, p, q, n, a, t);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}