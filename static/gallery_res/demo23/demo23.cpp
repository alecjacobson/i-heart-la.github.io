/*
`x(θ,v)` =  (`r bar`⋅`D_A`(`θ`, v)+`k_r``δ`(`θ`, v))/(`n bar`⋅`D_A`(`θ`, v)+`k_n``δ`(`θ`, v))
`y(θ,v)` =  (`s bar`⋅`D_A`(`θ`, v)+`k_s``δ`(`θ`, v))/(`n bar`⋅`D_A`(`θ`, v)+`k_n``δ`(`θ`, v))

where

`r bar`: ℝ^3
`s bar`: ℝ^3
`n bar`: ℝ^3
`θ`: ℝ 
v: ℝ 
`k_r`: ℝ 
`k_s`: ℝ 
`k_n`: ℝ 
`D_A`: ℝ,ℝ->ℝ^3
`δ`: ℝ,ℝ->ℝ  
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo23
 *
 * @param r_bar  ℝ^3
 * @param s_bar  ℝ^3
 * @param n_bar  ℝ^3
 * @param θ  ℝ
 * @param v  ℝ
 * @param k_r  ℝ
 * @param k_s  ℝ
 * @param k_n  ℝ
 * @param D_A  ℝ,ℝ->ℝ^3
 * @param δ  ℝ,ℝ->ℝ
 * @return y_left_parenthesis_θ_comma_v_right_parenthesis
 */
double demo23(
    const Eigen::Matrix<double, 3, 1> & r_bar,
    const Eigen::Matrix<double, 3, 1> & s_bar,
    const Eigen::Matrix<double, 3, 1> & n_bar,
    const double & θ,
    const double & v,
    const double & k_r,
    const double & k_s,
    const double & k_n,
    const std::function<Eigen::Matrix<double, 3, 1>(double, double)> & D_A,
    const std::function<double(double, double)> & δ)
{
    double x_left_parenthesis_θ_comma_v_right_parenthesis = ((r_bar).dot(D_A(θ, v)) + k_r * δ(θ, v)) / ((n_bar).dot(D_A(θ, v)) + k_n * δ(θ, v));

    double y_left_parenthesis_θ_comma_v_right_parenthesis = ((s_bar).dot(D_A(θ, v)) + k_s * δ(θ, v)) / ((n_bar).dot(D_A(θ, v)) + k_n * δ(θ, v));

    return y_left_parenthesis_θ_comma_v_right_parenthesis;
}


void generateRandomData(Eigen::Matrix<double, 3, 1> & r_bar,
    Eigen::Matrix<double, 3, 1> & s_bar,
    Eigen::Matrix<double, 3, 1> & n_bar,
    double & θ,
    double & v,
    double & k_r,
    double & k_s,
    double & k_n,
    std::function<Eigen::Matrix<double, 3, 1>(double, double)> & D_A,
    std::function<double(double, double)> & δ)
{
    θ = rand() % 10;
    v = rand() % 10;
    k_r = rand() % 10;
    k_s = rand() % 10;
    k_n = rand() % 10;
    r_bar = Eigen::VectorXd::Random(3);
    s_bar = Eigen::VectorXd::Random(3);
    n_bar = Eigen::VectorXd::Random(3);
    D_A = [](double, double)->Eigen::Matrix<double, 3, 1>{
        return Eigen::VectorXd::Random(3);
    };
    δ = [](double, double)->double{
        return rand() % 10;
    };
}


int main(int argc, char *argv[])
{
    Eigen::Matrix<double, 3, 1> r_bar;
    Eigen::Matrix<double, 3, 1> s_bar;
    Eigen::Matrix<double, 3, 1> n_bar;
    double θ;
    double v;
    double k_r;
    double k_s;
    double k_n;
    std::function<Eigen::Matrix<double, 3, 1>(double, double)> D_A;
    std::function<double(double, double)> δ;
    generateRandomData(r_bar, s_bar, n_bar, θ, v, k_r, k_s, k_n, D_A, δ);
    double func_value = demo23(r_bar, s_bar, n_bar, θ, v, k_r, k_s, k_n, D_A, δ);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}