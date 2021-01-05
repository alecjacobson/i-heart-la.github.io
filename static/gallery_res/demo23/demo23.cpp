/*
r̄ = v̄×ō
s̄ = ō×ū
n̄ = ū×v̄

`k_r` = r̄⋅(`C̄ₐ`-V̄)
`k_s` = s̄⋅(`C̄ₐ`-V̄)
`k_n` = n̄⋅(`C̄ₐ`-V̄)

`x(θ,v)` =  (r̄⋅`D_A`(θ, v)+`k_r`δ(θ, v))/(n̄⋅`D_A`(θ, v)+`k_n`δ(θ, v))
`y(θ,v)` =  (s̄⋅`D_A`(θ, v)+`k_s`δ(θ, v))/(n̄⋅`D_A`(θ, v)+`k_n`δ(θ, v))

where

v̄: ℝ^3
ō: ℝ^3
ū: ℝ^3
V̄: ℝ^3
`C̄ₐ`: ℝ^3
θ: ℝ 
v: ℝ 
`D_A`: ℝ,ℝ->ℝ^3
δ: ℝ,ℝ->ℝ 
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo23
 *
 * @param v̄  ℝ^3
 * @param ō  ℝ^3
 * @param ū  ℝ^3
 * @param V̄  ℝ^3
 * @param C_combining_macron_ₐ  ℝ^3
 * @param θ  ℝ
 * @param v  ℝ
 * @param D_A  ℝ,ℝ->ℝ^3
 * @param δ  ℝ,ℝ->ℝ
 * @return y_left_parenthesis_θ_comma_v_right_parenthesis
 */
double demo23(
    const Eigen::Matrix<double, 3, 1> & v̄,
    const Eigen::Matrix<double, 3, 1> & ō,
    const Eigen::Matrix<double, 3, 1> & ū,
    const Eigen::Matrix<double, 3, 1> & V̄,
    const Eigen::Matrix<double, 3, 1> & C_combining_macron_ₐ,
    const double & θ,
    const double & v,
    const std::function<Eigen::Matrix<double, 3, 1>(double, double)> & D_A,
    const std::function<double(double, double)> & δ)
{
    Eigen::Matrix<double, 3, 1> r̄ = (v̄).cross(ō);

    Eigen::Matrix<double, 3, 1> s̄ = (ō).cross(ū);

    Eigen::Matrix<double, 3, 1> n̄ = (ū).cross(v̄);

    double k_r = (r̄).dot((C_combining_macron_ₐ - V̄));

    double k_s = (s̄).dot((C_combining_macron_ₐ - V̄));

    double k_n = (n̄).dot((C_combining_macron_ₐ - V̄));

    double x_left_parenthesis_θ_comma_v_right_parenthesis = ((r̄).dot(D_A(θ, v)) + k_r * δ(θ, v)) / ((n̄).dot(D_A(θ, v)) + k_n * δ(θ, v));

    double y_left_parenthesis_θ_comma_v_right_parenthesis = ((s̄).dot(D_A(θ, v)) + k_s * δ(θ, v)) / ((n̄).dot(D_A(θ, v)) + k_n * δ(θ, v));

    return y_left_parenthesis_θ_comma_v_right_parenthesis;
}


void generateRandomData(Eigen::Matrix<double, 3, 1> & v̄,
    Eigen::Matrix<double, 3, 1> & ō,
    Eigen::Matrix<double, 3, 1> & ū,
    Eigen::Matrix<double, 3, 1> & V̄,
    Eigen::Matrix<double, 3, 1> & C_combining_macron_ₐ,
    double & θ,
    double & v,
    std::function<Eigen::Matrix<double, 3, 1>(double, double)> & D_A,
    std::function<double(double, double)> & δ)
{
    θ = rand() % 10;
    v = rand() % 10;
    v̄ = Eigen::VectorXd::Random(3);
    ō = Eigen::VectorXd::Random(3);
    ū = Eigen::VectorXd::Random(3);
    V̄ = Eigen::VectorXd::Random(3);
    C_combining_macron_ₐ = Eigen::VectorXd::Random(3);
    D_A = [](double, double)->Eigen::Matrix<double, 3, 1>{
        return Eigen::VectorXd::Random(3);
    };
    δ = [](double, double)->double{
        return rand() % 10;
    };
}


int main(int argc, char *argv[])
{
    Eigen::Matrix<double, 3, 1> v̄;
    Eigen::Matrix<double, 3, 1> ō;
    Eigen::Matrix<double, 3, 1> ū;
    Eigen::Matrix<double, 3, 1> V̄;
    Eigen::Matrix<double, 3, 1> C_combining_macron_ₐ;
    double θ;
    double v;
    std::function<Eigen::Matrix<double, 3, 1>(double, double)> D_A;
    std::function<double(double, double)> δ;
    generateRandomData(v̄, ō, ū, V̄, C_combining_macron_ₐ, θ, v, D_A, δ);
    double func_value = demo23(v̄, ō, ū, V̄, C_combining_macron_ₐ, θ, v, D_A, δ);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}