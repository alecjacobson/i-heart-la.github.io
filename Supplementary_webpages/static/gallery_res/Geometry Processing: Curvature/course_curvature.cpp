/*
`H(p)` = 1/(2π)int_[0, 2π] `kₙ`(φ, p) ∂φ

where 

p ∈ ℝ^3 : point on the surface
`kₙ`: ℝ,ℝ^3 → ℝ : normal curvature
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct course_curvatureResultType {
    double H_left_parenthesis_p_right_parenthesis;
    course_curvatureResultType(const double & H_left_parenthesis_p_right_parenthesis)
    : H_left_parenthesis_p_right_parenthesis(H_left_parenthesis_p_right_parenthesis)
    {}
};

/**
 * course_curvature
 *
 * @param p  point on the surface
 * @param kₙ  ℝ,ℝ^3 → ℝ : normal curvature
 * @return H_left_parenthesis_p_right_parenthesis
 */
course_curvatureResultType course_curvature(
    const Eigen::Matrix<double, 3, 1> & p,
    const std::function<double(double, Eigen::Matrix<double, 3, 1>)> & kₙ)
{
    double H_left_parenthesis_p_right_parenthesis = 1 / double((2 * M_PI)) * ;

    return course_curvatureResultType(H_left_parenthesis_p_right_parenthesis);
}


void generateRandomData(Eigen::Matrix<double, 3, 1> & p,
    std::function<double(double, Eigen::Matrix<double, 3, 1>)> & kₙ)
{
    p = Eigen::VectorXd::Random(3);
    kₙ = [](double, Eigen::Matrix<double, 3, 1>)->double{
        return rand() % 10;
    };
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    Eigen::Matrix<double, 3, 1> p;
    std::function<double(double, Eigen::Matrix<double, 3, 1>)> kₙ;
    generateRandomData(p, kₙ);
    course_curvatureResultType func_value = course_curvature(p, kₙ);
    std::cout<<"return value:\n"<<func_value.H_left_parenthesis_p_right_parenthesis<<std::endl;
    return 0;
}