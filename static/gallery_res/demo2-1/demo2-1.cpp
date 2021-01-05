/*
`xᵢ` = T_*,1
`xⱼ` = T_*,2
`xₖ` = T_*,3
`n(T)` = (`xⱼ`-`xᵢ`)×(`xₖ`-`xᵢ`)/||(`xⱼ`-`xᵢ`)×(`xₖ`-`xᵢ`)||

where
 
T: ℝ^(3×3)
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo2-1
 *
 * @param T  ℝ^(3×3)
 * @return n_left_parenthesis_T_right_parenthesis
 */
Eigen::Matrix<double, 3, 1> demo2-1(const Eigen::Matrix<double, 3, 3> & T)
{
    Eigen::Matrix<double, 3, 1> xᵢ = T.col(1-1);

    Eigen::Matrix<double, 3, 1> xⱼ = T.col(2-1);

    Eigen::Matrix<double, 3, 1> xₖ = T.col(3-1);

    Eigen::Matrix<double, 3, 1> n_left_parenthesis_T_right_parenthesis = ((xⱼ - xᵢ)).cross((xₖ - xᵢ)) / (((xⱼ - xᵢ)).cross((xₖ - xᵢ))).lpNorm<2>();

    return n_left_parenthesis_T_right_parenthesis;
}


void generateRandomData(Eigen::Matrix<double, 3, 3> & T)
{
    T = Eigen::MatrixXd::Random(3, 3);
}


int main(int argc, char *argv[])
{
    Eigen::Matrix<double, 3, 3> T;
    generateRandomData(T);
    Eigen::Matrix<double, 3, 1> func_value = demo2-1(T);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}