/*
`∂^2I_5/∂f^2` = 2[(A_1,1)I_3  (A_1,2)I_3  (A_1,3)I_3
                  (A_2,1)I_3  (A_2,2)I_3  (A_2,3)I_3
                  (A_3,1)I_3  (A_3,2)I_3  (A_3,3)I_3] 

where

A: ℝ^(3×3) 
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo17
 *
 * @param A  ℝ^(3×3)
 * @return _partial_differentialcircumflex_accent_2I_5_soliduspartial_differential_f_circumflex_accent_2
 */
Eigen::Matrix<double, 9, 9> demo17(const Eigen::Matrix<double, 3, 3> & A)
{
    Eigen::Matrix<double, 9, 9> __partial_differentialcircumflex_accent_2I_5_soliduspartial_differential_f_circumflex_accent_2_0;
    __partial_differentialcircumflex_accent_2I_5_soliduspartial_differential_f_circumflex_accent_2_0 << (A(1-1, 1-1)) * Eigen::MatrixXd::Identity(3, 3), (A(1-1, 2-1)) * Eigen::MatrixXd::Identity(3, 3), (A(1-1, 3-1)) * Eigen::MatrixXd::Identity(3, 3),
    (A(2-1, 1-1)) * Eigen::MatrixXd::Identity(3, 3), (A(2-1, 2-1)) * Eigen::MatrixXd::Identity(3, 3), (A(2-1, 3-1)) * Eigen::MatrixXd::Identity(3, 3),
    (A(3-1, 1-1)) * Eigen::MatrixXd::Identity(3, 3), (A(3-1, 2-1)) * Eigen::MatrixXd::Identity(3, 3), (A(3-1, 3-1)) * Eigen::MatrixXd::Identity(3, 3);
    Eigen::Matrix<double, 9, 9> _partial_differentialcircumflex_accent_2I_5_soliduspartial_differential_f_circumflex_accent_2 = 2 * __partial_differentialcircumflex_accent_2I_5_soliduspartial_differential_f_circumflex_accent_2_0;

    return _partial_differentialcircumflex_accent_2I_5_soliduspartial_differential_f_circumflex_accent_2;
}


void generateRandomData(Eigen::Matrix<double, 3, 3> & A)
{
    A = Eigen::MatrixXd::Random(3, 3);
}


int main(int argc, char *argv[])
{
    Eigen::Matrix<double, 3, 3> A;
    generateRandomData(A);
    Eigen::Matrix<double, 9, 9> func_value = demo17(A);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}