/*
min_(C ∈ ℝ^3) ∑_i ‖x_i + (R_i - I₃)C‖²

where

x_i ∈ ℝ^3
R_i ∈ ℝ^(3×3)
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct hand_modeling_3ResultType {
    double ret;
    hand_modeling_3ResultType(const double & ret)
    : ret(ret)
    {}
};

hand_modeling_3ResultType hand_modeling_3(
    const std::vector<Eigen::Matrix<double, 3, 1>> & x,
    const std::vector<Eigen::Matrix<double, 3, 3>> & R)
{
    const long dim_0 = x.size();
    assert( R.size() == dim_0 );

    double ret = ;
    return hand_modeling_3ResultType(ret);
}


void generateRandomData(std::vector<Eigen::Matrix<double, 3, 1>> & x,
    std::vector<Eigen::Matrix<double, 3, 3>> & R)
{
    const int dim_0 = rand()%10;
    x.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        x[i] = Eigen::VectorXd::Random(3);
    }
    R.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        R[i] = Eigen::MatrixXd::Random(3, 3);
    }
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    std::vector<Eigen::Matrix<double, 3, 1>> x;
    std::vector<Eigen::Matrix<double, 3, 3>> R;
    generateRandomData(x, R);
    hand_modeling_3ResultType func_value = hand_modeling_3(x, R);
    std::cout<<"return value:\n"<<func_value.ret<<std::endl;
    return 0;
}