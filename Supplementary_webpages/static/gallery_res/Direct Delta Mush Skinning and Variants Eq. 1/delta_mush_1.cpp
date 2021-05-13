/*
v_i = ∑_j w_i,j M_j u_i

where

w ∈ ℝ^(n×m)
M_j ∈ ℝ^(4×4)
u_i ∈ ℝ^4
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct delta_mush_1ResultType {
    std::vector<Eigen::Matrix<double, 4, 1>> v;
    delta_mush_1ResultType(const std::vector<Eigen::Matrix<double, 4, 1>> & v)
    : v(v)
    {}
};

delta_mush_1ResultType delta_mush_1(
    const Eigen::MatrixXd & w,
    const std::vector<Eigen::Matrix<double, 4, 4>> & M,
    const std::vector<Eigen::Matrix<double, 4, 1>> & u)
{
    const long n = w.rows();
    const long m = w.cols();
    const long dim_0 = M.size();
    const long dim_1 = u.size();
    assert( w.rows() == n );
    assert( dim_0 == m );
    assert( dim_1 == n );

    std::vector<Eigen::Matrix<double, 4, 1>> v(dim_1);
    for( int i=1; i<=dim_1; i++){
        Eigen::MatrixXd sum_0 = Eigen::MatrixXd::Zero(4, 1);
        for(int j=1; j<=w.cols(); j++){
            sum_0 += w(i-1, j-1) * M.at(j-1) * u.at(i-1);
        }
        v.at(i-1) = sum_0;
    }

    return delta_mush_1ResultType(v);
}


void generateRandomData(Eigen::MatrixXd & w,
    std::vector<Eigen::Matrix<double, 4, 4>> & M,
    std::vector<Eigen::Matrix<double, 4, 1>> & u)
{
    const int n = rand()%10;
    const int dim_1 = n;
    const int m = rand()%10;
    const int dim_0 = m;
    w = Eigen::MatrixXd::Random(n, m);
    M.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        M[i] = Eigen::MatrixXd::Random(4, 4);
    }
    u.resize(dim_1);
    for(int i=0; i<dim_1; i++){
        u[i] = Eigen::VectorXd::Random(4);
    }
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    Eigen::MatrixXd w;
    std::vector<Eigen::Matrix<double, 4, 4>> M;
    std::vector<Eigen::Matrix<double, 4, 1>> u;
    generateRandomData(w, M, u);
    delta_mush_1ResultType func_value = delta_mush_1(w, M, u);
    std::cout<<"vector return value:"<<std::endl;
    for(int i=0; i<func_value.v.size(); i++){
        std::cout<<"i:"<<i<<", value:\n"<<func_value.v.at(i)<<std::endl;
    }
    return 0;
}