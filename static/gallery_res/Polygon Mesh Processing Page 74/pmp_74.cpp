/*
`E_LSCM` = ∑_T A_T‖M_T v_T - [0 -1
                              1  0] M_T u_T‖²
where
 
v_i ∈ ℝ^3
u_i ∈ ℝ^3
M_i ∈ ℝ^(2×3)
A_i ∈ ℝ
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

double pmp_74(
    const std::vector<Eigen::Matrix<double, 3, 1>> & v,
    const std::vector<Eigen::Matrix<double, 3, 1>> & u,
    const std::vector<Eigen::Matrix<double, 2, 3>> & M,
    const std::vector<double> & A)
{
    const long _dim_0 = A.size();
    assert( M.size() == _dim_0 );
    assert( A.size() == _dim_0 );

    double _sum_0 = 0;
    for(int T=1; T<=u.size(); T++){
        Eigen::Matrix<double, 2, 2> _E_LSCM_0;
        _E_LSCM_0 << 0, -1,
        1, 0;
        _sum_0 += A.at(T-1) * pow((M.at(T-1) * v.at(T-1) - _E_LSCM_0 * M.at(T-1) * u.at(T-1)).lpNorm<2>(), 2);
    }
    double E_LSCM = _sum_0;

    return E_LSCM;
}


void generateRandomData(std::vector<Eigen::Matrix<double, 3, 1>> & v,
    std::vector<Eigen::Matrix<double, 3, 1>> & u,
    std::vector<Eigen::Matrix<double, 2, 3>> & M,
    std::vector<double> & A)
{
    const int _dim_0 = rand()%10;
    v.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        v[i] = Eigen::VectorXd::Random(3);
    }
    u.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        u[i] = Eigen::VectorXd::Random(3);
    }
    M.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        M[i] = Eigen::MatrixXd::Random(2, 3);
    }
    A.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        A[i] = rand() % 10;
    }
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    std::vector<Eigen::Matrix<double, 3, 1>> v;
    std::vector<Eigen::Matrix<double, 3, 1>> u;
    std::vector<Eigen::Matrix<double, 2, 3>> M;
    std::vector<double> A;
    generateRandomData(v, u, M, A);
    double func_value = pmp_74(v, u, M, A);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}