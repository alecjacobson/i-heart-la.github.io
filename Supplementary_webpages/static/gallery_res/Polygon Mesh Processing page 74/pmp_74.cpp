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

struct pmp_74ResultType {
    double E_LSCM;
    pmp_74ResultType(const double & E_LSCM)
    : E_LSCM(E_LSCM)
    {}
};

pmp_74ResultType pmp_74(
    const std::vector<Eigen::Matrix<double, 3, 1>> & v,
    const std::vector<Eigen::Matrix<double, 3, 1>> & u,
    const std::vector<Eigen::Matrix<double, 2, 3>> & M,
    const std::vector<double> & A)
{
    const long dim_0 = A.size();
    assert( v.size() == dim_0 );
    assert( u.size() == dim_0 );
    assert( M.size() == dim_0 );

    double sum_0 = 0;
    for(int T=1; T<=u.size(); T++){
        Eigen::Matrix<double, 2, 2> E_LSCM_0;
        E_LSCM_0 << 0, -1,
        1, 0;
        sum_0 += A.at(T-1) * pow((M.at(T-1) * v.at(T-1) - E_LSCM_0 * M.at(T-1) * u.at(T-1)).lpNorm<2>(), 2);
    }
    double E_LSCM = sum_0;

    return pmp_74ResultType(E_LSCM);
}


void generateRandomData(std::vector<Eigen::Matrix<double, 3, 1>> & v,
    std::vector<Eigen::Matrix<double, 3, 1>> & u,
    std::vector<Eigen::Matrix<double, 2, 3>> & M,
    std::vector<double> & A)
{
    const int dim_0 = rand()%10;
    v.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        v[i] = Eigen::VectorXd::Random(3);
    }
    u.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        u[i] = Eigen::VectorXd::Random(3);
    }
    M.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        M[i] = Eigen::MatrixXd::Random(2, 3);
    }
    A.resize(dim_0);
    for(int i=0; i<dim_0; i++){
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
    pmp_74ResultType func_value = pmp_74(v, u, M, A);
    std::cout<<"return value:\n"<<func_value.E_LSCM<<std::endl;
    return 0;
}