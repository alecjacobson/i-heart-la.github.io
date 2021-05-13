/*
`n(v)` = (∑_(i for i ∈ `N₁(v)`) α_i n(T_i))/‖∑_(i for i ∈ `N₁(v)`) α_i n(T_i)‖

where
 
T_i ∈ ℝ^(3×3)
α_i ∈ ℝ
`N₁(v)` ∈ {ℤ}
n: ℝ^(3×3) → ℝ^3
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct pmp_42ResultType {
    Eigen::Matrix<double, 3, 1> n_left_parenthesis_v_right_parenthesis;
    pmp_42ResultType(const Eigen::Matrix<double, 3, 1> & n_left_parenthesis_v_right_parenthesis)
    : n_left_parenthesis_v_right_parenthesis(n_left_parenthesis_v_right_parenthesis)
    {}
};

/**
 * pmp_42
 *
 * @param n  ℝ^(3×3) → ℝ^3
 * @return n_left_parenthesis_v_right_parenthesis
 */
pmp_42ResultType pmp_42(
    const std::vector<Eigen::Matrix<double, 3, 3>> & T,
    const std::vector<double> & α,
    const std::set<std::tuple< int > > & N₁_left_parenthesis_v_right_parenthesis,
    const std::function<Eigen::Matrix<double, 3, 1>(Eigen::Matrix<double, 3, 3>)> & n)
{
    const long dim_0 = α.size();
    assert( T.size() == dim_0 );

    Eigen::MatrixXd sum_0 = Eigen::MatrixXd::Zero(3, 1);
    for(int i=1; i<=α.size(); i++){
        if(N₁_left_parenthesis_v_right_parenthesis.find(std::tuple< int >(i)) != N₁_left_parenthesis_v_right_parenthesis.end()){
            sum_0 += α.at(i-1) * n(T.at(i-1));
        }
    }
    Eigen::MatrixXd sum_1 = Eigen::MatrixXd::Zero(3, 1);
    for(int i=1; i<=α.size(); i++){
        if(N₁_left_parenthesis_v_right_parenthesis.find(std::tuple< int >(i)) != N₁_left_parenthesis_v_right_parenthesis.end()){
            sum_1 += α.at(i-1) * n(T.at(i-1));
        }
    }
    Eigen::Matrix<double, 3, 1> n_left_parenthesis_v_right_parenthesis = (sum_0) / double((sum_1).lpNorm<2>());

    return pmp_42ResultType(n_left_parenthesis_v_right_parenthesis);
}


void generateRandomData(std::vector<Eigen::Matrix<double, 3, 3>> & T,
    std::vector<double> & α,
    std::set<std::tuple< int > > & N₁_left_parenthesis_v_right_parenthesis,
    std::function<Eigen::Matrix<double, 3, 1>(Eigen::Matrix<double, 3, 3>)> & n)
{
    const int dim_0 = rand()%10;
    T.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        T[i] = Eigen::MatrixXd::Random(3, 3);
    }
    α.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        α[i] = rand() % 10;
    }
    const int dim_2 = rand()%10;
    for(int i=0; i<dim_2; i++){
        N₁_left_parenthesis_v_right_parenthesis.insert(std::make_tuple(rand()%10));
    }
    n = [](Eigen::Matrix<double, 3, 3>)->Eigen::Matrix<double, 3, 1>{
        return Eigen::VectorXd::Random(3);
    };
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    std::vector<Eigen::Matrix<double, 3, 3>> T;
    std::vector<double> α;
    std::set<std::tuple< int > > N₁_left_parenthesis_v_right_parenthesis;
    std::function<Eigen::Matrix<double, 3, 1>(Eigen::Matrix<double, 3, 3>)> n;
    generateRandomData(T, α, N₁_left_parenthesis_v_right_parenthesis, n);
    pmp_42ResultType func_value = pmp_42(T, α, N₁_left_parenthesis_v_right_parenthesis, n);
    std::cout<<"return value:\n"<<func_value.n_left_parenthesis_v_right_parenthesis<<std::endl;
    return 0;
}