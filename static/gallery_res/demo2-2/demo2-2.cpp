/*
`n(v)` = (∑_(i for i ∈ `N₁(v)`) α_i n(T_i))/||∑_(i for i ∈ `N₁(v)`) α_i n(T_i)|| 

where
 
T_i: ℝ^(3×3)
α_i: ℝ
`N₁(v)`: {ℤ}
n: ℝ^(3×3) -> ℝ^3
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo2-2
 *
 * @param T  ℝ^(3×3)
 * @param α  ℝ
 * @param N₁_left_parenthesis_v_right_parenthesis  {ℤ}
 * @param n  ℝ^(3×3) -> ℝ^3
 * @return n_left_parenthesis_v_right_parenthesis
 */
Eigen::Matrix<double, 3, 1> demo2-2(
    const std::vector<Eigen::Matrix<double, 3, 3>> & T,
    const std::vector<double> & α,
    const std::set<std::tuple< int > > & N₁_left_parenthesis_v_right_parenthesis,
    const std::function<Eigen::Matrix<double, 3, 1>(Eigen::Matrix<double, 3, 3>)> & n)
{
    const long _dim_0 = α.size();
    assert( T.size() == _dim_0 );
    assert( α.size() == _dim_0 );

    Eigen::MatrixXd _sum_0 = Eigen::MatrixXd::Zero(3, 1);
    for(int i=1; i<=T.size(); i++){
        if(N₁_left_parenthesis_v_right_parenthesis.find(std::tuple< int >(i)) != N₁_left_parenthesis_v_right_parenthesis.end()){
            _sum_0 += α.at(i-1) * n(T.at(i-1));
        }
    }
    Eigen::MatrixXd _sum_1 = Eigen::MatrixXd::Zero(3, 1);
    for(int i=1; i<=T.size(); i++){
        if(N₁_left_parenthesis_v_right_parenthesis.find(std::tuple< int >(i)) != N₁_left_parenthesis_v_right_parenthesis.end()){
            _sum_1 += α.at(i-1) * n(T.at(i-1));
        }
    }
    Eigen::Matrix<double, 3, 1> n_left_parenthesis_v_right_parenthesis = (_sum_0) / double((_sum_1).lpNorm<2>());

    return n_left_parenthesis_v_right_parenthesis;
}


void generateRandomData(std::vector<Eigen::Matrix<double, 3, 3>> & T,
    std::vector<double> & α,
    std::set<std::tuple< int > > & N₁_left_parenthesis_v_right_parenthesis,
    std::function<Eigen::Matrix<double, 3, 1>(Eigen::Matrix<double, 3, 3>)> & n)
{
    const int _dim_0 = rand()%10;
    T.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        T[i] = Eigen::MatrixXd::Random(3, 3);
    }
    α.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        α[i] = rand() % 10;
    }
    const int N₁_left_parenthesis_v_right_parenthesis_0 = rand()%10;
    for(int i=0; i<N₁_left_parenthesis_v_right_parenthesis_0; i++){
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
    Eigen::Matrix<double, 3, 1> func_value = demo2-2(T, α, N₁_left_parenthesis_v_right_parenthesis, n);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}