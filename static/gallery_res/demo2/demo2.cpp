/*
n_i = (T_i,*,2-T_i,*,1)×(T_i,*,3-T_i,*,1)/||(T_i,*,2-T_i,*,1)×(T_i,*,3-T_i,*,1)||

`n(v)` = (∑_(i for i ∈ `N1(v)`) α_i n_i)/||∑_(i for i ∈ `N1(v)`) α_i n_i||
where
 
T_i: ℝ^(3×3)
α_i: ℝ
`N1(v)`: {ℤ}
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo2
 *
 * @param T  ℝ^(3×3)
 * @param α  ℝ
 * @param N1_left_parenthesis_v_right_parenthesis  {ℤ}
 * @return n_left_parenthesis_v_right_parenthesis
 */
Eigen::Matrix<double, 3, 1> demo2(
    const std::vector<Eigen::Matrix<double, 3, 3>> & T,
    const std::vector<double> & α,
    const std::set<std::tuple< int > > & N1_left_parenthesis_v_right_parenthesis)
{
    const long _dim_0 = α.size();
    assert( T.size() == _dim_0 );
    assert( α.size() == _dim_0 );

    std::vector<Eigen::Matrix<double, 3, 1>> n(_dim_0);
    for( int i=1; i<=_dim_0; i++){
        n.at(i-1) = ((T.at(i-1).col(2-1) - T.at(i-1).col(1-1))).cross((T.at(i-1).col(3-1) - T.at(i-1).col(1-1))) / (((T.at(i-1).col(2-1) - T.at(i-1).col(1-1))).cross((T.at(i-1).col(3-1) - T.at(i-1).col(1-1)))).lpNorm<2>();
    }


    Eigen::MatrixXd _sum_0 = Eigen::MatrixXd::Zero(3, 1);
    for(int i=1; i<=n.size(); i++){
        if(N1_left_parenthesis_v_right_parenthesis.find(std::tuple< int >(i)) != N1_left_parenthesis_v_right_parenthesis.end()){
            _sum_0 += α.at(i-1) * n.at(i-1);
        }
    }
    Eigen::MatrixXd _sum_1 = Eigen::MatrixXd::Zero(3, 1);
    for(int i=1; i<=n.size(); i++){
        if(N1_left_parenthesis_v_right_parenthesis.find(std::tuple< int >(i)) != N1_left_parenthesis_v_right_parenthesis.end()){
            _sum_1 += α.at(i-1) * n.at(i-1);
        }
    }
    Eigen::Matrix<double, 3, 1> n_left_parenthesis_v_right_parenthesis = (_sum_0) / (_sum_1).lpNorm<2>();

    return n_left_parenthesis_v_right_parenthesis;
}


void generateRandomData(std::vector<Eigen::Matrix<double, 3, 3>> & T,
    std::vector<double> & α,
    std::set<std::tuple< int > > & N1_left_parenthesis_v_right_parenthesis)
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
    const int N1_left_parenthesis_v_right_parenthesis_0 = rand()%10;
    for(int i=0; i<N1_left_parenthesis_v_right_parenthesis_0; i++){
        N1_left_parenthesis_v_right_parenthesis.insert(std::make_tuple(rand()%10));
    }
}


int main(int argc, char *argv[])
{
    std::vector<Eigen::Matrix<double, 3, 3>> T;
    std::vector<double> α;
    std::set<std::tuple< int > > N1_left_parenthesis_v_right_parenthesis;
    generateRandomData(T, α, N1_left_parenthesis_v_right_parenthesis);
    Eigen::Matrix<double, 3, 1> func_value = demo2(T, α, N1_left_parenthesis_v_right_parenthesis);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}