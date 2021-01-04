/*
`E(θ,σ²)` = -sum_n log(sum_m 1/||y|| 1/(2σ²)exp(-||x_n - y_m||²/(2σ²)) + 1/||x||)

where


x_i: ℝ^2
y_j: ℝ^2 
σ: ℝ
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo43
 *
 * @param x  ℝ^2
 * @param y  ℝ^2
 * @param σ  ℝ
 * @return E_left_parenthesis_θ_comma_σ²_right_parenthesis
 */
double demo43(
    const std::vector<Eigen::Matrix<double, 2, 1>> & x,
    const std::vector<Eigen::Matrix<double, 2, 1>> & y,
    const double & σ)
{
    const long _dim_0 = x.size();
    const long _dim_1 = y.size();
    double _sum_0 = 0;
    for(int n=1; n<=x.size(); n++){
        double _sum_1 = 0;
        for(int m=1; m<=y.size(); m++){
            _sum_1 += 1 /  * 1 / (2 * pow(σ, 2)) * exp(-pow((x.at(n-1) - y.at(m-1)).lpNorm<2>(), 2) / (2 * pow(σ, 2)));
        }
        _sum_0 += log10(_sum_1 + 1 / );
    }
    double E_left_parenthesis_θ_comma_σ²_right_parenthesis = -_sum_0;

    return E_left_parenthesis_θ_comma_σ²_right_parenthesis;
}


void generateRandomData(std::vector<Eigen::Matrix<double, 2, 1>> & x,
    std::vector<Eigen::Matrix<double, 2, 1>> & y,
    double & σ)
{
    σ = rand() % 10;
    const int _dim_0 = rand()%10;
    const int _dim_1 = rand()%10;
    x.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        x[i] = Eigen::VectorXd::Random(2);
    }
    y.resize(_dim_1);
    for(int i=0; i<_dim_1; i++){
        y[i] = Eigen::VectorXd::Random(2);
    }
}


int main(int argc, char *argv[])
{
    std::vector<Eigen::Matrix<double, 2, 1>> x;
    std::vector<Eigen::Matrix<double, 2, 1>> y;
    double σ;
    generateRandomData(x, y, σ);
    double func_value = demo43(x, y, σ);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}