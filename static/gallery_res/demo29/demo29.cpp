/*
sum_i α_i + 1/M sum_i sum_j (f(X_i,j)/`p_c`(X_i,j) - (sum_k α_k p_k X_i,j)/`p_c`(X_i,j))

where


α: ℝ^m
p: ℝ^m
X: ℝ^(m×n)
M: ℝ  
f: ℝ -> ℝ 
`p_c`: ℝ -> ℝ 
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo29
 *
 * @param α  ℝ^m
 * @param p  ℝ^m
 * @param X  ℝ^(m×n)
 * @param M  ℝ
 * @param f  ℝ -> ℝ
 * @param p_c  ℝ -> ℝ
 * @return ret
 */
double demo29(
    const Eigen::VectorXd & α,
    const Eigen::VectorXd & p,
    const Eigen::MatrixXd & X,
    const double & M,
    const std::function<double(double)> & f,
    const std::function<double(double)> & p_c)
{
    const long m = X.rows();
    const long n = X.cols();
    assert( α.size() == m );
    assert( p.size() == m );
    assert( X.rows() == m );
    assert( X.cols() == n );

    double _sum_0 = 0;
    for(int i=1; i<=α.size(); i++){
        _sum_0 += α(i-1);
    }
    double _sum_1 = 0;
    for(int i=1; i<=X.rows(); i++){
        double _sum_2 = 0;
        for(int j=1; j<=X.cols(); j++){
            double _sum_3 = 0;
            for(int k=1; k<=p.size(); k++){
                _sum_3 += α(k-1) * p(k-1) * X(i-1, j-1);
            }
            _sum_2 += (f(X(i-1, j-1)) / p_c(X(i-1, j-1)) - (_sum_3) / p_c(X(i-1, j-1)));
        }
        _sum_1 += _sum_2;
    }
    double ret = _sum_0 + 1 / M * _sum_1;
    return ret;
}


void generateRandomData(Eigen::VectorXd & α,
    Eigen::VectorXd & p,
    Eigen::MatrixXd & X,
    double & M,
    std::function<double(double)> & f,
    std::function<double(double)> & p_c)
{
    M = rand() % 10;
    const int m = rand()%10;
    const int n = rand()%10;
    α = Eigen::VectorXd::Random(m);
    p = Eigen::VectorXd::Random(m);
    X = Eigen::MatrixXd::Random(m, n);
    f = [](double)->double{
        return rand() % 10;
    };
    p_c = [](double)->double{
        return rand() % 10;
    };
}


int main(int argc, char *argv[])
{
    Eigen::VectorXd α;
    Eigen::VectorXd p;
    Eigen::MatrixXd X;
    double M;
    std::function<double(double)> f;
    std::function<double(double)> p_c;
    generateRandomData(α, p, X, M, f, p_c);
    double func_value = demo29(α, p, X, M, f, p_c);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}