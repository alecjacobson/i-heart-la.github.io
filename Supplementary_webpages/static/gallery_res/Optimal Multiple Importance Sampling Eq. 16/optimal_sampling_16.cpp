/*
∑_i α_i + 1/M ∑_i ∑_(j for j ⩽ n_i) (f(X_i,j)/`p_c`(X_i,j) - (∑_k α_k p_k(X_i,j))/`p_c`(X_i,j))

where


α ∈ ℝ^N
p_j ∈ ℝ → ℝ 
X ∈ ℝ^(N×m)
M ∈ ℝ  
n_i ∈ ℝ
f: ℝ → ℝ 
`p_c`: ℝ → ℝ 
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct optimal_sampling_16ResultType {
    double ret;
    optimal_sampling_16ResultType(const double & ret)
    : ret(ret)
    {}
};

/**
 * optimal_sampling_16
 *
 * @param f  ℝ → ℝ
 * @param p_c  ℝ → ℝ
 * @return ret
 */
optimal_sampling_16ResultType optimal_sampling_16(
    const Eigen::VectorXd & α,
    const std::vector<std::function<double(double)>> & p,
    const Eigen::MatrixXd & X,
    const double & M,
    const std::vector<double> & n,
    const std::function<double(double)> & f,
    const std::function<double(double)> & p_c)
{
    const long dim_0 = n.size();
    const long N = α.size();
    const long m = X.cols();
    const long dim_1 = p.size();
    assert( X.rows() == N );
    assert( N == dim_1 );

    double sum_0 = 0;
    for(int i=1; i<=α.size(); i++){
        sum_0 += α[i-1];
    }
    double sum_1 = 0;
    for(int i=1; i<=X.rows(); i++){
        double sum_2 = 0;
        for(int j=1; j<=X.cols(); j++){
            double sum_3 = 0;
            for(int k=1; k<=α.size(); k++){
                sum_3 += α[k-1] * p.at(k-1)(X(i-1, j-1));
            }
            if(j <= n.at(i-1)){
                sum_2 += (f(X(i-1, j-1)) / double(p_c(X(i-1, j-1))) - (sum_3) / double(p_c(X(i-1, j-1))));
            }
        }
        sum_1 += sum_2;
    }
    double ret = sum_0 + 1 / double(M) * sum_1;
    return optimal_sampling_16ResultType(ret);
}


void generateRandomData(Eigen::VectorXd & α,
    std::vector<std::function<double(double)>> & p,
    Eigen::MatrixXd & X,
    double & M,
    std::vector<double> & n,
    std::function<double(double)> & f,
    std::function<double(double)> & p_c)
{
    M = rand() % 10;
    const int dim_0 = rand()%10;
    const int N = rand()%10;
    const int dim_1 = N;
    const int m = rand()%10;
    α = Eigen::VectorXd::Random(N);
    p.resize(dim_1);
    for(int i=0; i<dim_1; i++){
        p[i] = [](double)->double{
            return rand() % 10;
        };
    }
    X = Eigen::MatrixXd::Random(N, m);
    n.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        n[i] = rand() % 10;
    }
    f = [](double)->double{
        return rand() % 10;
    };
    p_c = [](double)->double{
        return rand() % 10;
    };
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    Eigen::VectorXd α;
    std::vector<std::function<double(double)>> p;
    Eigen::MatrixXd X;
    double M;
    std::vector<double> n;
    std::function<double(double)> f;
    std::function<double(double)> p_c;
    generateRandomData(α, p, X, M, n, f, p_c);
    optimal_sampling_16ResultType func_value = optimal_sampling_16(α, p, X, M, n, f, p_c);
    std::cout<<"return value:\n"<<func_value.ret<<std::endl;
    return 0;
}