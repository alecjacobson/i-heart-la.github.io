/*
`G_σ(s_i^k)` = ∑_j l_j exp(-dist(`bᵢ`, b_j)/(2σ²))(s_j)^k

where

l_j: ℝ : the length of bj
dist: ℝ^n, ℝ^n -> ℝ : measures the geodesic distance between the centers of bi and bj along the boundary
σ: ℝ
`bᵢ`: ℝ^n
b_j: ℝ^n
s_j: ℝ : unit direction vector of bi
k: ℝ : iteration number
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo28
 *
 * @param l  ℝ : the length of bj
 * @param dist  ℝ^n, ℝ^n -> ℝ : measures the geodesic distance between the centers of bi and bj along the boundary
 * @param σ  ℝ
 * @param bᵢ  ℝ^n
 * @param b  ℝ^n
 * @param s  ℝ : unit direction vector of bi
 * @param k  ℝ : iteration number
 * @return G_σ_left_parenthesis_s_i_circumflex_accent_k_right_parenthesis
 */
double demo28(
    const std::vector<double> & l,
    const std::function<double(Eigen::VectorXd, Eigen::VectorXd)> & dist,
    const double & σ,
    const Eigen::VectorXd & bᵢ,
    const std::vector<Eigen::VectorXd> & b,
    const std::vector<double> & s,
    const double & k)
{
    const long _dim_0 = l.size();
    const long n = b[0].rows();
    assert( l.size() == _dim_0 );
    assert( bᵢ.size() == n );
    assert( b.size() == _dim_0 );
    for( const auto& el : b ) {
        assert( el.size() == n );
    }
    assert( s.size() == _dim_0 );

    double _sum_0 = 0;
    for(int j=1; j<=s.size(); j++){
        _sum_0 += l.at(j-1) * exp(-dist(bᵢ, b.at(j-1)) / double((2 * pow(σ, 2)))) * pow((s.at(j-1)), k);
    }
    double G_σ_left_parenthesis_s_i_circumflex_accent_k_right_parenthesis = _sum_0;

    return G_σ_left_parenthesis_s_i_circumflex_accent_k_right_parenthesis;
}


void generateRandomData(std::vector<double> & l,
    std::function<double(Eigen::VectorXd, Eigen::VectorXd)> & dist,
    double & σ,
    Eigen::VectorXd & bᵢ,
    std::vector<Eigen::VectorXd> & b,
    std::vector<double> & s,
    double & k)
{
    σ = rand() % 10;
    k = rand() % 10;
    const int _dim_0 = rand()%10;
    const int n = rand()%10;
    l.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        l[i] = rand() % 10;
    }
    dist = [](Eigen::VectorXd, Eigen::VectorXd)->double{
        return rand() % 10;
    };
    bᵢ = Eigen::VectorXd::Random(n);
    b.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        b[i] = Eigen::VectorXd::Random(n);
    }
    s.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        s[i] = rand() % 10;
    }
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    std::vector<double> l;
    std::function<double(Eigen::VectorXd, Eigen::VectorXd)> dist;
    double σ;
    Eigen::VectorXd bᵢ;
    std::vector<Eigen::VectorXd> b;
    std::vector<double> s;
    double k;
    generateRandomData(l, dist, σ, bᵢ, b, s, k);
    double func_value = demo28(l, dist, σ, bᵢ, b, s, k);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}