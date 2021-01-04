/*
`G_σ(s_i^k)` = sum_j l_j exp(-`dist`(`b_i`, b_j)/(2`σ`^2))(s_j)^k

where

l_j: ℝ : the length of b_j
`dist`: ℝ, ℝ -> ℝ : measures the geodesic distance between the centers of b_i and b_j along the boundary
`σ`: ℝ
`b_i`: ℝ
b_j: ℝ
s_j: ℝ : unit direction vector of b_i
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
 * @param l  ℝ : the length of b_j
 * @param dist  ℝ, ℝ -> ℝ : measures the geodesic distance between the centers of b_i and b_j along the boundary
 * @param σ  ℝ
 * @param b_i  ℝ
 * @param b  ℝ
 * @param s  ℝ : unit direction vector of b_i
 * @param k  ℝ : iteration number
 * @return G_σ_left_parenthesis_s_i_circumflex_accent_k_right_parenthesis
 */
double demo28(
    const std::vector<double> & l,
    const std::function<double(double, double)> & dist,
    const double & σ,
    const double & b_i,
    const std::vector<double> & b,
    const std::vector<double> & s,
    const double & k)
{
    const long _dim_0 = l.size();
    assert( l.size() == _dim_0 );
    assert( b.size() == _dim_0 );
    assert( s.size() == _dim_0 );

    double _sum_0 = 0;
    for(int j=1; j<=b.size(); j++){
        _sum_0 += l.at(j-1) * exp(-dist(b_i, b.at(j-1)) / (2 * pow(σ, 2))) * pow((s.at(j-1)), k);
    }
    double G_σ_left_parenthesis_s_i_circumflex_accent_k_right_parenthesis = _sum_0;

    return G_σ_left_parenthesis_s_i_circumflex_accent_k_right_parenthesis;
}


void generateRandomData(std::vector<double> & l,
    std::function<double(double, double)> & dist,
    double & σ,
    double & b_i,
    std::vector<double> & b,
    std::vector<double> & s,
    double & k)
{
    σ = rand() % 10;
    b_i = rand() % 10;
    k = rand() % 10;
    const int _dim_0 = rand()%10;
    l.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        l[i] = rand() % 10;
    }
    dist = [](double, double)->double{
        return rand() % 10;
    };
    b.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        b[i] = rand() % 10;
    }
    s.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        s[i] = rand() % 10;
    }
}


int main(int argc, char *argv[])
{
    std::vector<double> l;
    std::function<double(double, double)> dist;
    double σ;
    double b_i;
    std::vector<double> b;
    std::vector<double> s;
    double k;
    generateRandomData(l, dist, σ, b_i, b, s, k);
    double func_value = demo28(l, dist, σ, b_i, b, s, k);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}