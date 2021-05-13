/*
`G_σ(s^k_i)` = ∑_j l_j exp(-dist(`bᵢ`, b_j)²/(2σ²)) `s^k`_j

where

l_j ∈ ℝ : the length of bj
dist: ℝ^n, ℝ^n → ℝ : measures the geodesic distance 
σ ∈ ℝ
`bᵢ` ∈ ℝ^n
b_j ∈ ℝ^n
`s^k`_j ∈ ℝ^n : direction vector
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct atlas_refinement_3ResultType {
    Eigen::VectorXd G_σ_left_parenthesis_s_circumflex_accent_k_i_right_parenthesis;
    atlas_refinement_3ResultType(const Eigen::VectorXd & G_σ_left_parenthesis_s_circumflex_accent_k_i_right_parenthesis)
    : G_σ_left_parenthesis_s_circumflex_accent_k_i_right_parenthesis(G_σ_left_parenthesis_s_circumflex_accent_k_i_right_parenthesis)
    {}
};

/**
 * atlas_refinement_3
 *
 * @param l  the length of bj
 * @param dist  ℝ^n, ℝ^n → ℝ : measures the geodesic distance 
 * @param s_circumflex_accent_k  direction vector
 * @return G_σ_left_parenthesis_s_circumflex_accent_k_i_right_parenthesis
 */
atlas_refinement_3ResultType atlas_refinement_3(
    const std::vector<double> & l,
    const std::function<double(Eigen::VectorXd, Eigen::VectorXd)> & dist,
    const double & σ,
    const Eigen::VectorXd & bᵢ,
    const std::vector<Eigen::VectorXd> & b,
    const std::vector<Eigen::VectorXd> & s_circumflex_accent_k)
{
    const long dim_0 = l.size();
    const long n = bᵢ.size();
    assert( b.size() == dim_0 );
    for( const auto& el : b ) {
        assert( el.size() == n );
    }
    assert( s_circumflex_accent_k.size() == dim_0 );
    for( const auto& el : s_circumflex_accent_k ) {
        assert( el.size() == n );
    }

    Eigen::MatrixXd sum_0 = Eigen::MatrixXd::Zero(n, 1);
    for(int j=1; j<=s_circumflex_accent_k.size(); j++){
        sum_0 += l.at(j-1) * exp(-pow(dist(bᵢ, b.at(j-1)), 2) / double((2 * pow(σ, 2)))) * s_circumflex_accent_k.at(j-1);
    }
    Eigen::VectorXd G_σ_left_parenthesis_s_circumflex_accent_k_i_right_parenthesis = sum_0;

    return atlas_refinement_3ResultType(G_σ_left_parenthesis_s_circumflex_accent_k_i_right_parenthesis);
}


void generateRandomData(std::vector<double> & l,
    std::function<double(Eigen::VectorXd, Eigen::VectorXd)> & dist,
    double & σ,
    Eigen::VectorXd & bᵢ,
    std::vector<Eigen::VectorXd> & b,
    std::vector<Eigen::VectorXd> & s_circumflex_accent_k)
{
    σ = rand() % 10;
    const int dim_0 = rand()%10;
    const int n = rand()%10;
    l.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        l[i] = rand() % 10;
    }
    dist = [](Eigen::VectorXd, Eigen::VectorXd)->double{
        return rand() % 10;
    };
    bᵢ = Eigen::VectorXd::Random(n);
    b.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        b[i] = Eigen::VectorXd::Random(n);
    }
    s_circumflex_accent_k.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        s_circumflex_accent_k[i] = Eigen::VectorXd::Random(n);
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
    std::vector<Eigen::VectorXd> s_circumflex_accent_k;
    generateRandomData(l, dist, σ, bᵢ, b, s_circumflex_accent_k);
    atlas_refinement_3ResultType func_value = atlas_refinement_3(l, dist, σ, bᵢ, b, s_circumflex_accent_k);
    std::cout<<"return value:\n"<<func_value.G_σ_left_parenthesis_s_circumflex_accent_k_i_right_parenthesis<<std::endl;
    return 0;
}