/*

E = 1/`σ_N`^2`E_I` + sum_(j for j>1) α_j²/`σ_S`_j^2  + sum_(j for j>1) β_j²/`σ_T`_j^2   + sum_j (ρ_j-`ρ_bar`_j)²/`σ_ρ`_j^2 

where

`σ_N`: ℝ 
`E_I`: ℝ
α_i : ℝ
β_i : ℝ
`σ_S`_i: ℝ 
`σ_T`_i: ℝ 
ρ_i: ℝ 
`ρ_bar`_i: ℝ 
`σ_ρ`_i: ℝ 
ā_i: ℝ 
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo24
 *
 * @param σ_N  ℝ
 * @param E_I  ℝ
 * @param α  ℝ
 * @param β  ℝ
 * @param σ_S  ℝ
 * @param σ_T  ℝ
 * @param ρ  ℝ
 * @param ρ_bar  ℝ
 * @param σ_ρ  ℝ
 * @param ā  ℝ
 * @return E
 */
double demo24(
    const double & σ_N,
    const double & E_I,
    const std::vector<double> & α,
    const std::vector<double> & β,
    const std::vector<double> & σ_S,
    const std::vector<double> & σ_T,
    const std::vector<double> & ρ,
    const std::vector<double> & ρ_bar,
    const std::vector<double> & σ_ρ,
    const std::vector<double> & ā)
{
    const long _dim_0 = α.size();
    assert( α.size() == _dim_0 );
    assert( β.size() == _dim_0 );
    assert( σ_S.size() == _dim_0 );
    assert( σ_T.size() == _dim_0 );
    assert( ρ.size() == _dim_0 );
    assert( ρ_bar.size() == _dim_0 );
    assert( σ_ρ.size() == _dim_0 );
    assert( ā.size() == _dim_0 );

    double _sum_0 = 0;
    for(int j=1; j<=α.size(); j++){
        if(j > 1){
            _sum_0 += pow(α.at(j-1), 2) / pow(σ_S.at(j-1), 2);
        }
    }
    double _sum_1 = 0;
    for(int j=1; j<=β.size(); j++){
        if(j > 1){
            _sum_1 += pow(β.at(j-1), 2) / pow(σ_T.at(j-1), 2);
        }
    }
    double _sum_2 = 0;
    for(int j=1; j<=ρ.size(); j++){
        _sum_2 += pow((ρ.at(j-1) - ρ_bar.at(j-1)), 2) / pow(σ_ρ.at(j-1), 2);
    }
    double E = 1 / pow(σ_N, 2) * E_I + _sum_0 + _sum_1 + _sum_2;

    return E;
}


void generateRandomData(double & σ_N,
    double & E_I,
    std::vector<double> & α,
    std::vector<double> & β,
    std::vector<double> & σ_S,
    std::vector<double> & σ_T,
    std::vector<double> & ρ,
    std::vector<double> & ρ_bar,
    std::vector<double> & σ_ρ,
    std::vector<double> & ā)
{
    σ_N = rand() % 10;
    E_I = rand() % 10;
    const int _dim_0 = rand()%10;
    α.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        α[i] = rand() % 10;
    }
    β.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        β[i] = rand() % 10;
    }
    σ_S.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        σ_S[i] = rand() % 10;
    }
    σ_T.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        σ_T[i] = rand() % 10;
    }
    ρ.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        ρ[i] = rand() % 10;
    }
    ρ_bar.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        ρ_bar[i] = rand() % 10;
    }
    σ_ρ.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        σ_ρ[i] = rand() % 10;
    }
    ā.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        ā[i] = rand() % 10;
    }
}


int main(int argc, char *argv[])
{
    double σ_N;
    double E_I;
    std::vector<double> α;
    std::vector<double> β;
    std::vector<double> σ_S;
    std::vector<double> σ_T;
    std::vector<double> ρ;
    std::vector<double> ρ_bar;
    std::vector<double> σ_ρ;
    std::vector<double> ā;
    generateRandomData(σ_N, E_I, α, β, σ_S, σ_T, ρ, ρ_bar, σ_ρ, ā);
    double func_value = demo24(σ_N, E_I, α, β, σ_S, σ_T, ρ, ρ_bar, σ_ρ, ā);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}