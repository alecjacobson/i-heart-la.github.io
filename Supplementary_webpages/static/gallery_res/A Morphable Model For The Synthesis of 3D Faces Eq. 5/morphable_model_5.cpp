/*
E = 1/`σ_N`²`E_I` + ∑_j α_j²/`σ_S`_j² + ∑_j β_j²/`σ_T`_j²  + ∑_j (ρ_j-ρ̄_j)²/`σ_ρ`_j²

where

`σ_N` ∈ ℝ 
`E_I` ∈ ℝ
α_i ∈ ℝ
β_i ∈ ℝ
`σ_S`_i ∈ ℝ 
`σ_T`_i ∈ ℝ 
ρ_j ∈ ℝ 
ρ̄_j ∈ ℝ 
`σ_ρ`_j ∈ ℝ 
ā_i ∈ ℝ 
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct morphable_model_5ResultType {
    double E;
    morphable_model_5ResultType(const double & E)
    : E(E)
    {}
};

morphable_model_5ResultType morphable_model_5(
    const double & σ_N,
    const double & E_I,
    const std::vector<double> & α,
    const std::vector<double> & β,
    const std::vector<double> & σ_S,
    const std::vector<double> & σ_T,
    const std::vector<double> & ρ,
    const std::vector<double> & ρ̄,
    const std::vector<double> & σ_ρ,
    const std::vector<double> & ā)
{
    const long dim_0 = α.size();
    const long dim_1 = ρ.size();
    assert( β.size() == dim_0 );
    assert( σ_S.size() == dim_0 );
    assert( σ_T.size() == dim_0 );
    assert( ρ̄.size() == dim_1 );
    assert( σ_ρ.size() == dim_1 );
    assert( ā.size() == dim_0 );

    double sum_0 = 0;
    for(int j=1; j<=α.size(); j++){
        sum_0 += pow(α.at(j-1), 2) / double(pow(σ_S.at(j-1), 2));
    }
    double sum_1 = 0;
    for(int j=1; j<=β.size(); j++){
        sum_1 += pow(β.at(j-1), 2) / double(pow(σ_T.at(j-1), 2));
    }
    double sum_2 = 0;
    for(int j=1; j<=ρ.size(); j++){
        sum_2 += pow((ρ.at(j-1) - ρ̄.at(j-1)), 2) / double(pow(σ_ρ.at(j-1), 2));
    }
    double E = 1 / double(pow(σ_N, 2)) * E_I + sum_0 + sum_1 + sum_2;

    return morphable_model_5ResultType(E);
}


void generateRandomData(double & σ_N,
    double & E_I,
    std::vector<double> & α,
    std::vector<double> & β,
    std::vector<double> & σ_S,
    std::vector<double> & σ_T,
    std::vector<double> & ρ,
    std::vector<double> & ρ̄,
    std::vector<double> & σ_ρ,
    std::vector<double> & ā)
{
    σ_N = rand() % 10;
    E_I = rand() % 10;
    const int dim_0 = rand()%10;
    const int dim_1 = rand()%10;
    α.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        α[i] = rand() % 10;
    }
    β.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        β[i] = rand() % 10;
    }
    σ_S.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        σ_S[i] = rand() % 10;
    }
    σ_T.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        σ_T[i] = rand() % 10;
    }
    ρ.resize(dim_1);
    for(int i=0; i<dim_1; i++){
        ρ[i] = rand() % 10;
    }
    ρ̄.resize(dim_1);
    for(int i=0; i<dim_1; i++){
        ρ̄[i] = rand() % 10;
    }
    σ_ρ.resize(dim_1);
    for(int i=0; i<dim_1; i++){
        σ_ρ[i] = rand() % 10;
    }
    ā.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        ā[i] = rand() % 10;
    }
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    double σ_N;
    double E_I;
    std::vector<double> α;
    std::vector<double> β;
    std::vector<double> σ_S;
    std::vector<double> σ_T;
    std::vector<double> ρ;
    std::vector<double> ρ̄;
    std::vector<double> σ_ρ;
    std::vector<double> ā;
    generateRandomData(σ_N, E_I, α, β, σ_S, σ_T, ρ, ρ̄, σ_ρ, ā);
    morphable_model_5ResultType func_value = morphable_model_5(σ_N, E_I, α, β, σ_S, σ_T, ρ, ρ̄, σ_ρ, ā);
    std::cout<<"return value:\n"<<func_value.E<<std::endl;
    return 0;
}