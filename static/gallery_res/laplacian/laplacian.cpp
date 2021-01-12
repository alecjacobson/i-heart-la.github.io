/*
from trigonometry: cot

L_i,j = { cot(α_ij) + cot(β_ij) if j ∈ N(i)

L_i,i = -sum_(k for k != i) L_i,k

where

L: ℝ^(n×n)
α: ℝ^(n×n)
β: ℝ^(n×n)
N: ℤ -> {ℤ}
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * laplacian
 *
 * @param α  ℝ^(n×n)
 * @param β  ℝ^(n×n)
 * @param N  ℤ -> {ℤ}
 * @return L
 */
Eigen::SparseMatrix<double> laplacian(
    const Eigen::MatrixXd & α,
    const Eigen::MatrixXd & β,
    const std::function<std::set<std::tuple< int > >(int)> & N)
{
    const long n = β.cols();
    assert( α.rows() == n );
    assert( α.cols() == n );
    assert( β.rows() == n );
    assert( β.cols() == n );

    Eigen::SparseMatrix<double> L(n, n);
    std::vector<Eigen::Triplet<double> > tripletList_L;
    for( int i=1; i<=n; i++){
        for( int j=1; j<=n; j++){
            if(N(i).find(std::tuple< int >(j)) != N(i).end()){
                tripletList_L.push_back(Eigen::Triplet<double>(i-1, j-1, 1/tan(α(i-1, j-1)) + 1/tan(β(i-1, j-1))));
            }
        }
    }
    L.setFromTriplets(tripletList_L.begin(), tripletList_L.end());


    for( int i=1; i<=n; i++){
        double _sum_0 = 0;
        for(int k=1; k<=L.cols(); k++){
            if(k != i){
                _sum_0 += L.coeff(i-1, k-1);
            }
        }
        tripletList_L.push_back(Eigen::Triplet<double>(i-1, i-1, -_sum_0));
    }
    L.setFromTriplets(tripletList_L.begin(), tripletList_L.end());


    return L;
}


void generateRandomData(Eigen::MatrixXd & α,
    Eigen::MatrixXd & β,
    std::function<std::set<std::tuple< int > >(int)> & N)
{
    const int n = rand()%10;
    α = Eigen::MatrixXd::Random(n, n);
    β = Eigen::MatrixXd::Random(n, n);
    N = [](int)->std::set<std::tuple< int > >{
        std::set<std::tuple< int > > tmp;
        const int tmp_0 = rand()%10;
        for(int i=0; i<tmp_0; i++){
            tmp.insert(std::make_tuple(rand()%10));
        }
        return tmp;
    };
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    Eigen::MatrixXd α;
    Eigen::MatrixXd β;
    std::function<std::set<std::tuple< int > >(int)> N;
    generateRandomData(α, β, N);
    Eigen::SparseMatrix<double> func_value = laplacian(α, β, N);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}