/*
L_i,j = { w_i,j if (i,j) ∈ E
          0 otherwise

L_i,i = -sum_(l for l ≠ i) L_i,l

where

L ∈ ℝ^(n×n)
w ∈ ℝ^(n×n): edge weight matrix
E ∈ {ℤ²} index: edges
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct course_parameterizationResultType {
    Eigen::SparseMatrix<double> L;
    course_parameterizationResultType(const Eigen::SparseMatrix<double> & L)
    : L(L)
    {}
};

/**
 * course_parameterization
 *
 * @param w  edge weight matrix
 * @param E  edges
 * @return L
 */
course_parameterizationResultType course_parameterization(
    const Eigen::MatrixXd & w,
    const std::set<std::tuple< int, int > > & E)
{
    const long n = w.cols();
    assert( w.rows() == n );

    Eigen::SparseMatrix<double> L(n, n);
    std::vector<Eigen::Triplet<double> > tripletList_L;
    for( int i=1; i<=n; i++){
        for( int j=1; j<=n; j++){
            if(E.find(std::tuple< int, int >(i-1, j-1)) != E.end()){
                tripletList_L.push_back(Eigen::Triplet<double>(i-1, j-1, w(i-1, j-1)));
            }
        }
    }
    L.setFromTriplets(tripletList_L.begin(), tripletList_L.end());

    for( int i=1; i<=n; i++){
        double sum_0 = 0;
        for(int l=1; l<=L.cols(); l++){
            if(l != i){
                sum_0 += L.coeff(i-1, l-1);
            }
        }
        tripletList_L.push_back(Eigen::Triplet<double>(i-1, i-1, -sum_0));
    }
    L.setFromTriplets(tripletList_L.begin(), tripletList_L.end());

    return course_parameterizationResultType(L);
}


void generateRandomData(Eigen::MatrixXd & w,
    std::set<std::tuple< int, int > > & E)
{
    const int n = rand()%10;
    w = Eigen::MatrixXd::Random(n, n);
    const int dim_1 = rand()%10;
    for(int i=0; i<dim_1; i++){
        E.insert(std::make_tuple(rand()%10, rand()%10));
    }
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    Eigen::MatrixXd w;
    std::set<std::tuple< int, int > > E;
    generateRandomData(w, E);
    course_parameterizationResultType func_value = course_parameterization(w, E);
    std::cout<<"return value:\n"<<func_value.L<<std::endl;
    return 0;
}