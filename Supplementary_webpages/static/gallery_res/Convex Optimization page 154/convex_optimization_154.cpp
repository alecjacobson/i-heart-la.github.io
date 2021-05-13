/*
∑_i f_i²p_i - (∑_i f_i p_i)²

where

f_i ∈ ℝ  
p_i ∈ ℝ  
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct convex_optimization_154ResultType {
    double ret;
    convex_optimization_154ResultType(const double & ret)
    : ret(ret)
    {}
};

convex_optimization_154ResultType convex_optimization_154(
    const std::vector<double> & f,
    const std::vector<double> & p)
{
    const long dim_0 = f.size();
    assert( p.size() == dim_0 );

    double sum_0 = 0;
    for(int i=1; i<=p.size(); i++){
        sum_0 += pow(f.at(i-1), 2) * p.at(i-1);
    }
    double sum_1 = 0;
    for(int i=1; i<=p.size(); i++){
        sum_1 += f.at(i-1) * p.at(i-1);
    }
    double ret = sum_0 - pow((sum_1), 2);
    return convex_optimization_154ResultType(ret);
}


void generateRandomData(std::vector<double> & f,
    std::vector<double> & p)
{
    const int dim_0 = rand()%10;
    f.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        f[i] = rand() % 10;
    }
    p.resize(dim_0);
    for(int i=0; i<dim_0; i++){
        p[i] = rand() % 10;
    }
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    std::vector<double> f;
    std::vector<double> p;
    generateRandomData(f, p);
    convex_optimization_154ResultType func_value = convex_optimization_154(f, p);
    std::cout<<"return value:\n"<<func_value.ret<<std::endl;
    return 0;
}