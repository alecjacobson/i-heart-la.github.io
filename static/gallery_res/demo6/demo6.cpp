/*
sum_i f_i²p_i - (sum_i f_i p_i)²

where

f_i: ℝ  
p_i: ℝ  
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo6
 *
 * @param f  ℝ
 * @param p  ℝ
 * @return ret
 */
double demo6(
    const std::vector<double> & f,
    const std::vector<double> & p)
{
    const long _dim_0 = f.size();
    assert( f.size() == _dim_0 );
    assert( p.size() == _dim_0 );

    double _sum_0 = 0;
    for(int i=1; i<=f.size(); i++){
        _sum_0 += pow(f.at(i-1), 2) * p.at(i-1);
    }
    double _sum_1 = 0;
    for(int i=1; i<=f.size(); i++){
        _sum_1 += f.at(i-1) * p.at(i-1);
    }
    double ret = _sum_0 - pow((_sum_1), 2);
    return ret;
}


void generateRandomData(std::vector<double> & f,
    std::vector<double> & p)
{
    const int _dim_0 = rand()%10;
    f.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        f[i] = rand() % 10;
    }
    p.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        p[i] = rand() % 10;
    }
}


int main(int argc, char *argv[])
{
    std::vector<double> f;
    std::vector<double> p;
    generateRandomData(f, p);
    double func_value = demo6(f, p);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}