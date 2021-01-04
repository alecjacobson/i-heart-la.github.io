/*
given
p_i: ℝ^3: points on lines
d_i: ℝ^3: unit directions along lines

k_i = (p_i - (p_i⋅d_i)d_i)
a_i = (1,0,0) - d_i,1 d_i
b_i = (0,1,0) - d_i,2 d_i
c_i = (0,0,1) - d_i,3 d_i

 
M = [ ∑_i( a_i,1 - d_i,1 (d_i⋅a_i) )    ∑_i( a_i,2 - d_i,2 (d_i⋅a_i) )    ∑_i( a_i,3 - d_i,3 (d_i⋅a_i) )
      ∑_i( b_i,1 - d_i,1 (d_i⋅b_i) )    ∑_i( b_i,2 - d_i,2 (d_i⋅b_i) )    ∑_i( b_i,3 - d_i,3 (d_i⋅b_i) )
      ∑_i( c_i,1 - d_i,1 (d_i⋅c_i) )    ∑_i( c_i,2 - d_i,2 (d_i⋅c_i) )    ∑_i( c_i,3 - d_i,3 (d_i⋅c_i) ) ]

r = [ ∑_i( k_i⋅a_i )
      ∑_i( k_i⋅b_i )
      ∑_i( k_i⋅c_i ) ]

q = M^(-1) r
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo30
 *
 * @param p  ℝ^3: points on lines
 * @param d  ℝ^3: unit directions along lines
 * @return q
 */
Eigen::Matrix<double, 3, 1> demo30(
    const std::vector<Eigen::Matrix<double, 3, 1>> & p,
    const std::vector<Eigen::Matrix<double, 3, 1>> & d)
{
    const long _dim_0 = p.size();
    std::vector<Eigen::Matrix<double, 3, 1>> k(_dim_0);
    for( int i=1; i<=_dim_0; i++){
        k.at(i-1) = (p.at(i-1) - ((p.at(i-1)).dot(d.at(i-1))) * d.at(i-1));
    }


    std::vector<Eigen::Matrix<double, 3, 1>> a(_dim_0);
    for( int i=1; i<=_dim_0; i++){
        Eigen::Matrix<double, 3, 1> _a_i_0;
        _a_i_0 << 1, 0, 0;
        a.at(i-1) = _a_i_0 - d.at(i-1)(1-1) * d.at(i-1);
    }


    std::vector<Eigen::Matrix<double, 3, 1>> b(_dim_0);
    for( int i=1; i<=_dim_0; i++){
        Eigen::Matrix<double, 3, 1> _b_i_0;
        _b_i_0 << 0, 1, 0;
        b.at(i-1) = _b_i_0 - d.at(i-1)(2-1) * d.at(i-1);
    }


    std::vector<Eigen::Matrix<double, 3, 1>> c(_dim_0);
    for( int i=1; i<=_dim_0; i++){
        Eigen::Matrix<double, 3, 1> _c_i_0;
        _c_i_0 << 0, 0, 1;
        c.at(i-1) = _c_i_0 - d.at(i-1)(3-1) * d.at(i-1);
    }


    double _sum_0 = 0;
    for(int i=1; i<=d.size(); i++){
        _sum_0 += (a.at(i-1)(1-1) - d.at(i-1)(1-1) * ((d.at(i-1)).dot(a.at(i-1))));
    }
    double _sum_1 = 0;
    for(int i=1; i<=d.size(); i++){
        _sum_1 += (a.at(i-1)(2-1) - d.at(i-1)(2-1) * ((d.at(i-1)).dot(a.at(i-1))));
    }
    double _sum_2 = 0;
    for(int i=1; i<=d.size(); i++){
        _sum_2 += (a.at(i-1)(3-1) - d.at(i-1)(3-1) * ((d.at(i-1)).dot(a.at(i-1))));
    }
    double _sum_3 = 0;
    for(int i=1; i<=d.size(); i++){
        _sum_3 += (b.at(i-1)(1-1) - d.at(i-1)(1-1) * ((d.at(i-1)).dot(b.at(i-1))));
    }
    double _sum_4 = 0;
    for(int i=1; i<=b.size(); i++){
        _sum_4 += (b.at(i-1)(2-1) - d.at(i-1)(2-1) * ((d.at(i-1)).dot(b.at(i-1))));
    }
    double _sum_5 = 0;
    for(int i=1; i<=b.size(); i++){
        _sum_5 += (b.at(i-1)(3-1) - d.at(i-1)(3-1) * ((d.at(i-1)).dot(b.at(i-1))));
    }
    double _sum_6 = 0;
    for(int i=1; i<=d.size(); i++){
        _sum_6 += (c.at(i-1)(1-1) - d.at(i-1)(1-1) * ((d.at(i-1)).dot(c.at(i-1))));
    }
    double _sum_7 = 0;
    for(int i=1; i<=d.size(); i++){
        _sum_7 += (c.at(i-1)(2-1) - d.at(i-1)(2-1) * ((d.at(i-1)).dot(c.at(i-1))));
    }
    double _sum_8 = 0;
    for(int i=1; i<=d.size(); i++){
        _sum_8 += (c.at(i-1)(3-1) - d.at(i-1)(3-1) * ((d.at(i-1)).dot(c.at(i-1))));
    }
    Eigen::Matrix<double, 3, 3> _M_0;
    _M_0 << _sum_0, _sum_1, _sum_2,
    _sum_3, _sum_4, _sum_5,
    _sum_6, _sum_7, _sum_8;
    Eigen::Matrix<double, 3, 3> M = _M_0;

    double _sum_9 = 0;
    for(int i=1; i<=k.size(); i++){
        _sum_9 += ((k.at(i-1)).dot(a.at(i-1)));
    }
    double _sum_10 = 0;
    for(int i=1; i<=k.size(); i++){
        _sum_10 += ((k.at(i-1)).dot(b.at(i-1)));
    }
    double _sum_11 = 0;
    for(int i=1; i<=k.size(); i++){
        _sum_11 += ((k.at(i-1)).dot(c.at(i-1)));
    }
    Eigen::Matrix<double, 3, 1> _r_0;
    _r_0 << _sum_9,
    _sum_10,
    _sum_11;
    Eigen::Matrix<double, 3, 1> r = _r_0;

    Eigen::Matrix<double, 3, 1> q = M.inverse() * r;

    return q;
}


void generateRandomData(std::vector<Eigen::Matrix<double, 3, 1>> & p,
    std::vector<Eigen::Matrix<double, 3, 1>> & d)
{
    const int _dim_0 = rand()%10;
    p.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        p[i] = Eigen::VectorXd::Random(3);
    }
    d.resize(_dim_0);
    for(int i=0; i<_dim_0; i++){
        d[i] = Eigen::VectorXd::Random(3);
    }
}


int main(int argc, char *argv[])
{
    std::vector<Eigen::Matrix<double, 3, 1>> p;
    std::vector<Eigen::Matrix<double, 3, 1>> d;
    generateRandomData(p, d);
    Eigen::Matrix<double, 3, 1> func_value = demo30(p, d);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}