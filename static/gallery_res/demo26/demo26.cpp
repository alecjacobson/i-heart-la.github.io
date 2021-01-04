/*
`C(x,y)` = (sum_n sum_i c_n,i w_n,i R_n) / (sum_n sum_i w_n,i R_n)

where

c: ℝ^(x×y): the value of the Bayer pixel
w: ℝ^(x×y): the local sample weight
R: ℝ^x: the local robustness
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

/**
 * demo26
 *
 * @param c  ℝ^(x×y): the value of the Bayer pixel
 * @param w  ℝ^(x×y): the local sample weight
 * @param R  ℝ^x: the local robustness
 * @return C_left_parenthesis_x_comma_y_right_parenthesis
 */
double demo26(
    const Eigen::MatrixXd & c,
    const Eigen::MatrixXd & w,
    const Eigen::VectorXd & R)
{
    const long x = R.size();
    const long y = w.cols();
    assert( c.rows() == x );
    assert( c.cols() == y );
    assert( w.rows() == x );
    assert( w.cols() == y );
    assert( R.size() == x );

    double _sum_0 = 0;
    for(int n=1; n<=R.size(); n++){
        double _sum_1 = 0;
        for(int i=1; i<=w.cols(); i++){
            _sum_1 += c(n-1, i-1) * w(n-1, i-1) * R(n-1);
        }
        _sum_0 += _sum_1;
    }
    double _sum_2 = 0;
    for(int n=1; n<=R.size(); n++){
        double _sum_3 = 0;
        for(int i=1; i<=w.cols(); i++){
            _sum_3 += w(n-1, i-1) * R(n-1);
        }
        _sum_2 += _sum_3;
    }
    double C_left_parenthesis_x_comma_y_right_parenthesis = (_sum_0) / (_sum_2);

    return C_left_parenthesis_x_comma_y_right_parenthesis;
}


void generateRandomData(Eigen::MatrixXd & c,
    Eigen::MatrixXd & w,
    Eigen::VectorXd & R)
{
    const int x = rand()%10;
    const int y = rand()%10;
    c = Eigen::MatrixXd::Random(x, y);
    w = Eigen::MatrixXd::Random(x, y);
    R = Eigen::VectorXd::Random(x);
}


int main(int argc, char *argv[])
{
    Eigen::MatrixXd c;
    Eigen::MatrixXd w;
    Eigen::VectorXd R;
    generateRandomData(c, w, R);
    double func_value = demo26(c, w, R);
    std::cout<<"func_value:\n"<<func_value<<std::endl;
    return 0;
}