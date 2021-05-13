/*
`C(x,y)` = (∑_n ∑_i c_n,i w_n,i R̂_n) / (∑_n ∑_i w_n,i R̂_n)

where

c ∈ ℝ^(f×s): the value of the Bayer pixel
w ∈ ℝ^(f×s): the local sample weight
R̂ ∈ ℝ^f: the local robustness
*/
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>

struct multi_frame_1ResultType {
    double C_left_parenthesis_x_comma_y_right_parenthesis;
    multi_frame_1ResultType(const double & C_left_parenthesis_x_comma_y_right_parenthesis)
    : C_left_parenthesis_x_comma_y_right_parenthesis(C_left_parenthesis_x_comma_y_right_parenthesis)
    {}
};

/**
 * multi_frame_1
 *
 * @param c  the value of the Bayer pixel
 * @param w  the local sample weight
 * @param R̂  the local robustness
 * @return C_left_parenthesis_x_comma_y_right_parenthesis
 */
multi_frame_1ResultType multi_frame_1(
    const Eigen::MatrixXd & c,
    const Eigen::MatrixXd & w,
    const Eigen::VectorXd & R̂)
{
    const long f = c.rows();
    const long s = c.cols();
    assert( c.rows() == f );
    assert( w.rows() == f );
    assert( w.cols() == s );
    assert( R̂.size() == f );

    double sum_0 = 0;
    for(int n=1; n<=w.rows(); n++){
        double sum_1 = 0;
        for(int i=1; i<=w.cols(); i++){
            sum_1 += c(n-1, i-1) * w(n-1, i-1) * R̂[n-1];
        }
        sum_0 += sum_1;
    }
    double sum_2 = 0;
    for(int n=1; n<=w.rows(); n++){
        double sum_3 = 0;
        for(int i=1; i<=w.cols(); i++){
            sum_3 += w(n-1, i-1) * R̂[n-1];
        }
        sum_2 += sum_3;
    }
    double C_left_parenthesis_x_comma_y_right_parenthesis = (sum_0) / double((sum_2));

    return multi_frame_1ResultType(C_left_parenthesis_x_comma_y_right_parenthesis);
}


void generateRandomData(Eigen::MatrixXd & c,
    Eigen::MatrixXd & w,
    Eigen::VectorXd & R̂)
{
    const int f = rand()%10;
    const int s = rand()%10;
    c = Eigen::MatrixXd::Random(f, s);
    w = Eigen::MatrixXd::Random(f, s);
    R̂ = Eigen::VectorXd::Random(f);
}


int main(int argc, char *argv[])
{
    srand((int)time(NULL));
    Eigen::MatrixXd c;
    Eigen::MatrixXd w;
    Eigen::VectorXd R̂;
    generateRandomData(c, w, R̂);
    multi_frame_1ResultType func_value = multi_frame_1(c, w, R̂);
    std::cout<<"return value:\n"<<func_value.C_left_parenthesis_x_comma_y_right_parenthesis<<std::endl;
    return 0;
}