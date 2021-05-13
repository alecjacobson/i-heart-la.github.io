#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <igl/readMSH.h>


void circumcenter(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, Eigen::Vector3d& cc)
{
    Eigen::Matrix<double, 3, 1> _l_0;
    _l_0 << pow((b - c).lpNorm<2>(), 2), pow((a - c).lpNorm<2>(), 2), pow((a - b).lpNorm<2>(), 2);
    Eigen::Matrix<double, 3, 1> l = _l_0;

    Eigen::Matrix<double, 3, 1> _ba_0;
    _ba_0 << l[1-1] * (l[2-1] + l[3-1] - l[1-1]), l[2-1] * (l[3-1] + l[1-1] - l[2-1]), l[3-1] * (l[1-1] + l[2-1] - l[3-1]);
    Eigen::Matrix<double, 3, 1> ba = _ba_0;

    cc = 1 / double((ba[1-1] + ba[2-1] + ba[3-1])) * (ba[1-1] * a + ba[2-1] * b + ba[3-1] * c);
}

void circumcenter(const Eigen::Matrix<double, 4, 3>& t, Eigen::Vector3d& c)
{
    Eigen::Matrix3d A;
    Eigen::Vector3d b;
    
    const double n0 = t.row(0).squaredNorm();
    
    for(int k = 0; k < 3; ++k)
    {
        A.row(k) = t.row(k + 1) - t.row(0);
        b(k) = t.row(k + 1).squaredNorm() - n0;
    }
    
    c = 0.5 * A.fullPivHouseholderQr().solve(b);
}

double volume(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, const Eigen::Vector3d& d)
{
    Eigen::Matrix3d A;
    A.col(0) = b - a;
    A.col(1) = c - a;
    A.col(2) = d - a;
    
    return A.determinant() / 6.;
}

void dualLaplace(const Eigen::MatrixXd& V, const Eigen::MatrixXi& T, Eigen::SparseMatrix<double>& L, Eigen::SparseMatrix<double>& M)
{
    const size_t nt = T.rows();
    const size_t nv = V.rows();
    
    const int turn[4][4]
    {
        {-1, 2, 3, 1},
        {3, -1, 0, 2},
        {1, 3, -1, 0},
        {2, 0, 1, -1}
    };

    auto getTet = [&](const int i, Eigen::Matrix<double, 4, 3>& t)
    {
        for(int k = 0; k < 4; ++k)
        {
            t.row(k) = V.row(T(i, k));
        }
    };
    
    std::vector<Eigen::Triplet<double>> tripL, tripM;
    
    Eigen::Vector3d cc;
    Eigen::Matrix<double, 4, 3> t;
    
    for(int k = 0; k < nt; ++k)
    {
        getTet(k, t);
        circumcenter(t, cc);
        
        for(int i = 0; i < 4; ++i)
        {
            for(int j = 0; j < 4; ++j)
            {
                if(i != j)
                {
                    Eigen::Vector3d cf;
                    circumcenter(t.row(i), t.row(j), t.row(turn[i][j]), cf);
                    
                    const Eigen::Vector3d ce = 0.5 * (t.row(i) + t.row(j));
                    
                    const double vol = volume(t.row(i), ce, cf, cc);
                    const double wij = 6. * vol / (t.row(i) - t.row(j)).squaredNorm();
                    
                    tripL.emplace_back(T(k, i), T(k, j), wij);
                    tripL.emplace_back(T(k, j), T(k, i), wij);
                    
                    tripL.emplace_back(T(k, i), T(k, i), -wij);
                    tripL.emplace_back(T(k, j), T(k, j), -wij);
                    
                    tripM.emplace_back(T(k, i), T(k, i), vol);
                    tripM.emplace_back(T(k, j), T(k, j), vol);
                }
            }
        }
    }
    
    L.resize(nv, nv);
    M.resize(nv, nv);
    
    L.setFromTriplets(tripL.begin(), tripL.end());
    M.setFromTriplets(tripM.begin(), tripM.end());
}

int main(int argc, const char * argv[])
{
    if(argc == 2)
    {
        Eigen::SparseMatrix<double> L, M;
        Eigen::MatrixXd V;
        Eigen::MatrixXi T;

        if(igl::readMSH(argv[1], V, T))
        {
            dualLaplace(V, T, L, M);
            Eigen::saveMarket(L, "L.mtx");
            Eigen::saveMarket(M, "M.mtx");
        }
    }
    
    return 0;
}
