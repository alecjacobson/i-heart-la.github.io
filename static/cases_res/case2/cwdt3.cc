#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Regular_triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Regular_triangulation_face_base_2.h>
#include <CGAL/Regular_triangulation_vertex_base_2.h>

#include <CGAL/linear_least_squares_fitting_3.h>

#include <CGAL/algorithm.h>

#include <Eigen/Dense>

#include <igl/read_triangle_mesh.h>
#include <igl/bfs_orient.h>
#include <igl/winding_number.h>

#include <cassert>
#include <vector>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cstring>

#include <chrono>

#include "gurobi_c++.h"

#include "CLI11.hpp"

#define NODEBUG_OUTPUT
#ifdef DEBUG_OUTPUT
#define debug(x) std::cerr << x << std::endl
#else
#define debug(x) \
    {            \
        ;        \
    }
#endif

struct FaceInfo2
{
    FaceInfo2() {}
    int nesting_level;
    bool in_domain()
    {
        return nesting_level % 2 == 1;
    }
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT Weight;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Weighted_point_3 Weighted_point;

typedef CGAL::Regular_triangulation_vertex_base_3<K> Vb0;
typedef CGAL::Triangulation_vertex_base_with_info_3<int, K, Vb0> Vb;
typedef CGAL::Regular_triangulation_cell_base_3<K> Cb0;
typedef CGAL::Triangulation_cell_base_with_info_3<int, K, Cb0> Cb;

typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Regular_triangulation_3<K, Tds> Rt;

typedef Rt::Vertex_iterator Vertex_iterator;
typedef Rt::Finite_vertices_iterator Finite_vertices_iterator;
typedef Rt::Finite_cells_iterator Finite_cells_iterator;
typedef Rt::Vertex_handle VH;
typedef Rt::Cell_handle CH;

typedef CGAL::Triangulation_vertex_base_with_info_2<int, K> Vb2;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K> Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<K, Fbb> Fb2;
typedef CGAL::Triangulation_data_structure_2<Vb2, Fb2> Tds2;
typedef CGAL::No_intersection_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds2, Itag> CDT;

typedef CGAL::Regular_triangulation_vertex_base_2<K> Vbr0;
typedef CGAL::Triangulation_vertex_base_with_info_2<int, K, Vbr0> Vb2r;
typedef CGAL::Regular_triangulation_face_base_2<K> Cbr0;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K, Cbr0> Cb2r;

typedef CGAL::Triangulation_data_structure_2<Vb2r, Cb2r> Tds2r;
typedef CGAL::Regular_triangulation_2<K, Tds2r> Rt2;


Rt T;

GRBEnv *env;

// edge connectivity
struct edge
{
    std::array<int, 2> v;
    bool operator==(edge const &o) const
    {
        return ((v[0] == o.v[0]) || (v[0] == o.v[1])) && ((v[1] == o.v[0]) || (v[1] == o.v[1]));
    }
    int &operator[](int idx)
    {
        return v[idx];
    }
    const int &operator[](int idx) const
    {
        return v[idx];
    }
};

template <>
struct std::hash<edge>
{
    size_t operator()(edge const &key) const
    {
        return (size_t)key[0] * (size_t)key[1];
    }
};

struct triangle // assumes index elements are different!
{
    std::array<int, 3> v;
    bool operator==(triangle const &o) const
    {
        return ((v[0] == o.v[0]) || (v[0] == o.v[1]) || (v[0] == o.v[2])) &&
               ((v[1] == o.v[0]) || (v[1] == o.v[1]) || (v[1] == o.v[2])) &&
               ((v[2] == o.v[0]) || (v[2] == o.v[1]) || (v[2] == o.v[2]));
    }
    int &operator[](int idx)
    {
        return v[idx];
    }
    const int &operator[](int idx) const
    {
        return v[idx];
    }
};

template <>
struct std::hash<triangle>
{
    size_t operator()(triangle const &key) const
    {
        return (size_t)key[0] * (size_t)key[1] * (size_t)key[2];
    }
};

struct poly
{
    std::list<edge> be;
    std::vector<int> iv;
    std::vector<triangle> tri;
    K::Plane_3 plane;
};

std::vector<Point> V;
std::vector<double> W;
std::vector<poly> P; // constrained poly faces
int nv, np, splits;
bool minimize_heights,write_missing,write_max_constraints,write_witness;

std::unordered_map<edge, std::vector<int>> E2P;

std::vector<VH> id2vh;
std::vector<bool> ve;

bool gabriel;

// data strcture for cleaning up the tet mesh before writing
std::vector<std::array<int, 4>> tets;
std::vector<std::array<int, 4>> tetNeighbours;

void buildTetNeighbours()
{
    std::cerr << "Building tet neigbors" << std::endl;
    using namespace std;

    const int NT = tets.size();

    /** adjacent tets **/
    vector<vector<int>> adjacentTets(nv);

    int i = 0;
    for (const auto &t : tets)
    {
        for (int j = 0; j < 4; ++j)
            adjacentTets[t[j]].push_back(i);

        ++i;
    }

    for (auto &x : adjacentTets)
        sort(x.begin(), x.end());

    /** neighbor tets **/
    tetNeighbours.resize(NT);
    vector<int> candidates;

    for (int i = 0; i < NT; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            // Find Tet opposite to vertex j

            triangle cft({tets[i][(j + 1) % 4], tets[i][(j + 2) % 4], tets[i][(j + 3) % 4]});

            candidates.clear();

            int bc = cft[0];
            int mcs = adjacentTets[bc].size();
            if (adjacentTets[cft[1]].size() < mcs)
            {
                mcs = adjacentTets[cft[1]].size();
                bc = cft[1];
            }
            if (adjacentTets[cft[2]].size() < mcs)
                bc = cft[2];

            int curr = -1;
            for (auto at : adjacentTets[bc])
            {
                if (at != i)
                {
                    for (int k = 0; k < 4; k++)
                    {
                        triangle aft({tets[at][k], tets[at][(k + 1) % 4], tets[at][(k + 2) % 4]});
                        if (aft == cft)
                        {
                            curr = at;
                            break;
                        }
                        if (curr >= 0)
                            break;
                    }
                }
            }

            //            assert(curr != -1);
            tetNeighbours[i][j] = curr;
        }
    }
    std::cerr << "Done! " << std::endl;
}

int peel_by_winding_number(double thr)
{
    std::cerr << "Peeling... " << std::flush;
    Eigen::MatrixXd BC(T.number_of_finite_cells(), 3);

    Finite_cells_iterator cit = T.finite_cells_begin();
    for (int ti = 0; cit != T.finite_cells_end(); ++cit, ++ti)
    {
        Point c = CGAL::centroid(cit->vertex(0)->point().point(),
                                 cit->vertex(1)->point().point(),
                                 cit->vertex(2)->point().point(),
                                 cit->vertex(3)->point().point());
        BC(ti, 0) = c[0];
        BC(ti, 1) = c[1];
        BC(ti, 2) = c[2];
    }

    Eigen::VectorXd WN;
    Eigen::VectorXi CO;
    int nf = 0;
    for (int i = 0; i < np; i++)
        nf += P[i].tri.size();
    Eigen::MatrixXi F(nf, 3);
    nf = 0;
    for (int i = 0; i < np; i++)
        for (int j = 0; j < P[i].tri.size(); j++)
            F.row(nf++) << P[i].tri[j][0], P[i].tri[j][1], P[i].tri[j][2];
    igl::bfs_orient(F, F, CO);
    Eigen::MatrixXd EV(nv, 3);
    for (int i = 0; i < nv; i++)
        for (int j = 0; j < 3; j++)
            EV(i, j) = V[i][j];
    igl::winding_number(EV, F, BC, WN);

    int npt = 0;
    cit = T.finite_cells_begin();
    for (int ti = 0; cit != T.finite_cells_end(); ++cit, ++ti)
        if (fabs(WN(ti)) > thr)
        {
            cit->info() = 0;
        }
        else
        {
            cit->info() = -1;
            npt++;
        }
    std::cerr << "done - " << npt << " tets removed" << std::endl;
    return npt;
}

void write_cgal(std::string filename)
{
    std::cerr << "Wrtiting tets to " << filename << std::endl;

    CGAL::Unique_hash_map<VH, int> VH;
    Finite_vertices_iterator vit = T.finite_vertices_begin();
    for (int i = 0; vit != T.finite_vertices_end(); ++vit, ++i)
        VH[vit] = i;

    tets.resize(0);

    Finite_cells_iterator cit = T.finite_cells_begin();
    for (; cit != T.finite_cells_end(); ++cit)
        if (cit->info() != -1)
            tets.push_back({VH[cit->vertex(0)],
                            VH[cit->vertex(1)],
                            VH[cit->vertex(2)],
                            VH[cit->vertex(3)]});

    nv = T.number_of_vertices();
    buildTetNeighbours();

    constexpr std::array<std::array<unsigned int, 3>, 4> faceIds{
        1, 2, 3,
        0, 3, 2,
        0, 1, 3,
        0, 2, 1};

    for (auto &t : tets)
        for (int j = 0; j < 4; j++)
            t[j]++;
    nv++;

    int ti = 0;
    for (auto &nbh : tetNeighbours)
    {
        for (int j = 0; j < 4; ++j)
        {
            if (nbh[j] == -1)
                tets.push_back({tets[ti][faceIds[j][0]],
                                tets[ti][faceIds[j][1]],
                                tets[ti][faceIds[j][2]],
                                0});
        }
        ++ti;
    }

    buildTetNeighbours();

    std::ofstream cof(filename);
    cof << "3\n";

    cof << T.number_of_vertices() << "\n";

    cof << std::setprecision(std::numeric_limits<double>::max_digits10);

    for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
        cof << vit->point()[0]
            << " " << vit->point()[1]
            << " " << vit->point()[2]
            << " " << vit->point().weight()
            << "\n";

    cof << tets.size() << "\n";
    for (auto t : tets)
        cof << t[0] << " " << t[1] << " " << t[2] << " " << t[3] << "\n";

    for (auto n : tetNeighbours)
        cof << n[0] << " " << n[1] << " " << n[2] << " " << n[3] << "\n";

    cof.close();
}

void write_cgal_ordered(std::string filename)
{
    std::cerr << "Wrtiting tets to " << filename << std::endl;

    std::ofstream cof(filename);
    cof << "3\n";

    cof << nv << "\n";

    cof << std::setprecision(std::numeric_limits<double>::max_digits10);

    Finite_vertices_iterator vit = T.finite_vertices_begin();
    for (int i = 0; vit != T.finite_vertices_end(); ++vit, i++)
    {
        cof << V[i][0] << " " << V[i][1] << " " << V[i][2] << " " << W[i] << "\n";
        vit->info()++;
    }
    T.infinite_vertex()->info() = 0;

    cof << T.number_of_cells() << "\n";

    Rt::All_cells_iterator cit = T.all_cells_begin();
    for (int i = 0; cit != T.all_cells_end(); ++cit, i++)
    {
        cof << cit->vertex(0)->info() << " "
            << cit->vertex(1)->info() << " "
            << cit->vertex(2)->info() << " "
            << cit->vertex(3)->info() << "\n";
        cit->info() = i;
    }

    for (cit = T.all_cells_begin(); cit != T.all_cells_end(); ++cit)
        cof << cit->neighbor(0)->info() << " "
            << cit->neighbor(1)->info() << " "
            << cit->neighbor(2)->info() << " "
            << cit->neighbor(3)->info() << "\n";

    cof.close();
}
    


void write_woff(std::string filename)
{
    std::ofstream vf(filename);
    vf << "WOFF" << std::endl;
    vf << nv << " " << np << " 0" << std::endl;
    vf << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (int i = 0; i < nv; i++)
    {
        vf << V[i][0] << " " << V[i][1] << " " << V[i][2] << " " << W[i] << std::endl;
    }
    for (int i = 0; i < np; i++)
    {
        vf << P[i].be.size();
        for (auto e : P[i].be)
            vf << " " << e[0];
        vf << std::endl;
    }
    vf.close();
}


int ntri()
{
    int nt = 0;
    for (int i = 0; i < np; i++)
        nt += P[i].tri.size();
    return nt;
}


void write_off(std::string filename)
{
    std::ofstream vf(filename);
    vf << "OFF" << std::endl;
    vf << nv << " " << ntri() << " 0\n";
    for (int i = 0; i < nv; i++)
    {
        vf << V[i][0] << " " << V[i][1] << " " << V[i][2] << "\n";
    }
    for (int i = 0; i < np; i++)
        for (int ii = 0; ii < P[i].tri.size(); ii++)
        {
            int v0 = P[i].tri[ii][0];
            int v1 = P[i].tri[ii][1];
            int v2 = P[i].tri[ii][2];

            vf << "3 " << v0 << " " << v1 << " " << v2 << "\n";
        }
    vf.close();
}

void write_fc_off(std::string filename, std::vector<int> r)
{
    std::ofstream vf(filename);
    vf << "COFF" << std::endl;
    vf << nv << " " << ntri() << " 0" << std::endl;
    for (int i = 0; i < nv; i++)
    {
        vf << V[i][0] << " " << V[i][1] << " " << V[i][2] << std::endl;
    }
    for (int i = 0; i < np; i++)
        for (int ii = 0; ii < P[i].tri.size(); ii++)
        {
            int v0 = P[i].tri[ii][0];
            int v1 = P[i].tri[ii][1];
            int v2 = P[i].tri[ii][2];

            vf << "3 " << v0 << " " << v1 << " " << v2;
            switch (r[i])
            {
                case 0 : vf << " 0 255 0\n"; break;
                case 1 : vf << " 255 255 0\n"; break;
                case 2 : vf << " 255 0 0\n"; break;
                default : vf << " 0 0 255\n"; break;
            }
        }
    vf.close();
}

void mark_domains(CDT &ct,
                  CDT::Face_handle start,
                  int index,
                  std::list<CDT::Edge> &border)
{
    if (start->info().nesting_level != -1)
    {
        return;
    }
    std::list<CDT::Face_handle> queue;
    queue.push_back(start);
    while (!queue.empty())
    {
        CDT::Face_handle fh = queue.front();
        queue.pop_front();
        if (fh->info().nesting_level == -1)
        {
            fh->info().nesting_level = index;
            for (int i = 0; i < 3; i++)
            {
                CDT::Edge e(fh, i);
                CDT::Face_handle n = fh->neighbor(i);
                if (n->info().nesting_level == -1)
                {
                    if (ct.is_constrained(e))
                        border.push_back(e);
                    else
                        queue.push_back(n);
                }
            }
        }
    }
}
void mark_domains(CDT &cdt)
{
    for (CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it)
    {
        it->info().nesting_level = -1;
    }
    std::list<CDT::Edge> border;
    mark_domains(cdt, cdt.infinite_face(), 0, border);
    while (!border.empty())
    {
        CDT::Edge e = border.front();
        border.pop_front();
        CDT::Face_handle n = e.first->neighbor(e.second);
        if (n->info().nesting_level == -1)
        {
            mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
        }
    }
}



int triangulate_poly(int pid)
{
    poly &p = P[pid];
    int nf = p.be.size() + 2 * p.iv.size() - 2;
    p.tri.resize(nf);
    if (nf == 1)
    {
        auto eit = p.be.begin();
        for (int v = 0; v < 3; v++, ++eit)
        {
            p.tri[0][v] = (*eit)[0];
        }
        return 1;
    }
    CDT cdt;
    std::vector<std::pair<CDT::Point, int>> points;
    for (auto e : p.be)
    {
         points.push_back(std::make_pair(P[pid].plane.to_2d(V[e[0]]), e[0]));
    }
    if (!p.iv.empty())
        for (auto v : p.iv)
            points.push_back(std::make_pair(P[pid].plane.to_2d(V[v]), v));

    cdt.insert(points.begin(), points.end());
    std::unordered_map<int, CDT::Vertex_handle> vhm;
    CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin();
    for (; vit != cdt.finite_vertices_end(); ++vit)
        vhm[vit->info()] = vit;
    for (auto e : p.be)
        cdt.insert_constraint(vhm[e[0]], vhm[e[1]]);
    mark_domains(cdt);
    CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
    int fc = 0;
    for (; fit != cdt.finite_faces_end(); ++fit)
        if (fit->info().in_domain())
        {
            for (int v = 0; v < 3; v++)
            {
                p.tri[fc][v] = fit->vertex(v)->info();
            }
            fc++;

        }
    return fc;
}

int triangulate_cpoly_w(int pid)
{
    poly &p = P[pid];
    int nf = p.be.size() + 2 * p.iv.size() - 2;
    p.tri.resize(nf);
    if (nf == 1) // it's a triangle, so just create the single triangle from the edges
    {
        auto eit = p.be.begin();
        for (int v = 0; v < 3; v++, ++eit)
        {
            p.tri[0][v] = (*eit)[0];
        }
        return 1;
    }
    Rt2 rt2;
    //CDT cdt;
    std::vector<std::pair<Rt2::Point, int>> points;
    for (auto e : p.be)
    {
        points.push_back(std::make_pair(Rt2::Point(P[pid].plane.to_2d(V[e[0]]),W[e[0]])
                                        , e[0]));
    }
    if (!p.iv.empty())
        for (auto v : p.iv)
            points.push_back(std::make_pair(Rt2::Point(P[pid].plane.to_2d(V[v]),W[v]), v));
    rt2.insert(points.begin(), points.end());
    std::unordered_map<int, Rt2::Vertex_handle> vhm;
    
    Rt2::Finite_vertices_iterator vit = rt2.finite_vertices_begin();
    for (; vit != rt2.finite_vertices_end(); ++vit)
        vhm[vit->info()] = vit;

//    for (auto e : p.be)
//        cdt.insert_constraint(vhm[e[0]], vhm[e[1]]);
// instead probalby check that constraints are there
    
    //mark_domains_rt(rt2);
    Rt2::Finite_faces_iterator fit = rt2.finite_faces_begin();
    int fc = 0;
    for (; fit != rt2.finite_faces_end(); ++fit)
        //if (fit->info().in_domain())
        {
            for (int v = 0; v < 3; v++)
            {
                p.tri[fc][v] = fit->vertex(v)->info();
            }
            fc++;
        }
    if (fc != nf)
    {
        std::cerr << "something very wrong?" << fc << " = " << nf << std::endl;
        
        exit(1);
    }
    return fc;
}


int triangulate()
{
    int nf = 0;
    for (int i = 0; i < np; i++)
        nf += triangulate_poly(i);
    return nf;
}

void compute_planes()
{
    for (int i = 0; i < np; i++)
    {
        std::vector<Point> v;
        for (auto e : P[i].be)
            v.push_back(V[e[0]]);
        if (v.size() == 3)
            P[i].plane = K::Plane_3(v[0], v[1], v[2]);
        else
            linear_least_squares_fitting_3(v.begin(), v.end(), P[i].plane, CGAL::Dimension_tag<0>());
    }
}

// reads off and woff
bool read_off(std::string filename)
{
    std::ifstream pf(filename);
    std::string line;
    std::getline(pf, line);
    std::istringstream isf(line);

    std::string ft;
    isf >> ft;

//    std::cerr << "File type: X" << ft << "X" << std::endl;
    if (ft.size() < 3 || (ft.substr(ft.size() - 3)) != "OFF")
        return false;
    bool woff = (ft == "WOFF");

    int ne;
    std::getline(pf, line);
    std::istringstream iss(line);
    iss >> nv >> np >> ne;

    std::cerr << "Reading " << nv << " vertices ... " << std::flush;

    V.resize(nv);
    W.resize(nv);

    for (int i = 0; i < nv; i++)
    {
        std::getline(pf, line);
        std::istringstream iss(line);
        double x, y, z, w = 0.0;
        if (woff)
            iss >> x >> y >> z >> w;
        else
            iss >> x >> y >> z;
        V[i] = Point(x, y, z);
        W[i] = w;
    }
 
    std::cerr << "done\nReading " << np << " constrained faces ... " << std::flush;

    P = std::vector<poly>(np);

    double mina = 180.0;
    int nf = 0;
    for (int i = 0; i < np; i++)
    {
        std::getline(pf, line);
        std::istringstream iss(line);
        int d;
        iss >> d;
        if (d < 3)
            std::cerr << "Degree < 3" << std::endl;
        int id[d];
        for (int v = 0; v < d; v++)
            iss >> id[v];
        for (int v = 0; v < d; v++)
        {
            //edge e = edge{id[v],id[(v+1)%d]};
            edge e{id[v], id[(v + 1) % d]};
            P[i].be.push_back(e);
            E2P[e].push_back(i);
        }
        for (int v = 0; v < d; v++)
        {
            double a = CGAL::approximate_angle(V[id[v]],V[id[(v+1)%d]],V[id[(v+2)%d]]);
            if (a < mina) mina = a;
        }
    }

    std::cerr << "done!" << std::endl;
    std::cerr << "Smallest interior angle is " << mina << std::endl;
    
    std::vector<int> degree(nv, 0);
    for (auto ep : E2P)
    {
        degree[ep.first[0]]++;
        degree[ep.first[1]]++;
    }
    for (int i = 0; i < nv; i++)
        if (degree[i] < 3)
            std::cerr << "Degree of vertex " << i << " is " << degree[i] << std::endl;
    return true;
}

void read_igl(std::string filename)
{
    Eigen::MatrixXd EV;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh(filename, EV, F);
    nv = EV.rows();
    V.resize(nv);
    W.resize(nv, 0.0);
    for (int i = 0; i < nv; i++)
        V[i] = Point(EV(i, 0), EV(i, 1), EV(i, 2));

    np = F.rows();
    P = std::vector<poly>(np);
    for (int i = 0; i < np; i++)
    {
        int d = 3;
        for (int v = 0; v < d; v++)
        {
            edge e{F(i, v), F(i, (v + 1) % d)};
            P[i].be.push_back(e);
            E2P[e].push_back(i);
        }
    }
    std::cerr << "Read " << nv << " vertices and " << np << " triangular polygons" << std::endl;
    std::vector<bool> bb(np, false);
}

bool remove_edge(edge e)
{
    auto pit = E2P.find(e);
    if (pit == E2P.end())
        return false;
    std::vector<int> pid = pit->second;
    if (pid.size() != 2)
        return false;
    poly p0 = P[pid[0]];
    poly p1 = P[pid[1]];
    std::vector<int> nbv;
    return true;
}

bool split_edge(edge e)
{
    auto pit = E2P.find(e);
    if (pit == E2P.end())
        return false;

    Point em = CGAL::midpoint(V[e[0]], V[e[1]]);
    V.push_back(em);
    W.push_back(0.0);

    for (auto pid : pit->second)
    {
        auto eit = std::find(P[pid].be.begin(), P[pid].be.end(), e);
        if (eit == P[pid].be.end())
            return false;

        int v = eit->v[0];
        eit->v[0] = nv;
        P[pid].be.insert(eit, edge{v, nv});

        E2P[edge{e[0], nv}].push_back(pid);
        E2P[edge{e[1], nv}].push_back(pid);
    }
    E2P.erase(e);
    nv++;
    return true;
}

int split_missing_edges()
{
    std::cerr << "Checking edges... " << std::endl;
    std::vector<std::pair<edge, std::vector<int>>> me;
    for (auto it : E2P)
    {
        edge e = it.first;
        CH ch;
        int c0, c1;
        if (!ve[e[0]] || !ve[e[1]] ||
            !T.is_edge(id2vh[e[0]], id2vh[e[1]], ch, c0, c1))
            me.push_back(it);
    }
    int nme = me.size();
    std::cerr << "Edges missing (encroached): " << nme << std::endl;
    if (nme == 0)
        return 0;

    std::set<int> mp;
    for (int i = 0; i < nme; i++)
    {
        edge e = me[i].first;
        Point em = CGAL::midpoint(V[e[0]], V[e[1]]);
        V.push_back(em);
        W.push_back(0.0);

        for (auto pid : me[i].second)
        {
            auto eit = std::find(P[pid].be.begin(), P[pid].be.end(), e);
            if (eit == P[pid].be.end())
                std::cerr << "Edge " << e[0] << ", " << e[1] << " not found in poly " << pid << std::endl;

            int v = eit->v[0];
            eit->v[0] = nv;
            P[pid].be.insert(eit, edge{v, nv});

            E2P[edge{e[0], nv}].push_back(pid);
            E2P[edge{e[1], nv}].push_back(pid);
            mp.insert(pid);
        }
        E2P.erase(e);
        nv++;
    }
    std::cerr << "Modified " << mp.size() << std::endl;
    std::vector<bool> af(np, true);

    for (auto i : mp)
    {
        triangulate_poly(i);
        af[i] = false;
    }

    //write_fc_off("splits" + std::to_string(splits) + ".off", af);

    return nme;
}


void missing_faces(std::vector<int> &cmt)
{
    cmt.resize(0);
    for (int i = 0; i < np; i++)
    {
        for (int ii = 0; ii < P[i].tri.size(); ii++)
        {
            triangle t = P[i].tri[ii];
            CH ch;
            if (!(ve[t[0]] && ve[t[1]] && ve[t[2]]))
            {
                cmt.push_back(i);
                continue;
            }
            int c0, c1, c2;
            if (!T.is_facet(id2vh[t[0]], id2vh[t[1]], id2vh[t[2]],
                            ch, c0, c1, c2))
                cmt.push_back(i);
        }
    }
}

int split_missing_faces()
{

    std::unordered_set<int> mp;
    std::unordered_set<edge> ebs;
    std::vector<std::pair<Point, int>> ccp;
    int ni = 0;

    for (int i = 0; i < np; i++)
    {
        std::unordered_set<edge> eis;

        for (int ii = 0; ii < P[i].tri.size(); ii++)
        {
            triangle t = P[i].tri[ii];
            std::cerr << "Checking triangle " << t[0] << ", " << t[1] << ", " << t[2] << std::endl;
            // we can assume all vertices are there, because all edges are there
            CH ch;
            int c0, c1, c2;
            if (!T.is_facet(id2vh[t[0]], id2vh[t[1]], id2vh[t[2]],
                            ch, c0, c1, c2))
            {
                bool icc = true;
                Point cc = CGAL::circumcenter(V[t[0]], V[t[1]], V[t[2]]);
                std::cerr << "cc: " << cc << std::endl;
                for (int j = 0; j < 3; j++)
                {
                    // consider edge j, j+2
                    edge e = edge{t[j], t[(j + 2) % 3]};
                    // is it part of the boundary of P[i]?
                    auto eit = std::find(P[i].be.begin(), P[i].be.end(), e);
                    bool ob = (eit != P[i].be.end());
                    // is the angle opposite e obtuse?
                    bool wa = (CGAL::angle(V[t[j]], V[t[(j + 1) % 3]], V[t[(j + 2) % 3]]) != CGAL::ACUTE);
                    // if angle is obtuse split edge in any case
                    // if edge is on boundary and encroached by cc also split
                    if (wa)
                    {
                        if (ob)
                            ebs.insert(e);
                        else
                        {
                            std::cerr << "Should insert on interior edge " << e[0] << " - " << e[1] << " in poly " << i << std::endl;
                            auto ieit = std::find(eis.begin(), eis.end(), e);
                            if (ieit == eis.end())
                            {
                                std::cerr << "Not found for polygon " << i << "  -> Inserting" << std::endl;
                                Point em = CGAL::midpoint(V[e[0]], V[e[1]]);
                                V.push_back(em);
                                W.push_back(0.0);
                                P[i].iv.push_back(nv);
                                mp.insert(i);
                                nv++;
                                ni++;
                                eis.insert(e);
                            }
                        }
                        icc = false;
                        break;
                    }
                    if (ob)
                    {
                        Point em = CGAL::midpoint(V[e[0]], V[e[1]]);
                        if ((cc - em).squared_length() <=
                            (V[e[0]] - em).squared_length())
                        {
                            ebs.insert(e);
                            icc = false;
                            break;
                        }
                    }
                }
                if (icc)
                {
                    std::cerr << "Potentially inserting cc" << std::endl;
                    ccp.push_back(std::make_pair(cc, i));
                }
                /*
                    // insert cc only if it doesn't encroach an edge of P[i]
                bool icc = true;    
                for (auto e : P[i].be)
                {
                    Point em = CGAL::midpoint(V[e[0]],V[e[1]]);
                    if ((cc-em).squared_length() <=
                        (V[e[0]]-em).squared_length())
                    {
                        ees.insert(e);
                        icc = false;
                        std::cerr << "Encroaching!" << std::endl;
                        break;
                    }
                }
                if (icc)
                    ccp.push_back(std::make_pair(cc,i));
                 */
            }
        }
    }
    std::cerr << "Interior insertions: " << ni << std::endl;
    std::cerr << "Encroached edges by ccs: " << (int)ebs.size() << std::endl;

    for (auto e : ebs)
    {
        Point em = CGAL::midpoint(V[e[0]], V[e[1]]);
        V.push_back(em);
        W.push_back(0.0);

        for (auto pid : E2P[e])
        {
            auto eit = std::find(P[pid].be.begin(), P[pid].be.end(), e);
            if (eit == P[pid].be.end())
                std::cerr << "Edge " << e[0] << ", " << e[1] << " not found in poly " << pid << std::endl;

            int v = eit->v[0];
            eit->v[0] = nv;
            P[pid].be.insert(eit, edge{v, nv});

            E2P[edge{e[0], nv}].push_back(pid);
            E2P[edge{e[1], nv}].push_back(pid);
            mp.insert(pid);
        }
        E2P.erase(e);
        nv++;
    }
    std::cerr << "Modified " << mp.size() << "by edge splitting" << std::endl;
    std::vector<bool> af(np, true);

    if (mp.size() == 0)
    {
        for (auto cp : ccp)
        {
            V.push_back(cp.first);
            W.push_back(0.0);
            P[cp.second].iv.push_back(nv);
            mp.insert(cp.second);
            //std::cerr << "Interior point " << nv << " in poly " << cp.second << std::endl;
            nv++;
        }
    }

    for (auto i : mp)
    {
        triangulate_poly(i);
        af[i] = false;
    }
//    write_fc_off("splits" + std::to_string(splits) + ".off", af);
//    std::cerr << "wrote splits" + std::to_string(splits) + ".off" << std::endl;
    return ni + ebs.size() + ccp.size();
}

int solve(int nr, double mintol, double maxtol, double height_factor, double weight_factor)
{
    int nic = 0;

    id2vh.resize(0);
    id2vh.resize(nv);

    std::vector<std::pair<Weighted_point, int>> Pi(nv);

    int nf = ntri();
    std::cerr << "Number of constrained triangles: " << nf << std::endl;

    std::vector<bool> dve(nf, false), ignore(nf, false);
    std::vector<std::unordered_map<int, double>> constr_tol(nf);
    std::vector<std::unordered_set<int>> constr(nf);

    Eigen::VectorXd H;
    Eigen::MatrixXd G;
    if (!gabriel)
    {
        H = Eigen::VectorXd::Zero(nf);
        G = Eigen::MatrixXd::Zero(nv, 3);
    }

    try
    {
        GRBModel model = GRBModel(*env);
        model.set(GRB_IntParam_OutputFlag, 0);

        // Create variables
        //GRBVar* w = model.addVars(&lb[0], NULL, NULL, NULL, NULL, nv);
        GRBVar *w = model.addVars(NULL, NULL, NULL, NULL, NULL, nv);

        for (int i = 0; i < nv; i++)
            if (W[i] > 0.0)
                w[i].set(GRB_DoubleAttr_Start, W[i]);

        GRBVar *h, *g;
        if (!gabriel)
        {
            std::vector<double> lbf(nf, -GRB_INFINITY);
            std::vector<double> lbv(3 * nv, -GRB_INFINITY);
            h = model.addVars(&lbf[0], NULL, NULL, NULL, NULL, nf);
            g = model.addVars(&lbv[0], NULL, NULL, NULL, NULL, 3 * nv);
        }

        // Set objective
        GRBQuadExpr qe;
        GRBLinExpr le;
        std::vector<double> onev(nv, weight_factor); //,b(nv);
        std::vector<double> onep(nf, height_factor); //,b(nv);
        // squared norm of height vector
        if (!gabriel)
            qe.addTerms(&onep[0], h, h, nf);
        // squared norm of weight vector
        // qe.addTerms(&one[0],w,w,nv);
        qe.addTerms(&onev[0], w, w, nv);
        le.addTerms(&onev[0], w, nv);

        if (minimize_heights)
            model.setObjective(qe, GRB_MINIMIZE);
        else
            model.setObjective(le, GRB_MINIMIZE);
        int ntc = 0;
        int nutc = 0;

        Eigen::MatrixXd Z(nf,3);
        for (int r = 0; r < nr; r++)
        {
            std::cout << "Round " << r << std::endl;
 
           // for (int i = 0; i < np; i++)
           //     triangulate_cpoly_w(i);
            
            for (int i = 0; i < nv; i++)
                Pi[i] = std::make_pair(Weighted_point(V[i], W[i]), i);

            T.clear();
            T.insert(Pi.begin(), Pi.end());
            //T.is_valid(true);

            // set up integer ids for vertices and mark vertices in the triangulation
            ve.resize(0);
            ve.resize(nv, false);
            Finite_vertices_iterator vit = T.finite_vertices_begin();
            for (; vit != T.finite_vertices_end(); ++vit)
            {
                id2vh[vit->info()] = vit;
                ve[vit->info()] = true;
            }

            int nacv = 0;
            int nmv = 0;
            int nmt = 0;
            int ti = 0;

            nmv = nv - T.number_of_vertices();
            std::cout << "Number of missing vertices: " << nmv << " / " << nv << std::endl;
            if (nmv < 0)
            {
                for (int i = 0; i < nv; i++)
                    if (!ve[i])
                    {
                        Eigen::RowVector3d p(V[i][0],
                                             V[i][1],
                                             V[i][2]);

                        //double xs = p0.squaredNorm();
                        Point z(G(i, 0) + V[i][0], G(i, 1) + V[i][1], G(i, 2) + V[i][2]);
                        VH vh = T.nearest_power_vertex(z);
                        int vi = vh->info();
                        if (vi == i)
                            std::cerr << "Weird!" << std::endl;
                        std::cerr << "Adding " << i << " - " << vi << std::endl;

                        Eigen::RowVector3d pi(V[vi][0],
                                              V[vi][1],
                                              V[vi][2]);

                        Eigen::Vector3d a = 2.0 * (pi - p);
                        double rhs = pi.squaredNorm() - p.squaredNorm() - p * a - 1e-6;

                        GRBLinExpr l = 0;
                        if (!gabriel)
                        {
                            l += a[0] * g[i];
                            l += a[1] * g[i + nv];
                            l += a[2] * g[i + nv + nv];
                        }
                        l -= w[i];
                        l += w[vi];
                        model.addConstr(l, GRB_LESS_EQUAL, rhs, std::to_string(-i));
                         ntc++;
                    }
            }

            std::vector<int> ac(np,0);
            for (int i = 0; i < np; i++)
                for (int ii = 0; ii < P[i].tri.size(); ii++, ti++)
                    if (!ignore[ti])
                    {
                        int v0 = P[i].tri[ii][0];
                        int v1 = P[i].tri[ii][1];
                        int v2 = P[i].tri[ii][2];

                        // check if triangle v0,v1,v2 is part of the triangulation
                        // if so, continue with next triangle
                        // if (!(ve[v0] &&ve[v1] && ve[v2]))
                        //     continue;

                        CH ch;
                        int c0, c1, c2, c3 = -1;
                        if (ve[v0] && ve[v1] && ve[v2] &&
                            T.is_facet(id2vh[v0], id2vh[v1], id2vh[v2],
                                       ch, c0, c1, c2))
                        {
                            c3 = 6 - c0 - c1 - c2;
                            if (!gabriel && !minimize_heights)
                                continue;
                            if (T.is_Gabriel(ch, c3))
                                continue;
                        }

                        Eigen::RowVector3d p0(V[v0][0],
                                              V[v0][1],
                                              V[v0][2]);
                        Eigen::RowVector3d p1(V[v1][0],
                                              V[v1][1],
                                              V[v1][2]);
                        Eigen::RowVector3d p2(V[v2][0],
                                              V[v2][1],
                                              V[v2][2]);

                        Eigen::Vector3d n(P[i].plane.orthogonal_vector()[0],
                                          P[i].plane.orthogonal_vector()[1],
                                          P[i].plane.orthogonal_vector()[2]);

                        double x0s = p0.squaredNorm();
                        double x1s = p1.squaredNorm();
                        double x2s = p2.squaredNorm();

                        Eigen::Matrix3d M;
                        M.row(0) = 2.0 * (p1 - p0);
                        M.row(1) = 2.0 * (p2 - p0);
                        M.row(2) = n;

                        Eigen::Matrix3d Mi = M.inverse();

                        Eigen::Vector3d r(x1s - x0s,
                                          x2s - x0s,
                                          p0 * n);

                        Eigen::Vector3d m = Mi * r;

                        int vi;
                        double tol = mintol;
                        if (c3 < 0 || gabriel)
                        {
                            Eigen::Vector3d b(W[v0] - W[v1],
                                              W[v0] - W[v2],
                                              0);
                            if (!gabriel)
                                b[2] = H[ti];

                            Z.row(ti) = m + Mi * b;
                            Point p(Z(ti,0), Z(ti,1), Z(ti,2));
                            //std::cerr << p << std::endl;
                            // because of what's going on behind the scenes in CGAL
                            // the following is highly inefficient and should be optimized
 
                            vi = T.nearest_power_vertex(p)->info();

                            auto vit = constr_tol[ti].find(vi);
                            if (vit != constr_tol[ti].end())
                            {
                                tol = vit->second;
                                tol *= 2.0;
                                if (tol > maxtol)
                                {
                                    ignore[ti] = true;
                                    nic++;
                                    std::cerr << "Ignoring " << ti << " because it has reach max tolerance" << std::endl;
                                    GRBConstr *c = model.getConstrs();
                                    for (int ci = 0; ci < model.get(GRB_IntAttr_NumConstrs); ++ci)
                                    {
                                        int ctid = atoi(c[ci].get(GRB_StringAttr_ConstrName).c_str());
                                        if (ctid == ti)
                                            model.remove(c[ci]);
                                    }
                                    model.update();
                                    continue;
                                }
                                vit->second = tol;
                            }
                            else
                                constr_tol[ti][vi] = mintol;
                        }
                        else
                        {
                            VH vh = ch->vertex(c3);
                            if (T.is_infinite(vh) || constr[ti].count(vh->info()) > 0)
                                vh = T.mirror_vertex(ch, c3);
                            if (T.is_infinite(vh) || constr[ti].count(vh->info()) > 0)
                                continue;
                            vi = vh->info();
                            if (constr_tol[ti].find(vi) != constr_tol[ti].end())
                                continue;
                            constr_tol[ti][vi] = mintol;
                        }
 
                        Eigen::RowVector3d pi(V[vi][0],
                                              V[vi][1],
                                              V[vi][2]);

                        Eigen::RowVector3d a = 2.0 * (pi - p0);
                        double rhs = -a * m + pi.squaredNorm() - x0s - tol;

                        a *= Mi;
                        GRBLinExpr l = 0;
                        l += (a[0] + a[1] - 1.0) * w[v0];
                        l -= a[0] * w[v1];
                        l -= a[1] * w[v2];
                        if (!gabriel)
                            l += a[2] * h[ti];
                        l += w[vi];
                        model.addConstr(l, GRB_LESS_EQUAL, rhs, std::to_string(ti));
                        constr[ti].insert(vi);
                        ntc++;
                        //}
                        nmt++;
                        ac[i] = 1;
                        if (c3 < 0) ac[i] = 2;
                    }

            if (nmt == 0)
            {
                std::cout << "No constraint violations anymore, breaking in round: " << r << std::endl;
                int mnc = 0, mnci = -1, tncv = 0;
                std::unordered_map<int,int> chist;
                for (int csi = 0; csi < ti; csi++)
                {
                    int cts = constr_tol[csi].size();
                    if (cts > mnc)
                    {
                        mnci = csi;
                        mnc = constr_tol[csi].size();
                    }
                    tncv += cts;
                    chist[cts]++;
                }
                std::cerr << "Constraint triangle " << mnci << " has max # of constraints: " << mnc << std::endl;
                std::cerr << "Total constraints: " << tncv << ", per constrained triangle: " << ((double)tncv / (double)ntri()) << std::endl;
                std::cerr << "Histogram" << std::endl;
                for ( int csi = 0; csi <= mnc; csi++)
                    std::cerr << csi << ": " << chist[csi] << std::endl;
                return nic;
            }
            else
                std::cout << "Number of missing triangles: " << nmt << " / " << ti << std::endl
                          << "Total number of constraints added: " << ntc << std::endl;

            model.optimize();

            int status = model.get(GRB_IntAttr_Status);
            std::cerr << "Status: " << status << std::endl;

            if (status == GRB_INFEASIBLE || status == GRB_INF_OR_UNBD)
            {
                model.computeIIS();
                std::cerr << "\nThe following constraint cannot be satisfied: " << std::flush;
                GRBConstr *c = model.getConstrs();
                int tid = -1;
                for (int i = 0; i < model.get(GRB_IntAttr_NumConstrs); ++i)
                {
                    if (c[i].get(GRB_IntAttr_IISConstr) == 1)
                    {
                        tid = atoi(c[i].get(GRB_StringAttr_ConstrName).c_str());
                        std::cerr << tid << std::endl;
                        if (tid >= 0)
                            break;
                    }
                }
                if (ignore[tid])
                    std::cerr << "Already ignoring tid: " << tid << std::endl;
                ignore[tid] = true;
                nic++;
                for (int i = 0; i < model.get(GRB_IntAttr_NumConstrs); ++i)
                {
                    int ctid = atoi(c[i].get(GRB_StringAttr_ConstrName).c_str());
                    if (ctid == tid)
                        model.remove(c[i]);
                }
                model.update();
            }
            else
            {
                for (int i = 0; i < nv; i++)
                    W[i] = w[i].get(GRB_DoubleAttr_X);
                if (!gabriel)
                {
                    for (int i = 0; i < nf; i++)
                        H[i] = h[i].get(GRB_DoubleAttr_X);
                    for (int i = 0; i < nv; i++)
                    {
                        G(i, 0) = g[i].get(GRB_DoubleAttr_X);
                        G(i, 1) = g[i + nv].get(GRB_DoubleAttr_X);
                        G(i, 2) = g[i + nv + nv].get(GRB_DoubleAttr_X);
                    }
                }
            }
        }
    }
    catch (GRBException e)
    {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    }
    W.resize(0);
    W.resize(nv, 0.0);
    return -1;
}

    

int main(int argc, char *argv[])
{
    CLI::App app{"CRT3"};

    std::string infilename = "default";
    app.add_option("-m,--mesh", infilename, "Input mesh, libigl supported format")->required();

    int nr = 100;
    app.add_option("-r,--rounds", nr, "# rounds for LP solving");

    int ns = 100;
    app.add_option("-s,--splits", ns, "# splitting rounds");

    bool w = false;
    app.add_flag("-w,--write", w, "Write tet and woff mesh (off in case of splits");

    bool peel = false;
    app.add_flag("-p,--peel", peel, "Peel outer tets before writing");
    double thr = 0.5;
    app.add_option("--pt", thr, "Threshold for winding number when peeling");

    gabriel = false;
    app.add_flag("-g,--gabriel", gabriel, "Ensure edges have gabriel property");

    minimize_heights = true;
    double height_factor = 1.0;
    CLI::Option *mh = app.add_option("--hf", height_factor, "Factor for minimizing heights");
    double weight_factor = 1.0;
    CLI::Option *bh = app.add_option("--wf", weight_factor, "Factor for minimizing weights");

    double tolerance = 1e-9;
    app.add_option("-t,--tolerance", tolerance, "minimal distance to constraints");

    CLI11_PARSE(app, argc, argv);

    size_t lastindex = infilename.find_last_of(".");
    std::string rawname = infilename.substr(0, lastindex);
    std::string suffix = infilename.substr(lastindex + 1);
    std::transform(suffix.begin(), suffix.end(), suffix.begin(), ::toupper);
    if (suffix == "OFF" || suffix == "WOFF")
        read_off(infilename);
    else
        read_igl(infilename);

    compute_planes();
    triangulate();

    int onv = nv;
    int onf = ntri();

    env = new GRBEnv();
    env->start();

    splits = 0;
    for (; splits < ns; splits++)
    {
        std::cout << "Loop " << splits << std::endl;
        int nif = solve(nr, tolerance, 1e2*tolerance, height_factor, weight_factor);
        std::cout << "Number of ignored triangles: " << nif << std::endl;
        if (nif == 0)
            break;
        int nse = split_missing_edges();
        if (nse == 0)
            split_missing_faces();
        env->resetParams();
    }
    std::vector<int> mf;
    missing_faces(mf);
    if (mf.size() == 0)
    {
        std::cout << "---\nFound regular triangulation of constraints after " << splits << " splits" << std::endl;
    }

    int npt = 0;
    if (peel)
        npt = peel_by_winding_number(thr);

    std::cerr << "Vertices: " << nv << " ( " << onv << " + " << (nv-onv) << " )" << std::endl;
    std::cout << "Surface Triangles: " << ntri() << " ( " << onf << " + " << (ntri()-onf) << " )" << std::endl;
    std::cout << "Tets: " << T.number_of_finite_cells() << " - " << npt << " (peeled) = " << (T.number_of_finite_cells() - npt) << std::endl;


    if (w)
    {
        write_woff(rawname + ".woff");
        if (peel)
            write_cgal(rawname + ".rct");
        else
            write_cgal_ordered(rawname + ".rct");
        if (splits > 0)
            write_off(rawname + "_s.off");
    }

    return 0;
}
