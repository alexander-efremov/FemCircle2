#ifndef FEM_CIRCLE_GRAPH_UTILS_H
#define FEM_CIRCLE_GRAPH_UTILS_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <ostream>

using namespace boost;

struct VertexPropertyDouble {
    double m_value;
    int m_level = 0;
};

struct GraphProperty {
    GraphProperty(unsigned int nx, unsigned int ny, double a, double b, double c, double d,
                  double v, double u, double tau, int r_lvl, double hx_smallest, double hy_smallest) : nx(nx), ny(ny), a(a), b(b), c(c), d(d), v(v), u(u), tau(tau),
                                                                                                       r_lvl(r_lvl), hx_smallest(hx_smallest), hy_smallest(hy_smallest) {}

    unsigned int nx;
    unsigned int ny;
    double a;
    double b;
    double c;
    double d;
    double v;
    double u;
    double tau;
    int r_lvl;
    double hx_smallest;
    double hy_smallest;
};

typedef adjacency_list<vecS, vecS, undirectedS, VertexPropertyDouble, no_property, GraphProperty> GraphDouble;

GraphDouble *create_graph_as_grid(
        unsigned int nx_1,
        unsigned int ny_1,
        double a,
        double b,
        double c,
        double d,
        double v,
        double u,
        double tau,
        int r_lev,
        double hx_smallest,
        double hy_smallest,
        double defaultValue = -1.);

bool is_graph_connected(const GraphDouble &g);

void print_graph(const char *filename, const GraphDouble &g);

double calc_graph_sum(const GraphDouble &g, bool isAbs);

void generate_png(const char *dotFilepath, const char *pngFilename);

#endif // FEM_CIRCLE_GRAPH_UTILS_H