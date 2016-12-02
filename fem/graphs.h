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
                  double v, double u, double tau, int r_lvl, double hx_min, double hy_min, double hx, double hy,
                  unsigned int ideal_square_size_nx, unsigned int ideal_square_size_ny,unsigned int current_tl) : nx(
            nx), ny(ny), a(a), b(b), c(c), d(d), v(v), u(u), tau(tau),
            r_lvl(r_lvl),
            hx_min(hx_min),
            hy_min(hy_min),
            hx(hx),
            hy(hy),
            ideal_square_size_nx(ideal_square_size_nx),
            ideal_square_size_ny(ideal_square_size_ny),
            current_tl(current_tl)
    {}

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
    double hx_min;
    double hy_min;
    double hx;
    double hy;
    unsigned int ideal_square_size_nx;
    unsigned int ideal_square_size_ny;
    unsigned int current_tl;
};

typedef adjacency_list<vecS, vecS, undirectedS, VertexPropertyDouble, no_property, GraphProperty> Graph;

Graph *create_graph_as_grid(
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
        double hx_min,
        double hy_min,
        double hx,
        double hy,
        unsigned int ideal_square_size_nx,
        unsigned int ideal_square_size_ny,
        unsigned int current_tl,
        double defaultValue = -1.);

bool is_graph_connected(const Graph &g);

void print_graph(const char *filename, const Graph &g);

double calc_graph_sum(const Graph &g, bool isAbs);

void generate_png(const char *dotFilepath, const char *pngFilename);

#endif // FEM_CIRCLE_GRAPH_UTILS_H