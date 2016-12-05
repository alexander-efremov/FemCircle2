#ifndef FEM_CIRCLE_GRAPH_UTILS_H
#define FEM_CIRCLE_GRAPH_UTILS_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <ostream>

using namespace boost;

struct VertexPropertyDouble {
    double m_value = 1000;
    int m_i;
    int m_j;
    int m_level = 0;

    bool IsInnerPoint(int nx1, int ny1) {
        return !IsCornerPoint(nx1, ny1) && !IsBorderPoint(nx1, ny1);
    }

    bool IsCornerPoint(int nx1, int ny1) {
        return is_corner_node(m_i, m_j, nx1, ny1);
    }

    bool IsBorderPoint(int nx1, int ny1) {
        return is_border_node(m_i, m_j, nx1, ny1);
    }

private:
    bool is_corner_node(int i, int j, int nx, int ny) {
        return (i == 0 && j == 0) || (i == 0 && j == ny - 1) || (i == nx - 1 && j == 0) || (i == nx - 1 && j == nx - 1);
    }

    bool is_border_node(int i, int j, int nx, int ny) {
        if (is_corner_node(i, j, nx, ny)) return false;
        return (i > 0 && j == 0) || (i == 0 && j > 0) || (i == nx - 1 && j > 0) || (i > 0 && j == nx - 1);
    }
};

struct GraphProperty {
    GraphProperty(unsigned int nx, unsigned int ny, double a, double b, double c, double d,
                  double v, double u, double tau, int r_lvl, double hx_min, double hy_min, double hx, double hy,
                  unsigned int ideal_square_size_nx, unsigned int ideal_square_size_ny) : nx_1(
            nx), ny_1(ny), a(a), b(b), c(c), d(d), v(v), u(u), tau(tau),
                                                                                                                   r_lvl(r_lvl),
                                                                                                                   hx_min(hx_min),
                                                                                                                   hy_min(hy_min),
                                                                                                                   hx(hx),
                                                                                                                   hy(hy),
                                                                                                                   ideal_square_size_nx(
                                                                                                                           ideal_square_size_nx),
                                                                                                                   ideal_square_size_ny(
                                                                                                                           ideal_square_size_ny),
                                                                                                                   current_tl(0) {}

    unsigned int nx_1;
    unsigned int ny_1;
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

Graph *create_graph_as_grid(const Graph &graph, double defaultValue = -1);

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
        double defaultValue = -1.);

bool is_graph_connected(const Graph &g);

void print_graph(const char *filename, const Graph &g);

double calc_graph_sum(const Graph &g, bool isAbs);

void generate_png(const char *dotFilepath, const char *pngFilename);

#endif // FEM_CIRCLE_GRAPH_UTILS_H