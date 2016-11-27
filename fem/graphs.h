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

struct VertexPropertyInt {
    int m_value;
};

struct GraphProperty {
    GraphProperty(unsigned int nx, unsigned int ny) : nx(nx), ny(ny) {}
    unsigned int nx;
    unsigned int ny;
};

typedef adjacency_list<vecS, vecS, undirectedS, VertexPropertyDouble, no_property, GraphProperty> GraphDouble;
typedef adjacency_list<vecS, vecS, undirectedS, VertexPropertyInt, no_property, GraphProperty> GraphInt;

GraphDouble* create_graph_as_grid(int nx, int ny, double defaultValue = -1.);
GraphInt* create_graph_as_grid(int nx, int ny, int defaultValue = -1);
bool is_graph_connected(const GraphDouble &g);
bool is_graph_connected(const GraphInt &g);
void print_graph(const char *filename, const GraphDouble &g);
void print_graph(const char *filename, const GraphInt &g);
double calc_graph_sum(const GraphDouble &a, int ox_len, int oy_len, bool isAbs);
void generate_png(const char *dotFilepath, const char *pngFilename);

#endif // FEM_CIRCLE_GRAPH_UTILS_H