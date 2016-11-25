#ifndef FEM_CIRCLE_GRAPH_UTILS_H
#define FEM_CIRCLE_GRAPH_UTILS_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <ostream>

using namespace boost;

struct VertexProperty {
    bool is_corner;
    bool is_border;
};

struct GraphProperty {
    GraphProperty(unsigned int n, unsigned int m);

    unsigned int nx;
    unsigned int ny;
};

typedef adjacency_list<vecS, vecS, undirectedS, VertexProperty, no_property, GraphProperty> Graph;

bool is_graph_connected(Graph g);

void generate_png(const char *dotFilepath, const char *pngFilename);

void print_graph(const char *filename, Graph g);

#endif // FEM_CIRCLE_GRAPH_UTILS_H