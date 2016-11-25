#ifndef FEM_CIRCLE_GRAPH_UTILS_H
#define FEM_CIRCLE_GRAPH_UTILS_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <ostream>

using namespace boost;

struct VertexProperty {
};

struct GraphProperty {
    GraphProperty(unsigned int n, unsigned int m);

    unsigned int n;
    unsigned int m;
};

typedef adjacency_list<vecS, vecS, undirectedS, VertexProperty, no_property, GraphProperty> Graph;


class VertexPropertyWriter {
public:
    VertexPropertyWriter(int nx_1, int ny_1) : _nx_1(nx_1), _ny_1(ny_1) {}

    template<class T>
    void operator()(std::ostream &out, const T &v) const {

        int x = v % _ny_1;
        int y = v / _nx_1;
        y = _ny_1-y;
        out << " [pos=\"" << x << "," << y << "!\"]";
    }

private:
    int _nx_1;
    int _ny_1;
};

bool is_graph_connected(Graph g);

void generate_png(const char *dotFilepath, const char *pngFilename);

void print_graph(const char *filename, Graph g);

#endif // FEM_CIRCLE_GRAPH_UTILS_H