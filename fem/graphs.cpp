#include <string>
#include <iostream>
#include "graphs.h"
#include "utils.h"
#include <boost/graph/connected_components.hpp>

using namespace std;

bool is_graph_connected(Graph g) {
    std::vector<int> component(num_vertices(g));
    int component_count = connected_components(g, &component[0]);
    return component_count == 1;
}

void generate_png(const char *dotFilepath, const char *pngFilename) {
    char cmd[256] = {0};
    // добавить функцию чтобы можно было путь этот задавать
    strcat(cmd, "neato -Tpng ");
    strcat(cmd, dotFilepath);
    strcat(cmd, " > ");
    strcat(cmd, pngFilename);
    strcat(cmd, " 2>&1 ");
    system(cmd);
}

class VertexPropertyWriter {
public:
    VertexPropertyWriter(const Graph &g) : _g(g) {}

    template<class T>
    void operator()(std::ostream &out, const T &v) const {
        auto t = _g.m_property.get();
        int x = v % t->nx;
        int y = v / t->ny;
        out << " ["
            << " shape=point"
            << " xlabel=\"(" << x << "," << y << ")\""
            << "pos=\"" << x << "," << t->ny - y << "!\""
            << (is_corner_node(x, y, t->nx, t->ny) ? " fillcolor=\"red\" style=\"filled\"" : "")
            << (is_border_node(x, y, t->nx, t->ny) ? " fillcolor=\"blue\" style=\"filled\"" : "")
            << "]";
    }

private:
    Graph _g;
};

struct GraphPropertyWriter {
    void operator()(std::ostream &out) const {
        out << "node [shape=point, width=.15, height=.15  style=\"filled\", fillcolor=\"grey\"]" << std::endl;
        out << "edge [style=dashed, fillcolor=\"grey\"]" << std::endl;
    }
};

void print_graph(const char *filename, Graph g) {
    std::ofstream out;
    out.open(filename);
    write_graphviz(out, g, VertexPropertyWriter(g), default_writer(), GraphPropertyWriter());
    out.close();
}

GraphProperty::GraphProperty(unsigned int n, unsigned int m) : nx(n), ny(m) {}
