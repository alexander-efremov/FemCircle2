#include <string>
#include <iostream>
#include "graphs.h"
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
    strcat(cmd, "\"D:\\work\\Files\\Tools\\graphviz\\bin\\neato.exe\" -Tpng ");
    strcat(cmd, dotFilepath);
    strcat(cmd, " > ");
    strcat(cmd, pngFilename);
    strcat(cmd, " 2>&1 ");
    system(cmd);
}

void print_graph(const char *filename, Graph g) {
    std::ofstream out;
    out.open(filename);
    auto t = g.m_property.get();
    write_graphviz(out, g, VertexPropertyWriter(t->n, t->m));
    out.close();
}

GraphProperty::GraphProperty(unsigned int n, unsigned int m) : n(n), m(m) {}
