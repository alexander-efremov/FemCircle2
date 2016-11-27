#include <string>
#include <iostream>
#include "graphs.h"
#include "utils.h"
#include <boost/graph/connected_components.hpp>

using namespace std;

bool is_graph_connected(const GraphDouble &g) {
    std::vector<int> component(num_vertices(g));
    int component_count = connected_components(g, &component[0]);
    return component_count == 1;
}

bool is_graph_connected(const GraphInt &g) {
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

class VertexPropertyWriterDouble {
public:
    VertexPropertyWriterDouble(const GraphDouble &g) : _g(g) {}

    template<class T1>
    void operator()(std::ostream &out, const T1 &v) const {
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
    GraphDouble _g;
};

class VertexPropertyWriterInt {
public:
    VertexPropertyWriterInt(const GraphInt &g) : _g(g) {}

    template<class T1>
    void operator()(std::ostream &out, const T1 &v) const {
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
    GraphInt _g;
};

struct GraphPropertyWriter {
    void operator()(std::ostream &out) const {
        out << "node [shape=point, width=.15, height=.15  style=\"filled\", fillcolor=\"grey\"]" << std::endl;
        out << "edge [style=dashed, fillcolor=\"grey\"]" << std::endl;
    }
};

void print_graph(const char *filename, const GraphDouble &g) {
    std::ofstream out;
    out.open(filename);
    write_graphviz(out, g, VertexPropertyWriterDouble(g), default_writer(), GraphPropertyWriter());
    out.close();
}

void print_graph(const char *filename, const GraphInt &g) {
    std::ofstream out;
    out.open(filename);
    write_graphviz(out, g, VertexPropertyWriterInt(g), default_writer(), GraphPropertyWriter());
    out.close();
}

// nx = NX_1, ny = NY_1

GraphDouble* create_graph_as_grid(int nx, int ny, double defaultValue) {
    assert(nx > 0);
    assert(ny > 0);

    GraphDouble *g = new GraphDouble(nx * ny, GraphProperty(nx, ny));
    for (int i = 0; i < nx; ++i) {
        int stride = i * nx;
        for (int j = 0; j < ny - 1; ++j) {
            auto e = add_edge((unsigned long long int) (j + stride), (unsigned long long int) (j + stride + 1), *g);
            g->m_vertices[e.first.m_source].m_property.m_value = defaultValue;
        }
    }

    for (size_t j = 0; j < ny - 1; ++j) {
        size_t stride = j * ny;
        for (size_t i = 0; i < nx; i++)
            add_edge(i + stride, i + stride + nx, *g);
    }
    return g;
}

GraphInt* create_graph_as_grid(int nx, int ny, int defaultValue) {
    assert(nx > 0);
    assert(ny > 0);

    GraphInt *g = new GraphInt(nx * ny, GraphProperty(nx, ny));
    for (int i = 0; i < nx; ++i) {
        int stride = i * nx;
        for (int j = 0; j < ny - 1; ++j) {
            auto e = add_edge((unsigned long long int) (j + stride), (unsigned long long int) (j + stride + 1), *g);
            g->m_vertices[e.first.m_source].m_property.m_value = defaultValue;
        }
    }

    for (size_t j = 0; j < ny - 1; ++j) {
        size_t stride = j * ny;
        for (size_t i = 0; i < nx; i++)
            add_edge(i + stride, i + stride + nx, *g);
    }
    return g;
}

double calc_graph_sum(const GraphDouble &a, int ox_len, int oy_len, bool isAbs) {
    double res = 0;
    for (int i = 0; i < ox_len; i++)
    {
        for (int j = 0; j < oy_len; j++)
        {
            if (isAbs)
                res += fabs(a.m_vertices[i * oy_len + j].m_property.m_value);
            else
                res += a.m_vertices[i * oy_len + j].m_property.m_value;
        }
    }
    return res;
}

