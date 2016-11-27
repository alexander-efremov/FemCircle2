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
// nx = NX_1, ny = NY_1

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
        int r_lvl,
        double hx_smallest,
        double hy_smallest,
        double defaultValue) {
    assert(nx_1 > 0);
    assert(ny_1 > 0);

    GraphDouble *g = new GraphDouble(nx_1 * ny_1,
                                     GraphProperty(nx_1, ny_1, a, b, c, d, u, v, tau, r_lvl, hx_smallest, hy_smallest));
    for (int i = 0; i < nx_1; ++i) {
        int stride = i * nx_1;
        for (int j = 0; j < ny_1 - 1; ++j) {
            auto e = add_edge((unsigned long long int) (j + stride), (unsigned long long int) (j + stride + 1), *g);
            g->m_vertices[e.first.m_source].m_property.m_value = defaultValue;
        }
    }

    for (int j = 0; j < ny_1 - 1; ++j) {
        int stride = j * ny_1;
        for (int i = 0; i < nx_1; i++) {
            add_edge((unsigned long long int) (i + stride), (unsigned long long int) (i + stride + nx_1), *g);
        }
    }
    return g;
}

double calc_graph_sum(const GraphDouble &g, bool isAbs) {
    double res = 0.;
    GraphProperty *gp = g.m_property.get();
    int nx = gp->nx;
    int ny = gp->ny;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double t = g.m_vertices[i * ny + j].m_property.m_value;
            res += isAbs ? fabs(t) : t;
        }
    }
    return res;
}