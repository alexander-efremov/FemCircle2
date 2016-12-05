#include <string>
#include <iostream>
#include "graphs.h"
#include "utils.h"
#include <boost/graph/connected_components.hpp>

using namespace std;

bool is_graph_connected(const Graph &g) {
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
        int px = v % t->nx_1;
        int py = v / t->ny_1;
        int x = _g.m_vertices[v].m_property.m_i;
        int y = _g.m_vertices[v].m_property.m_j;
        double value = _g.m_vertices[v].m_property.m_value;
        out << " ["
            << " shape=point"
            << " xlabel=\"(" << x << "," << y << "," << value << ")\""
            << "pos=\"" << px << "," << t->ny_1 - py << "!\""
            << (is_corner_node(px, py, t->nx_1, t->ny_1) ? " fillcolor=\"red\" style=\"filled\"" : "")
            << (is_border_node(px, py, t->nx_1, t->ny_1) ? " fillcolor=\"blue\" style=\"filled\"" : "")
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

void print_graph(const char *filename, const Graph &g) {
    std::ofstream out;
    out.open(filename);
    write_graphviz(out, g, VertexPropertyWriter(g), default_writer(), GraphPropertyWriter());
    out.close();
}

Graph *create_graph_as_grid(const Graph &graph, double defaultValue) {
    GraphProperty *gp = graph.m_property.get();
    return create_graph_as_grid(gp->nx_1, gp->ny_1, gp->a, gp->b, gp->c, gp->d, gp->v, gp->u, gp->tau, gp->r_lvl,
                                gp->hx_min, gp->hy_min, gp->hx, gp->hy, gp->ideal_square_size_nx,
                                gp->ideal_square_size_ny,
                                defaultValue);
}

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
        int r_lvl,
        double hx_min,
        double hy_min,
        double hx,
        double hy,
        unsigned int ideal_square_size_nx,
        unsigned int ideal_square_size_ny,
        double defaultValue) {
    assert(nx_1 > 0);
    assert(ny_1 > 0);

    Graph *g = new Graph(nx_1 * ny_1, GraphProperty(nx_1, ny_1, a, b, c, d, u, v, tau, r_lvl, hx_min, hy_min, hx, hy,
                                                    ideal_square_size_nx, ideal_square_size_ny));
    for (int i = 0; i < nx_1; ++i) {
        int stride = i * nx_1;
        pair<detail::adj_list_gen<boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexPropertyDouble, boost::no_property, GraphProperty>, boost::vecS, boost::vecS, boost::undirectedS, VertexPropertyDouble, boost::no_property, GraphProperty, boost::listS>::config::edge_descriptor, bool> e;
        for (int j = 0; j < ny_1 - 1; ++j) {
            e = add_edge((unsigned long long int) (j + stride), (unsigned long long int) (j + stride + 1), *g);
            g->m_vertices[e.first.m_source].m_property.m_value = defaultValue;
            g->m_vertices[e.first.m_source].m_property.m_i = i;
            g->m_vertices[e.first.m_source].m_property.m_j = j;
        }
        g->m_vertices[e.first.m_target].m_property.m_i = i;
        g->m_vertices[e.first.m_target].m_property.m_j = ny_1 - 1;
    }

    for (int j = 0; j < ny_1 - 1; ++j) {
        int stride = j * ny_1;
        for (int i = 0; i < nx_1; i++) {
            add_edge((unsigned long long int) (i + stride), (unsigned long long int) (i + stride + nx_1), *g);
        }
    }
    return g;
}

double calc_graph_sum(const Graph &g, bool isAbs) {
    double res = 0.;
    GraphProperty *gp = g.m_property.get();
    int nx = gp->nx_1;
    int ny = gp->ny_1;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double t = g.m_vertices[i * ny + j].m_property.m_value;
            res += isAbs ? fabs(t) : t;
        }
    }
    return res;
}