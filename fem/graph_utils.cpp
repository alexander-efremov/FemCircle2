#include "graph_utils.h"

void print_graph(const char *filename, Graph g, bool generate_png)
{
    std::ofstream out;
    out.open(filename);
    write_graphviz(out, g);
    out.close();

    if(generate_png)
    {

    }
}
