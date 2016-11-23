#ifndef FEM_CIRCLE_GRAPH_UTILS_H
#define FEM_CIRCLE_GRAPH_UTILS_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

using namespace boost;

typedef adjacency_list<vecS, vecS, undirectedS> Graph;

void print_graph(const char *filename, Graph g);
void generate_png(const char *filename);

#endif //FEM_CIRCLE_GRAPH_UTILS_H