//
// Created by efremov_aa on 05.12.2016.
//
#include <utils.h>
#include <graphs.h>
#include <common.h>
#include "gtest/gtest.h"
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/copy.hpp>

class GraphFixture : public ::testing::Test {
protected:
    virtual void SetUp() {
    }

public:
    GraphFixture() : Test() {
    }
};

using namespace boost;

template<typename TGraph>
class bfs_visitor : public default_bfs_visitor {
public:
    bfs_visitor(TGraph& g) : _g(g) {}

    template<typename Vertex>
    void discover_vertex(Vertex u, const TGraph &g) const {
        _g.m_vertices[u].m_property.m_value = 99;
    }
private:
    TGraph& _g;
};

TEST_F(GraphFixture, graph_bfs) {
    double d = 10.;
    NX = (unsigned int) d;
    NY = (unsigned int) d;
    NX_1 = NX + 1;
    NY_1 = NY + 1;
    XY = NX_1 * NY_1;
    auto *g = create_graph_as_grid(NX_1, NY_1, A, B, C, D, V, U, TAU, R_LVL, HX_SMALLEST, HY_SMALLEST, HX, HY, IDEAL_SQ_SIZE_X, IDEAL_SQ_SIZE_Y, 1);
    auto &gr = *g;
    //breadth_first_search(gr, vertex(0, gr), visitor(bfs_solver_visitor1<Graph>(gr)));
    auto it = boost::vertices(*g);
    // найдем узел графа, который содержит эту точку
    for (auto iti = it.first; iti != it.second; ++iti) {
        std::cout << *iti << std::endl;
    }
//    print_graph("graph_bfs.dot", gr);
//    generate_png("graph_bfs.dot", "graph_bfs.png");
    delete g;
}

// We need a bool to store visibility of each edge
struct EdgeVisibility {
    bool visible;
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, EdgeVisibility> MyGraph;

void TestCopy(const MyGraph &g) {
    MyGraph copiedMyGraph = g;
    unsigned int counter = 0;
    auto copyIter = boost::edges(copiedMyGraph);
    for (auto it = copyIter.first; it != copyIter.second; ++it) {
        if (counter == 10) copiedMyGraph[*it].visible = false;
        else counter++;
    }

    std::cout << "Input MyGraph edge visiblity: " << std::endl;
    auto gEdgeIteratorRange = boost::edges(g);
    for (MyGraph::edge_iterator it1 = gEdgeIteratorRange.first; it1 != gEdgeIteratorRange.second; ++it1) {
        std::cout << g[*it1].visible << " ";
    }
    std::cout << std::endl;
    std::cout << "Copied MyGraph edge visiblity: " << std::endl;

    auto copyIt2 = boost::edges(copiedMyGraph);
    for (auto it2 = copyIt2.first; it2 != copyIt2.second; ++it2) {
        std::cout << copiedMyGraph[*it2].visible << " ";
    }
    std::cout << std::endl;
}

MyGraph CreateMyGraph() {
    MyGraph g(19);

    boost::add_edge(0, 1, g);
    boost::add_edge(1, 2, g);
    boost::add_edge(2, 3, g);
    boost::add_edge(3, 4, g);
    boost::add_edge(4, 5, g);
    boost::add_edge(5, 6, g);
    boost::add_edge(6, 7, g);
    boost::add_edge(7, 8, g);
    boost::add_edge(7, 9, g);
    boost::add_edge(9, 10, g);
    boost::add_edge(7, 11, g);
    boost::add_edge(11, 12, g);
    boost::add_edge(12, 13, g);
    boost::add_edge(13, 14, g);
    boost::add_edge(13, 15, g);
    boost::add_edge(15, 16, g);
    boost::add_edge(16, 17, g);
    boost::add_edge(17, 18, g);

    // Set all edges to visible
    std::pair<MyGraph::edge_iterator, MyGraph::edge_iterator> edgeIteratorRange = boost::edges(g);
    for (MyGraph::edge_iterator edgeIterator = edgeIteratorRange.first;
         edgeIterator != edgeIteratorRange.second; ++edgeIterator) {
        g[*edgeIterator].visible = true;
    }

    return g;
}

TEST_F(GraphFixture, GraphCopy) {
    MyGraph gt = CreateMyGraph();
    TestCopy(gt);
}