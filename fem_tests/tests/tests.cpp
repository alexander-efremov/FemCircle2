#include <utils.h>
#include <graphs.h>
#include <common.h>
#include "gtest/gtest.h"
#include <boost/graph/breadth_first_search.hpp>

#include <iostream>


class FemFixture : public ::testing::Test {
protected:
    virtual void TearDown() {
        if (G1 != NULL)
            delete[] G1;
        if (G2 != NULL)
            delete[] G2;
        if (G3 != NULL)
            delete[] G3;
        if (G4 != NULL)
            delete[] G4;
    }

    virtual void SetUp() {
    }

public:
    FemFixture() : Test() {
    }
};

void print_params() {
    printf("\nNXxNY = %dx%d\n", NX, NY);
    printf("(U, V) = (%le, %le)\n", U, V);
    printf("(HX, HY) = (%le, %le)\n", HX, HY);
    printf("TAU = %le\n", TAU);
    printf("TIME_STEP_CNT = %d\n", TIME_STEP_CNT);
    printf("INTEGR_TYPE = %d\n", INTEGR_TYPE);
    printf("IDEAL_SQ_SIZE = %dx%d\n", IDEAL_SQ_SIZE_X, IDEAL_SQ_SIZE_Y);
    printf("CENTER_OFFSET = %le, %le\n", CENTER_OFFSET_X, CENTER_OFFSET_Y);
    printf("NX3 = %d\n", NX3);
    printf("NX3_1 = %d\n", NX3_1);
    printf("NY3 = %d\n", NY3);
    printf("NY3_1 = %d\n", NY3_1);
    printf("XY = %d\n", XY);
    printf("R_LVL = %d\n", R_LVL);
    printf("R = %d\n", R);
    printf("EPS_GRID = %e\n", EPS_GRID);
    printf("RES_EPS = %e\n", RES_EPS);
    printf("APPROX_TYPE = %d\n", APPROX_TYPE);
    fflush(stdout);
}

void print_params_1() {
    printf("\nNXxNY = %dx%d\n", NX, NY);
    printf("(U, V) = (%le, %le)\n", U, V);
    printf("(HX, HY) = (%le, %le)\n", HX, HY);
    printf("TAU = %le\n", TAU);
    printf("TIME_STEP_CNT = %d\n", TIME_STEP_CNT);
    printf("INTEGR_TYPE = %d\n", INTEGR_TYPE);
    printf("IDEAL_SQ_SIZE = %dx%d\n", IDEAL_SQ_SIZE_X, IDEAL_SQ_SIZE_Y);
    printf("CENTER_OFFSET = %le, %le\n", CENTER_OFFSET_X, CENTER_OFFSET_Y);
    printf("XY = %d\n", XY);
    printf("R_LVL = %d\n", R_LVL);
    printf("R = %d\n", R);
    printf("EPS_GRID = %e\n", EPS_GRID);
    printf("RES_EPS = %e\n", RES_EPS);
    printf("APPROX_TYPE = %d\n", APPROX_TYPE);
    fflush(stdout);
}

void init_boundary_arrays_and_cp(int nx, int ny) {
    G1 = new int[nx];
    G2 = new int[ny];
    G3 = new int[nx];
    G4 = new int[ny];
    for (int i = 0; i < nx; ++i) {
        G1[i] = 0;
        G3[i] = 0;
    }
    for (int j = 0; j < ny; ++j) {
        G2[j] = 0;
        G4[j] = 0;
    }
    CP00 = 0;
    CP10 = 0;
    CP01 = 0;
    CP11 = 0;
}

void run_solver_1(unsigned int d) {
    assert(d >= 50);

    A = 0.;
    B = 1.;
    C = 0.;
    D = 1.;
    R_SQ = 0.099 * 0.099;
    INN_DENSITY = 1.;
    OUT_DENSITY = 0.;
    R_LVL = 2;

    NX = d;
    NY = d;
    NX_1 = NX + 1;
    NY_1 = NY + 1;
    HX = (B - A) / NX;
    HY = (D - C) / NY;
    IDEAL_SQ_SIZE_X = 64;
    IDEAL_SQ_SIZE_Y = 64;
    EPS_GRID = 0.5;
    RES_EPS = 1.e-9;

    int sz = d;
    sz = sz * ((int) std::pow(3., (double) R_LVL));
    NX3 = sz;
    NY3 = sz;
    NX3_1 = NX3 + 1;
    NY3_1 = NY3 + 1;
    XY = NX3_1 * NY3_1;
    R = (int) std::pow(3., (double) R_LVL);
    HX_SMALLEST = (B - A) / (NX * std::pow(3., (double) R_LVL));
    HY_SMALLEST = (D - C) / (NY * std::pow(3., (double) R_LVL));

    CENTER_OFFSET_X = 0.3;
    CENTER_OFFSET_Y = 0.3;

    U = 1.;
    V = 1.;
    TAU = 1.e-3;

    TIME_STEP_CNT = 5;

    APPROX_TYPE = 1;

    init_boundary_arrays_and_cp(NX3_1, NY3_1);
    print_params();

    int *grid = new int[XY];
    int *gridPr = new int[XY];

    double *density = solve_1(grid, gridPr);
    double *exact0 = calc_exact_1(grid, 0, NX3_1, NY3_1, HX_SMALLEST, HY_SMALLEST, R_LVL);
    double *exactT = calc_exact_1(grid, TAU * TIME_STEP_CNT, NX3_1, NY3_1, HX_SMALLEST, HY_SMALLEST, R_LVL);

    double x0 = get_center_x();
    double y0 = get_center_y();
    print_surface("exact", NX, NY, HX, HY, 0, A, C, x0, y0, TAU, U, V, exact0);
    print_surface("exact", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U, V, exactT);

    delete[] density;
    delete[] exact0;
    delete[] exactT;
    delete[] grid;
    delete[] gridPr;
}

// убрано притягивание сетки
// убрана рекурсия
// деревья вместо массивов?
//
TEST_F(FemFixture, test1) {
    for (int i = 1; i < 2; ++i) {
        double d = 0;
        switch (i) {
            case 0:
                d = 50.;
                break;
            case 1:
                d = 100.;
                break;
            case 2:
                d = 200.;
                break;
            case 3:
                d = 400.;
                break;
            case 4:
                d = 800.;
                break;
            case 5:
                d = 1600.;
                break;
            default:
                return;
        }
        run_solver_1((unsigned int) d);
    }
}


void run_solver_2(unsigned int d) {
    assert(d >= 50);

    A = 0.;
    B = 1.;
    C = 0.;
    D = 1.;
    R_SQ = 0.099 * 0.099;
    INN_DENSITY = 1.;
    OUT_DENSITY = 0.;
    U = 1.;
    V = 1.;
    TAU = 1.e-3;
    TIME_STEP_CNT = 5;
    APPROX_TYPE = 1;
    NX = d;
    NY = d;
    NX_1 = NX + 1;
    NY_1 = NY + 1;
    HX = (B - A) / NX;
    HY = (D - C) / NY;
    IDEAL_SQ_SIZE_X = 64;
    IDEAL_SQ_SIZE_Y = 64;
    RES_EPS = 1.e-9;

    R_LVL = 2;
    HX_SMALLEST = (B - A) / (NX * std::pow(3., (double) R_LVL));
    HY_SMALLEST = (D - C) / (NY * std::pow(3., (double) R_LVL));
    EPS_GRID = 0.5;

    int sz = d;
    sz = sz * ((int) std::pow(3., (double) R_LVL));
    NX3 = sz;
    NY3 = sz;
    NX3_1 = NX3 + 1;
    NY3_1 = NY3 + 1;
    XY = NX3_1 * NY3_1;
    R = (int) std::pow(3., (double) R_LVL);

    CENTER_OFFSET_X = 0.3;
    CENTER_OFFSET_Y = 0.3;

    //init_boundary_arrays_and_cp(NX3_1, NY3_1);
    print_params_1();

    Graph *density = solve_2();
//    Graph* exact0 = calc_exact_2(gridRef, 0, NX3_1, NY3_1, HX_SMALLEST, HY_SMALLEST, R_LVL);
//    Graph* exactT = calc_exact_2(gridRef, TAU * TIME_STEP_CNT, NX3_1, NY3_1, HX_SMALLEST, HY_SMALLEST, R_LVL);

//    double x0 = get_center_x();
//    double y0 = get_center_y();
//    print_surface("exact", NX, NY, HX, HY, 0, A, C, x0, y0, TAU, U, V, exact0);
//    print_surface("exact", NX, NY, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U, V, exactT);

    delete density;
//    delete exact0;
//    delete exactT;
}

// убрано притягивание сетки
// убрана рекурсия
// деревья вместо массивов?
//
TEST_F(FemFixture, test2) {
    for (int i = 1; i < 2; ++i) {
        unsigned int d = 0;
        switch (i) {
            case 0:
                d = 50;
                break;
            case 1:
                d = 100;
                break;
            case 2:
                d = 200;
                break;
            case 3:
                d = 400;
                break;
            case 4:
                d = 800;
                break;
            case 5:
                d = 1600;
                break;
            default:
                return;
        }
        run_solver_2(d);
    }
}

TEST_F(FemFixture, graph) {
    double d = 10.;
    NX = (unsigned int) d;
    NY = (unsigned int) d;
    NX_1 = NX + 1;
    NY_1 = NY + 1;
    XY = NX_1 * NY_1;

    Graph *g = create_graph_as_grid(NX_1, NY_1, A, B, C, D, V, U, TAU, R_LVL, HX_SMALLEST, HY_SMALLEST, HX, HY, 1);
    Graph &gr = *g;
    bool b = is_graph_connected(gr);

    if (b) std::cout << "\nIs connected graph" << std::endl;
    else std::cout << "\nIs NOT connected graph" << std::endl;

    print_graph("grid.dot", gr);
    generate_png("grid.dot", "grid.png");

    delete g;
}

using namespace boost;

template<typename TGraph>
class bfs_solver_visitor : public default_bfs_visitor {
public:
    bfs_solver_visitor(TGraph& g) : _g(g) {}

    template<typename Vertex>
    void discover_vertex(Vertex u, const TGraph &g) const {
        _g.m_vertices[u].m_property.m_value = 99;
    }
private:
    TGraph& _g;
};


TEST_F(FemFixture, graph_bfs) {
    using namespace boost;
    double d = 10.;
    NX = (unsigned int) d;
    NY = (unsigned int) d;
    NX_1 = NX + 1;
    NY_1 = NY + 1;
    XY = NX_1 * NY_1;
    auto *g = create_graph_as_grid(NX_1, NY_1, A, B, C, D, V, U, TAU, R_LVL, HX_SMALLEST, HY_SMALLEST, HX, HY, 1);
    auto &gr = *g;
    breadth_first_search(gr, vertex(0, gr), visitor(bfs_solver_visitor<Graph>(gr)));
    print_graph("graph_bfs.dot", gr);
    generate_png("graph_bfs.dot", "graph_bfs.png");
    delete g;
}