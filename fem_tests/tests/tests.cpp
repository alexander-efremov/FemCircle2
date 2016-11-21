#include <utils.h>
#include <cmath>
#include <common.h>
#include <fstream>
#include "gtest/gtest.h"
#include <algorithm>
#include <locale>

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

// убрано притягивание сетки
// убрана рекурсия
// деревья вместо массивов?
//
TEST_F(FemFixture, test1) {
    double tme = 0.;
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

        A = 0.;
        B = 1.;
        C = 0.;
        D = 1.;
        R_SQ = 0.099 * 0.099;
        INN_DENSITY = 1.;
        OUT_DENSITY = 0.;
        R_LVL = 2;

        int sz = (int) d;
        sz = sz * ((int) std::pow(3., R_LVL));

        NX = (int) d;
        NY = (int) d;
        NX_1 = NX + 1;
        NY_1 = NY + 1;
        HX = (B - A) / NX;
        HY = (D - C) / NY;
        IDEAL_SQ_SIZE_X = 64;
        IDEAL_SQ_SIZE_Y = 64;
        EPS_GRID = 0.5;
        RES_EPS = 1.e-9;

        NX3 = sz;
        NY3 = sz;
        NX3_1 = NX3 + 1;
        NY3_1 = NY3 + 1;
        R = (int) std::pow(3., R_LVL);
        HX_SMALLEST = (B - A) / (NX * std::pow(3., R_LVL));
        HY_SMALLEST = (D - C) / (NY * std::pow(3., R_LVL));

        CENTER_OFFSET_X = 0.3;
        CENTER_OFFSET_Y = 0.3;

        U = 1.;
        V = 1.;
        TAU = 1.e-3;

        TIME_STEP_CNT = 5;
        XY = NX3_1 * NY3_1;

        APPROX_TYPE = 2;

        init_boundary_arrays_and_cp(NX3_1, NY3_1);
        print_params();

        int *grid = new int[XY];
        int *gridPr = new int[XY];

        double *density = solve_1(tme, grid, gridPr);
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
}