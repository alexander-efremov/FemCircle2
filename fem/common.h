#ifndef FEM_CIRCLE_COMMON_H
#define FEM_CIRCLE_COMMON_H

#include "consts.h"
#include "graphs.h"

double func_u(double t, double x, double y);
double func_v(double t, double x, double y);

// ссылка на тест https://ideone.com/dRQkab
// откуда взято http://www.sunshine2k.de/coding/java/PointOnLine/PointOnLine.html
static void get_projected_point_on_line_fast(double start_x, double start_y, double end_x, double end_y, double px1,
                                             double py1,
                                             double &x, double &y) {
    // get dot product of e1, e2
    double e1x = end_x - start_x;
    double e1y = end_y - start_y;
    double e2x = px1 - start_x;
    double e2y = py1 - start_y;
    double valDp = e1x * e2x - e1y * e2y;

    // get squared length of e1
    double len2 = e1x * e1x + e1y * e1y;

    x = start_x + (valDp * e1x) / len2;
    y = start_y + (valDp * e1y) / len2;
//    x = px1;
//    y = py1;
}

inline static void get_coord_new_hor_1(double x11_orig, double y11_orig, double x22_orig, double y22_orig, double x33_orig,
                                double y33_orig, double x44_orig, double y44_orig, double x55_orig, double y55_orig,
                                double x66_orig, double y66_orig,
                                double x11, double y11, double x22, double y22, double x33, double y33, double x44,
                                double y44, double x55, double y55, double x66, double y66,
                                int p, int q, double (&points)[4][7][2],
                                double &x1, double &y1, double &x2, double &y2,
                                double &x3, double &y3, double &x4, double &y4, double hx_lev, double hy_lev, int tl,
                                double tau) {
    double time = tl * tau;

    if (p == -1 && q == -1) { // cell -1-1
        // 1
        x1 = x11;
        y1 = y11;
        points[0][0][0] = x1;
        points[0][0][1] = y1;
        // 2
        double t_x2 = x11_orig + hx_lev;
        double t_y2 = y11_orig;
        t_x2 = t_x2 - tau * func_u(time, t_x2, t_y2);
        t_y2 = t_y2 - tau * func_v(time, t_x2, t_y2);
        get_projected_point_on_line_fast(x11, y11, x22, y22, t_x2, t_y2, x2, y2);
        points[0][1][0] = x2;
        points[0][1][1] = y2;
        // 3
        double t_x3 = x11_orig + hx_lev;
        double t_y3 = y11_orig + hy_lev;
        x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        points[1][1][0] = x3;
        points[1][1][1] = y3;
        // 4
        double t_x4 = x11_orig;
        double t_y4 = y11_orig + hy_lev;
        t_x4 = t_x4 - tau * func_u(time, t_x4, t_y4);
        t_y4 = t_y4 - tau * func_v(time, t_x4, t_y4);
        get_projected_point_on_line_fast(x11, y11, x66, y66, t_x4, t_y4, x4, y4);
        points[1][0][0] = x4;
        points[1][0][1] = y4;
    }
    else if (p == -1 && q == 0) { // cell -10
        // 1
        x1 = points[1][0][0];
        y1 = points[1][0][1];
        // 2
        x2 = points[1][1][0];
        y2 = points[1][1][1];
        // 3
        double t_x3 = x11_orig + hx_lev;
        double t_y3 = y11_orig + 2. * hy_lev;
        x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        points[2][1][0] = x3;
        points[2][1][1] = y3;
        // 4
        double t_x4 = x11_orig;
        double t_y4 = y11_orig + 2. * hy_lev;
        t_x4 = t_x4 - tau * func_u(time, t_x4, t_y4);
        t_y4 = t_y4 - tau * func_v(time, t_x4, t_y4);
        get_projected_point_on_line_fast(x11, y11, x66, y66, t_x4, t_y4, x4, y4);
        points[2][0][0] = x4;
        points[2][0][1] = y4;
    }
    else if (p == -1 && q == 1) { // cell -11
        // 1
        x1 = points[2][0][0];
        y1 = points[2][0][1];
        // 2
        x2 = points[2][1][0];
        y2 = points[2][1][1];
        // 3
        double t_x3 = x66_orig + hx_lev;
        double t_y3 = y66_orig;
        t_x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        t_y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        get_projected_point_on_line_fast(x66, y66, x55, y55, t_x3, t_y3, x3, y3);
        points[3][1][0] = x3;
        points[3][1][1] = y3;
        // 4
        x4 = x66;
        y4 = y66;
        points[3][0][0] = x4;
        points[3][0][1] = y4;
    }
    else if (p == 0 && q == -1) { // cell 0-1
        // 1
        x1 = points[0][1][0];
        y1 = points[0][1][1];
        // 2
        double t_x2 = x11_orig + 2. * hx_lev;
        double t_y2 = y11_orig;
        t_x2 = t_x2 - tau * func_u(time, t_x2, t_y2);
        t_y2 = t_y2 - tau * func_v(time, t_x2, t_y2);
        get_projected_point_on_line_fast(x11, y11, x22, y22, t_x2, t_y2, x2, y2);
        points[0][2][0] = x2;
        points[0][2][1] = y2;
        // 3
        double t_x3 = x11_orig + 2. * hx_lev;
        double t_y3 = y11_orig + hy_lev;
        x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        points[1][2][0] = x3;
        points[1][2][1] = y3;
        // 4
        x4 = points[1][1][0];
        y4 = points[1][1][1];
    }
    else if (p == 0 && q == 0) {  // cell 00
        // 1
        x1 = points[1][1][0];
        y1 = points[1][1][1];
        // 2
        x2 = points[1][2][0];
        y2 = points[1][2][1];
        // 3
        double t_x3 = x11_orig + 2. * hx_lev;
        double t_y3 = y11_orig + 2. * hy_lev;
        x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        points[2][2][0] = x3;
        points[2][2][1] = y3;
        // 4
        x4 = points[2][1][0];
        y4 = points[2][1][1];
    }
    else if (p == 0 && q == 1) { // cell 01
        // 1
        x1 = points[2][1][0];
        y1 = points[2][1][1];
        // 2
        x2 = points[2][2][0];
        y2 = points[2][2][1];
        // 3
        double t_x3 = x66_orig + 2. * hx_lev;
        double t_y3 = y66_orig;
        t_x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        t_y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        get_projected_point_on_line_fast(x66, y66, x55, y55, t_x3, t_y3, x3, y3);
        points[3][2][0] = x3;
        points[3][2][1] = y3;
        // 4
        x4 = points[3][1][0];
        y4 = points[3][1][1];
    }
    else if (p == 1 && q == -1) { // cell 1-1
        // 1
        x1 = points[0][2][0];
        y1 = points[0][2][1];
        // 2
        x2 = x22;
        y2 = y22;
        points[0][3][0] = x2;
        points[0][3][1] = y2;
        // 3
        double t_x3 = x22_orig;
        double t_y3 = y22_orig + hy_lev;
        x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        points[1][3][0] = x3;
        points[1][3][1] = y3;
        // 4
        x4 = points[1][2][0];
        y4 = points[1][2][1];
    }
    else if (p == 1 && q == 0) { // cell 10
        // 1
        x1 = points[1][2][0];
        y1 = points[1][2][1];
        // 2
        x2 = points[1][3][0];
        y2 = points[1][3][1];
        // 3
        double t_x3 = x22_orig;
        double t_y3 = y22_orig + 2. * hy_lev;
        x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        points[2][3][0] = x3;
        points[2][3][1] = y3;
        // 4
        x4 = points[2][2][0];
        y4 = points[2][2][1];
    }
    else if (p == 1 && q == 1) { // cell 11
        // 1
        x1 = points[2][2][0];
        y1 = points[2][2][1];
        // 2
        x2 = points[2][3][0];
        y2 = points[2][3][1];
        // 3
        x3 = x33;
        y3 = y33;
        points[3][3][0] = x3;
        points[3][3][1] = y3;
        // 4
        x4 = points[3][2][0];
        y4 = points[3][2][1];
    }
    else if (p == 2 && q == -1) { // cell 2-1
        // 1
        x1 = points[0][3][0];
        y1 = points[0][3][1];
        // 2
        double t_x2 = x22_orig + hx_lev;
        double t_y2 = y22_orig;
        t_x2 = t_x2 - tau * func_u(time, t_x2, t_y2);
        t_y2 = t_y2 - tau * func_v(time, t_x2, t_y2);
        get_projected_point_on_line_fast(x22, y22, x33, y33, t_x2, t_y2, x2, y2);
        points[0][4][0] = x2;
        points[0][4][1] = y2;
        // 3
        double t_x3 = x22_orig + hx_lev;
        double t_y3 = y22_orig + hy_lev;
        x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        points[1][4][0] = x3;
        points[1][4][1] = y3;
        // 4
        x4 = points[1][3][0];
        y4 = points[1][3][1];
    }
    else if (p == 2 && q == 0) { // cell 20
        // 1
        x1 = points[1][3][0];
        y1 = points[1][3][1];
        // 2
        x2 = points[1][4][0];
        y2 = points[1][4][1];
        // 3
        double t_x3 = x22_orig + hx_lev;
        double t_y3 = y22_orig + 2. * hy_lev;
        x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        points[2][4][0] = x3;
        points[2][4][1] = y3;
        // 4
        x4 = points[2][3][0];
        y4 = points[2][3][1];
    }
    else if (p == 2 && q == 1) { // cell 21
        // 1
        x1 = points[2][3][0];
        y1 = points[2][3][1];
        // 2
        x2 = points[2][4][0];
        y2 = points[2][4][1];
        // 3
        double t_x3 = x55_orig + hx_lev;
        double t_y3 = y55_orig;
        t_x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        t_y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        get_projected_point_on_line_fast(x55, y55, x44, y44, t_x3, t_y3, x3, y3);
        points[3][4][0] = x3;
        points[3][4][1] = y3;
        // 4
        x4 = x55;
        y4 = y55;
        points[3][3][0] = x4;
        points[3][3][1] = y4;
    }
    else if (p == 3 && q == -1) { // cell 3-1
        // 1
        x1 = points[0][4][0];
        y1 = points[0][4][1];
        // 2
        double t_x2 = x22_orig + 2. * hx_lev;
        double t_y2 = y22_orig;
        t_x2 = t_x2 - tau * func_u(time, t_x2, t_y2);
        t_y2 = t_y2 - tau * func_v(time, t_x2, t_y2);
        get_projected_point_on_line_fast(x22, y22, x33, y33, t_x2, t_y2, x2, y2);
        points[0][5][0] = x2;
        points[0][5][1] = y2;
        // 3
        double t_x3 = x22_orig + 2. * hx_lev;
        double t_y3 = y22_orig + hy_lev;
        x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        points[1][5][0] = x3;
        points[1][5][1] = y3;
        // 4
        x4 = points[1][4][0];
        y4 = points[1][4][1];
    }
    else if (p == 3 && q == 0) { // cell 30
        // 1
        x1 = points[1][4][0];
        y1 = points[1][4][1];
        // 2
        x2 = points[1][5][0];
        y2 = points[1][5][1];
        // 3
        double t_x3 = x22_orig + 2. * hx_lev;
        double t_y3 = y22_orig + 2. * hy_lev;
        x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        points[2][5][0] = x3;
        points[2][5][1] = y3;
        // 4
        x4 = points[2][4][0];
        y4 = points[2][4][1];
    }
    else if (p == 3 && q == 1) { // cell 31
        // 1
        x1 = points[2][4][0];
        y1 = points[2][4][1];
        // 2
        x2 = points[2][5][0];
        y2 = points[2][5][1];
        // 3
        double t_x3 = x55_orig + 2. * hx_lev;
        double t_y3 = y55_orig;
        t_x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        t_y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        get_projected_point_on_line_fast(x55, y55, x44, y44, t_x3, t_y3, x3, y3);
        points[3][5][0] = x3;
        points[3][5][1] = y3;
        // 4
        x4 = points[3][4][0];
        y4 = points[3][4][1];
    }
    else if (p == 4 && q == -1) { // cell 4-1
        // 1
        x1 = points[0][5][0];
        y1 = points[0][5][1];
        // 2
        x2 = x33;
        y2 = y33;
        points[0][6][0] = x2;
        points[0][6][1] = y2;
        // 3
        double t_x3 = x33_orig;
        double t_y3 = y33_orig + hy_lev;
        t_x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        t_y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        get_projected_point_on_line_fast(x33, y33, x44, y44, t_x3, t_y3, x3, y3);
        points[1][6][0] = x3;
        points[1][6][1] = y3;
        // 4
        x4 = points[1][5][0];
        y4 = points[1][5][1];
    }
    else if (p == 4 && q == 0) { // cell 40
        // 1
        x1 = points[1][5][0];
        y1 = points[1][5][1];
        // 2
        x2 = points[1][6][0];
        y2 = points[1][6][1];
        // 3
        double t_x3 = x33_orig;
        double t_y3 = y33_orig + 2. * hy_lev;
        t_x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        t_y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        get_projected_point_on_line_fast(x33, y33, x44, y44, t_x3, t_y3, x3, y3);
        points[2][6][0] = x3;
        points[2][6][1] = y3;
        // 4
        x4 = points[2][5][0];
        y4 = points[2][5][1];
    }
    else if (p == 4 && q == 1) {  // cell 41
        // 1
        x1 = points[2][5][0];
        y1 = points[2][5][1];
        // 2
        x2 = points[2][6][0];
        y2 = points[2][6][1];
        // 3
        x3 = x44;
        y3 = y44;
        points[3][6][0] = x3;
        points[3][6][1] = y3;
        // 4
        x4 = points[3][5][0];
        y4 = points[3][5][1];
    }
}

inline static void get_coord_new_vert_1(double x11_orig, double y11_orig, double x22_orig, double y22_orig, double x33_orig,
                                 double y33_orig, double x44_orig, double y44_orig, double x55_orig, double y55_orig,
                                 double x66_orig, double y66_orig,
                                 double x11, double y11, double x22, double y22, double x33, double y33, double x44,
                                 double y44, double x55, double y55, double x66, double y66,
                                 int p, int q, double (&points)[7][4][2],
                                 double &x1, double &y1, double &x2, double &y2,
                                 double &x3, double &y3, double &x4, double &y4, double hx_lev, double hy_lev, int tl,
                                 double tau) {
    double time = tl * tau;

    if (p == -1 && q == -1) { // cell -1-1
        // 1
        x1 = x11;
        y1 = y11;
        points[0][0][0] = x1;
        points[0][0][1] = y1;
        // 2
        double t_x2 = x11_orig + hx_lev;
        double t_y2 = y11_orig;
        t_x2 = t_x2 - tau * func_u(time, t_x2, t_y2);
        t_y2 = t_y2 - tau * func_v(time, t_x2, t_y2);
        get_projected_point_on_line_fast(x11, y11, x22, y22, t_x2, t_y2, x2, y2);
        points[0][1][0] = x2;
        points[0][1][1] = y2;
        // 3
        double t_x3 = x11_orig + hx_lev;
        double t_y3 = y11_orig + hy_lev;
        x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        points[1][1][0] = x3;
        points[1][1][1] = y3;
        // 4
        double t_x4 = x11_orig;
        double t_y4 = y11_orig + hy_lev;
        t_x4 = t_x4 - tau * func_u(time, t_x4, t_y4);
        t_y4 = t_y4 - tau * func_v(time, t_x4, t_y4);
        get_projected_point_on_line_fast(x11, y11, x66, y66, t_x4, t_y4, x4, y4);
        points[1][0][0] = x4;
        points[1][0][1] = y4;
    }
    else if (p == -1 && q == 0) { // cell -10
        // 1
        x1 = points[1][0][0];
        y1 = points[1][0][1];
        // 2
        x2 = points[1][1][0];
        y2 = points[1][1][1];
        // 3
        double t_x3 = x11_orig + hx_lev;
        double t_y3 = y11_orig + 2. * hy_lev;
        x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        points[2][1][0] = x3;
        points[2][1][1] = y3;
        // 4
        double t_x4 = x11_orig;
        double t_y4 = y11_orig + 2. * hy_lev;
        t_x4 = t_x4 - tau * func_u(time, t_x4, t_y4);
        t_y4 = t_y4 - tau * func_v(time, t_x4, t_y4);
        get_projected_point_on_line_fast(x11, y11, x66, y66, t_x4, t_y4, x4, y4);
        points[2][0][0] = x4;
        points[2][0][1] = y4;
    }
    else if (p == -1 && q == 1) { // cell -11
        // 1
        x1 = points[2][0][0];
        y1 = points[2][0][1];
        // 2
        x2 = points[2][1][0];
        y2 = points[2][1][1];
        // 3
        double t_x3 = x66_orig + hx_lev;
        double t_y3 = y66_orig;
        x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        points[3][1][0] = x3;
        points[3][1][1] = y3;
        // 4
        x4 = x66;
        y4 = y66;
        points[3][0][0] = x4;
        points[3][0][1] = y4;
    }
    else if (p == -1 && q == 2) { // cell -12
        // 1
        x1 = points[3][0][0];
        y1 = points[3][0][1];
        // 2
        x2 = points[3][1][0];
        y2 = points[3][1][1];
        // 3
        double t_x3 = x66_orig + hx_lev;
        double t_y3 = y66_orig + hy_lev;
        x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        points[4][1][0] = x3;
        points[4][1][1] = y3;
        // 4
        double t_x4 = x66_orig;
        double t_y4 = y66_orig + hy_lev;
        t_x4 = t_x4 - tau * func_u(time, t_x4, t_y4);
        t_y4 = t_y4 - tau * func_v(time, t_x4, t_y4);
        get_projected_point_on_line_fast(x66, y66, x55, y55, t_x4, t_y4, x4, y4);
        points[4][0][0] = x4;
        points[4][0][1] = y4;
    }
    else if (p == -1 && q == 3) {  // cell -13
        // 1
        x1 = points[4][0][0];
        y1 = points[4][0][1];
        // 2
        x2 = points[4][1][0];
        y2 = points[4][1][1];
        // 3
        double t_x3 = x66_orig + hx_lev;
        double t_y3 = y66_orig + 2. * hy_lev;
        x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        points[5][1][0] = x3;
        points[5][1][1] = y3;
        // 4
        double t_x4 = x66_orig;
        double t_y4 = y66_orig + 2. * hy_lev;
        t_x4 = t_x4 - tau * func_u(time, t_x4, t_y4);
        t_y4 = t_y4 - tau * func_v(time, t_x4, t_y4);
        get_projected_point_on_line_fast(x66, y66, x55, y55, t_x4, t_y4, x4, y4);
        points[5][0][0] = x4;
        points[5][0][1] = y4;
    }
    else if (p == -1 && q == 4) { // cell -14
        // 1
        x1 = points[5][0][0];
        y1 = points[5][0][1];
        // 2
        x2 = points[5][1][0];
        y2 = points[5][1][1];
        // 3
        double t_x3 = x55_orig + hx_lev;
        double t_y3 = y55_orig;
        t_x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        t_y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        get_projected_point_on_line_fast(x55, y55, x44, y44, t_x3, t_y3, x3, y3);
        points[6][1][0] = x3;
        points[6][1][1] = y3;
        // 4
        x4 = x55;
        y4 = y55;
        points[6][0][0] = x4;
        points[6][0][1] = y4;
    }
    else if (p == 0 && q == -1) { // cell 0-1
        // 1
        x1 = points[0][1][0];
        y1 = points[0][1][1];
        // 2
        double t_x2 = x11_orig + 2. * hx_lev;
        double t_y2 = y11_orig;
        t_x2 = t_x2 - tau * func_u(time, t_x2, t_y2);
        t_y2 = t_y2 - tau * func_v(time, t_x2, t_y2);
        get_projected_point_on_line_fast(x11, y11, x22, y22, t_x2, t_y2, x2, y2);
        points[0][2][0] = x2;
        points[0][2][1] = y2;
        // 3
        double t_x3 = x11_orig + 2. * hx_lev;
        double t_y3 = y11_orig + hy_lev;
        x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        points[1][2][0] = x3;
        points[1][2][1] = y3;
        // 4
        x4 = points[1][1][0];
        y4 = points[1][1][1];
    }
    else if (p == 0 && q == 0) { // cell 00
        // 1
        x1 = points[1][1][0];
        y1 = points[1][1][1];
        // 2
        x2 = points[1][2][0];
        y2 = points[1][2][1];
        // 3
        double t_x3 = x11_orig + 2. * hx_lev;
        double t_y3 = y11_orig + 2. * hy_lev;
        x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        points[2][2][0] = x3;
        points[2][2][1] = y3;
        // 4
        x4 = points[2][1][0];
        y4 = points[2][1][1];
    }
    else if (p == 0 && q == 1) { // cell 01
        // 1
        x1 = points[2][1][0];
        y1 = points[2][1][1];
        // 2
        x2 = points[2][2][0];
        y2 = points[2][2][1];
        // 3
        double t_x3 = x66 + 2. * hx_lev;
        double t_y3 = y66;
        x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        points[3][2][0] = x3;
        points[3][2][1] = y3;
        // 4
        x4 = points[3][1][0];
        y4 = points[3][1][1];
    }
    else if (p == 0 && q == 2) { // cell 02
        // 1
        x1 = points[3][1][0];
        y1 = points[3][1][1];
        // 2
        x2 = points[3][2][0];
        y2 = points[3][2][1];
        // 3
        double t_x3 = x66_orig + 2. * hx_lev;
        double t_y3 = y66_orig + hy_lev;
        x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        points[4][2][0] = x3;
        points[4][2][1] = y3;
        // 4
        x4 = points[4][1][0];
        y4 = points[4][1][1];
    }
    else if (p == 0 && q == 3) { // cell 03
        // 1
        x1 = points[4][1][0];
        y1 = points[4][1][1];
        // 2
        x2 = points[4][2][0];
        y2 = points[4][2][1];
        // 3
        double t_x3 = x66_orig + 2. * hx_lev;
        double t_y3 = y66_orig + 2. * hy_lev;
        x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        points[5][2][0] = x3;
        points[5][2][1] = y3;
        // 4
        x4 = points[5][1][0];
        y4 = points[5][1][1];
    }
    else if (p == 0 && q == 4) { // cell 04
        // 1
        x1 = points[5][1][0];
        y1 = points[5][1][1];
        // 2
        x2 = points[5][2][0];
        y2 = points[5][2][1];
        // 3
        double t_x3 = x55_orig + 2. * hx_lev;
        double t_y3 = y55_orig;
        t_x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        t_y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        get_projected_point_on_line_fast(x55, y55, x44, y44, t_x3, t_y3, x3, y3);
        points[6][2][0] = x3;
        points[6][2][1] = y3;
        // 4
        x4 = points[6][1][0];
        y4 = points[6][1][1];
    }
    else if (p == 1 && q == -1) { // cell 1-1
        // 1
        x1 = points[0][2][0];
        y1 = points[0][2][1];
        // 2
        x2 = x22;
        y2 = y22;
        points[0][3][0] = x2;
        points[0][3][1] = y2;
        // 3
        double t_x3 = x22_orig;
        double t_y3 = y22_orig + hy_lev;
        x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        get_projected_point_on_line_fast(x22, y22, x33, y33, t_x3, t_y3, x3, y3);
        points[1][3][0] = x3;
        points[1][3][1] = y3;
        // 4
        x4 = points[1][2][0];
        y4 = points[1][2][1];
    }
    else if (p == 1 && q == 0) { // cell 10
        // 1
        x1 = points[1][2][0];
        y1 = points[1][2][1];
        // 2
        x2 = points[1][3][0];
        y2 = points[1][3][1];
        // 3
        double t_x3 = x22_orig;
        double t_y3 = y22_orig + 2. * hy_lev;
        x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        get_projected_point_on_line_fast(x22, y22, x33, y33, t_x3, t_y3, x3, y3);
        points[2][3][0] = x3;
        points[2][3][1] = y3;
        // 4
        x4 = points[2][2][0];
        y4 = points[2][2][1];
    }
    else if (p == 1 && q == 1) { // cell 11
        // 1
        x1 = points[2][2][0];
        y1 = points[2][2][1];
        // 2
        x2 = points[2][3][0];
        y2 = points[2][3][1];
        // 3
        x3 = x33;
        y3 = y33;
        points[3][3][0] = x3;
        points[3][3][1] = y3;
        // 4
        x4 = points[3][3][0];
        y4 = points[3][3][1];
    }
    else if (p == 1 && q == 2) { // cell 12
        // 1
        x1 = points[3][2][0];
        y1 = points[3][2][1];
        // 2
        x2 = points[3][3][0];
        y2 = points[3][3][1];
        // 3
        double t_x3 = x33_orig;
        double t_y3 = y33_orig + hy_lev;
        t_x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        t_y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        get_projected_point_on_line_fast(x33, y33, x44, y44, t_x3, t_y3, x3, y3);
        points[4][3][0] = x3;
        points[4][3][1] = y3;
        // 4
        x4 = points[4][2][0];
        y4 = points[4][2][1];
    }
    else if (p == 1 && q == 3) { // cell 13
        // 1
        x1 = points[4][2][0];
        y1 = points[4][2][1];
        // 2
        x2 = points[4][3][0];
        y2 = points[4][3][1];
        // 3
        double t_x3 = x33_orig;
        double t_y3 = y33_orig + 2. * hy_lev;
        t_x3 = t_x3 - tau * func_u(time, t_x3, t_y3);
        t_y3 = t_y3 - tau * func_v(time, t_x3, t_y3);
        get_projected_point_on_line_fast(x33, y33, x44, y44, t_x3, t_y3, x3, y3);
        points[5][3][0] = x3;
        points[5][3][1] = y3;
        // 4
        x4 = points[5][2][0];
        y4 = points[5][2][1];
    }
    else if (p == 1 && q == 4) {  // cell 14
        // 1
        x1 = points[5][2][0];
        y1 = points[5][2][1];
        // 2
        x2 = points[5][3][0];
        y2 = points[5][3][1];
        // 3
        x3 = x44;
        y3 = y44;
        points[6][3][0] = x3;
        points[6][3][1] = y3;
        // 4
        x4 = points[6][2][0];
        y4 = points[6][2][1];
    }
}

inline double get_center_x() {
    return A + CENTER_OFFSET_X;
}

inline double get_center_y() {
    return C + CENTER_OFFSET_Y;
}

inline void get_left_area(int ii, int jj, double &x1, double &y1, double &x2, double &y2, double &x3, double &y3, double &x4,
                   double &y4, double hy, double hx_add, double hy_add, double a, double b, double c,
                   double d) {
    // p1 (A, y_{j-1/2})
    x1 = a;
    y1 = c + jj * hy - hy_add / 2.;
    // p2 (x_{1/2}, y_{j-1/2})
    x2 = a + hx_add / 2.;
    y2 = c + jj * hy - hy_add / 2.;
    //p3 (x_{1/2}, y_{j+1/2})
    x3 = a + hx_add / 2.;
    y3 = c + jj * hy + hy_add / 2.;
    //p4 (A, y_{j+1/2})
    x4 = a;
    y4 = c + jj * hy + hy_add / 2.;
    if (x1 < a || x1 >= b || x2 <= a || x2 >= b || x3 <= a || x3 >= b || x4 < a || x4 >= b
        || y1 <= c || y1 >= d || y2 <= c || y2 >= d || y3 <= c || y3 >= d || y4 <= c || y4 >= d) {
        printf("9. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
}

inline void get_top_area(int ii, int jj, double &x1, double &y1, double &x2, double &y2, double &x3, double &y3, double &x4,
                  double &y4, double hx, double hx_add, double hy_add, double a, double b, double c,
                  double d) {
    // p1 (x_{i-1/2}, y_{NY-1/2})
    x1 = a + ii * hx - hx_add / 2.;
    y1 = d - hy_add / 2.;
    // p2 (x_{i+1/2}, y_{NY-1/2})
    x2 = a + ii * hx + hx_add / 2.;
    y2 = d - hy_add / 2.;
    //p3 (x_{i+1/2}, D)
    x3 = a + ii * hx + hx_add / 2.;
    y3 = d;
    //p4 (x_{i-1/2}, D)
    x4 = a + ii * hx - hx_add / 2.;
    y4 = d;
    if (x1 <= a || x1 >= b || x2 <= a || x2 >= b || x3 <= a || x3 >= b || x4 <= a || x4 >= b
        || y1 <= c || y1 >= d || y2 <= c || y2 >= d || y3 <= c || y3 > d || y4 <= c || y4 > d) {
        printf("8. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
}

inline void get_right_area(int ii, int jj, double &x1, double &y1, double &x2, double &y2, double &x3, double &y3,
                    double &x4, double &y4, double hy, double hx_add, double hy_add, double a, double b, double c,
                    double d) {
    // p1 (x_{NX-1/2}, y_{j-1/2})
    x1 = b - hx_add / 2.;
    y1 = c + jj * hy - hy_add / 2.;
    // p2 (B, y_{j-1/2})
    x2 = b;
    y2 = c + jj * hy - hy_add / 2.;
    // p3 (B, y_{j+1/2})
    x3 = b;
    y3 = c + jj * hy + hy_add / 2.;
    // p4 (x_{NX-1/2}, y_{j+1/2})
    x4 = b - hx_add / 2.;
    y4 = c + jj * hy + hy_add / 2.;
    if (x1 <= a || x1 >= b || x2 <= a || x2 > b || x3 <= a || x3 > b || x4 <= a || x4 >= b
        || y1 <= c || y1 >= d || y2 <= c || y2 >= d || y3 <= c || y3 >= d || y4 <= c || y4 >= d) {
        printf("7. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
}

inline void get_bottom_area(int ii, int jj, double &x1, double &y1, double &x2, double &y2, double &x3, double &y3, double &x4,
                     double &y4, double hx, double hx_add, double hy_add, double a, double b, double c,
                     double d) {
    // p1 (x_{i-1/2}, C)
    x1 = a + ii * hx - hx_add / 2.;
    y1 = c;
    // p2 (x_{i+1/2}, C)
    x2 = a + ii * hx + hx_add / 2.;
    y2 = c;
    // p3 (x_{i+1/2}, y_{1/2})
    x3 = a + ii * hx + hx_add / 2.;
    y3 = c + hy_add / 2.;
    // p4 (x_{i-1/2}, y_{1/2})
    x4 = a + ii * hx - hx_add / 2.;
    y4 = c + hy_add / 2.;
    if (x1 <= a || x1 >= b || x2 <= a || x2 >= b || x3 <= a || x3 >= b || x4 <= a || x4 >= b
        || y1 < c || y1 >= d || y2 < c || y2 >= d || y3 <= c || y3 >= d || y4 <= c || y4 >= d) {
        printf("6. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
}

inline void get_point_11_area(int ii, int jj, double &x1, double &y1, double &x2, double &y2, double &x3, double &y3,
                       double &x4, double &y4, double hx_add, double hy_add, double a, double b, double c,
                       double d) {
    // p1 (x_{NX-1/2}, y_{NY-1/2})
    x1 = b - hx_add / 2.;
    y1 = d - hy_add / 2.;
    // p2 (B, y_{NY-1/2})
    x2 = b;
    y2 = d - hy_add / 2.;
    // p3 (B, D)
    x3 = b;
    y3 = d;
    // p4 (x_{NX-1/2}, D)
    x4 = b - hx_add / 2.;
    y4 = d;
    if (x1 <= a || x1 >= b || x2 <= a || x2 > b || x3 <= a || x3 > b || x4 <= a || x4 >= b
        || y1 <= c || y1 >= d || y2 <= c || y2 >= d || y3 <= c || y3 > d || y4 <= c || y4 > d) {
        printf("5. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
}

inline void get_point_00_area(int ii, int jj, double &x1, double &y1, double &x2, double &y2, double &x3, double &y3,
                       double &x4, double &y4, double hx_add, double hy_add, double a, double b, double c,
                       double d) {
    // p1 (A, C)
    x1 = a;
    y1 = c;
    // p2 (x_{1/2}, C)
    x2 = a + hx_add / 2.;
    y2 = c;
    // p3 (x_{1/2}, y_{1/2})
    x3 = a + hx_add / 2.;
    y3 = c + hy_add / 2.;
    // p4 (A, y_{1/2})
    x4 = a;
    y4 = c + hy_add / 2.;
    if (x1 < a || x1 >= b || x2 <= a || x2 >= b || x3 <= a || x3 >= b || x4 < a || x4 >= b
        || y1 < c || y1 >= d || y2 < c || y2 >= d || y3 <= c || y3 >= d || y4 <= c || y4 >= d) {
        printf("4. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
}

inline void get_point_01_area(int ii, int jj, double &x1, double &y1, double &x2, double &y2, double &x3, double &y3,
                       double &x4, double &y4, double hx_add, double hy_add, double a, double b, double c,
                       double d) {
    // p1 (A, y_{NY-1/2})
    x1 = a;
    y1 = d - hy_add / 2.;
    // p2 (x_{1/2}, y_{NY-1/2})
    x2 = a + hx_add / 2.;
    y2 = d - hy_add / 2.;
    // p3 (x_{1/2}, D)
    x3 = a + hx_add / 2.;
    y3 = d;
    // p4 (A, D)
    x4 = a;
    y4 = d;
    if (x1 < a || x1 >= b || x2 <= a || x2 >= b || x3 <= a || x3 >= b || x4 < a || x4 >= b
        || y1 <= c || y1 >= d || y2 <= c || y2 >= d || y3 <= c || y3 > d || y4 <= c || y4 > d) {
        printf("3. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
}

inline void get_point_10_area(int ii, int jj, double &x1, double &y1, double &x2, double &y2, double &x3, double &y3,
                       double &x4, double &y4, double hx_add, double hy_add, double a, double b, double c,
                       double d) {
    // p1 (x_{NX-1/2}, C)
    x1 = b - hx_add / 2.;
    y1 = c;
    // p2 (B, C)
    x2 = b;
    y2 = c;
    // p3 (B, y_{1/2})
    x3 = b;
    y3 = c + hy_add / 2.;
    // p4 (x_{NX-1/2}, y_{1/2})
    x4 = b - hx_add / 2.;
    y4 = c + hy_add / 2.;
    if (x1 <= a || x1 >= b || x2 <= a || x2 > b || x3 <= a || x3 > b || x4 <= a || x4 >= b
        || y1 < c || y1 >= d || y2 < c || y2 >= d || y3 <= c || y3 >= d || y4 <= c || y4 >= d)
        printf("2. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
}

inline void get_inner_area(int ii, int jj, double &x1, double &y1, double &x2, double &y2, double &x3, double &y3, double &x4,
                    double &y4, double hx, double hy, double hx_add, double hy_add, double a, double b, double c,
                    double d) {
    // p1 (x_{i-1/2}, y_{j-1/2})
    x1 = a + ii * hx - hx_add / 2.; // a + ii * hx_small - hx_lev / 2.
    y1 = c + jj * hy - hy_add / 2.;
    // p2 (x_{i+1/2}, y_{j-1/2})
    x2 = a + ii * hx + hx_add / 2.;
    y2 = c + jj * hy - hy_add / 2.;
    // p3 (x_{i+1/2}, y_{j+1/2})
    x3 = a + ii * hx + hx_add / 2.;
    y3 = c + jj * hy + hy_add / 2.;
    // p4 (x_{i-1/2}, y_{j+1/2})
    x4 = a + ii * hx - hx_add / 2.;
    y4 = c + jj * hy + hy_add / 2.;
    if (x1 <= a || x1 >= b || x2 <= a || x2 >= b || x3 <= a || x3 >= b || x4 <= a || x4 >= b
        || y1 <= c || y1 >= d || y2 <= c || y2 >= d || y3 <= c || y3 >= d || y4 <= c || y4 >= d) {
        printf("1. Inner point, ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
}

inline void get_inner_area(double x, double y, double &x1, double &y1, double &x2, double &y2, double &x3, double &y3, double &x4,
                    double &y4, double hx_add, double hy_add, double a, double b, double c,
                    double d) {
    // p1 (x_{i-1/2}, y_{j-1/2})
    x1 = x - hx_add / 2.;
    y1 = y - hy_add / 2.;
    // p2 (x_{i+1/2}, y_{j-1/2})
    x2 = x + hx_add / 2.;
    y2 = y - hy_add / 2.;
    // p3 (x_{i+1/2}, y_{j+1/2})
    x3 = x + hx_add / 2.;
    y3 = y + hy_add / 2.;
    // p4 (x_{i-1/2}, y_{j+1/2})
    x4 = x - hx_add / 2.;
    y4 = y + hy_add / 2.;
    if (x1 <= a || x1 >= b || x2 <= a || x2 >= b || x3 <= a || x3 >= b || x4 <= a || x4 >= b
        || y1 <= c || y1 >= d || y2 <= c || y2 >= d || y3 <= c || y3 >= d || y4 <= c || y4 >= d) {
        printf("1. Inner point, ERROR INDEX x=%.8le y=%.8le : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", x, y, x1, y1, x2, y2, x3, y3, x4, y4);
    }
}

double *solve_1(int *grid, int *gridPr);

double *calc_error_1(int *grid, double *solution, double tt, int nx3_1, int ny3_1, double hx_smallest, double hy_smallest, int max_lvl);

double *calc_exact_1(int *grid, double t, int nx3_1, int ny3_1, double hx_smallest, double hy_smallest, int max_lvl);

GraphDouble* solve_2();

GraphDouble* calc_error_2(const GraphInt& grid, const GraphInt& solution, double tt, int nx3_1, int ny3_1, double hx_smallest, double hy_smallest, int max_lvl);

GraphDouble* calc_exact_2(const GraphInt& grid, double t, int nx3_1, int ny3_1, double hx_smallest, double hy_smallest, int max_lvl);

#endif //FEM_CIRCLE_COMMON_H