#include "consts.h"

#if defined _WIN32 || defined __CYGWIN__ || defined WIN32
long timerStart = 0;
#else
//struct timeval timerStart;
#endif

double A = 0.;
double B = 0.;
double C = 0.;
double D = 0.;
unsigned int NX = 0;
unsigned int NY = 0;
unsigned int NX_1 = 0;
unsigned int NY_1 = 0;
int NX3 = 0;
int NY3 = 0;
int NX3_1 = 0;
int NY3_1 = 0;
int XY = 0;
double TAU = 0.;
int TIME_STEP_CNT = 0;
int JAK_ITER_CNT = 0;
double HX = 0.;
double HY = 0.;
double R_SQ = 0.;
double INN_DENSITY = 0.;
double OUT_DENSITY = 0.;
double U = 0.;
double V = 0.;
double DBL_MIN_TRIM = 1.e-16;
double RES_EPS = 1.e-14;
double EPS = 1.e-8;
int INTEGR_TYPE = 1;
unsigned int IDEAL_SQ_SIZE_X = 0;
unsigned int IDEAL_SQ_SIZE_Y = 0;
double CENTER_OFFSET_X = 0.;
double CENTER_OFFSET_Y = 0.;
double OMEGA = 0.;
int* G1 = nullptr;
int* G2 = nullptr;
int* G3 = nullptr;
int* G4 = nullptr;
int CP00 = 0;
int CP10 = 0;
int CP11 = 0;
int CP01 = 0;
int R_LVL = 0;
int R = 0;
double EPS_GRID = 0.;
double HX_SMALLEST = 0.;
double HY_SMALLEST = 0.;
int APPROX_TYPE  = 1;