#include <math.h> 
#include "stdlib.h"
#define MAXSTP 1000000 
#define TINY 1.0e-30 
#define MAXSTEPS 50000

void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1, double hmin, int *nok, int *nbad, void (*derivs)(double, double [], double []), void (*rkqs)(double [], double [], int, double *, double, double, double [], double *, double *, void (*)(double, double [], double [])));

void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps, double yscal[], double *hdid, double *hnext, void (*derivs)(double, double [], double []));



void rkck(double y[], double dydx[], int n, double x, double h, double yout[], double yerr[], void (*derivs)(double, double [], double []));
