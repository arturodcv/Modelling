/* Runge-Kutta-Fehlberg-Simo 78 with adaptive stepsize
 * Implemented by Lluis Alseda
 * Version 2.2, June 26, 2020. */
#include <stdio.h>
#include <stdlib.h>
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define ABS(x)   ((x) < 0.0 ? -(x) : (x))

int RKF78( double *, double *,
           double *, double *,
           double, double, double,
           void *,
           void (*)(double, double, double *, void *));

int RKF78Sys( double *, double *, unsigned,
              double *, double *,
              double, double, double,
              void *,
              void (*)(double, double *, unsigned, double *, void *));

double eighthroot(double);
