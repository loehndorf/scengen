// Header file for cubic.c

#ifndef HKW_CUBIC_H
#define HKW_CUBIC_H

#include "dll_export_def.h"

#define N 4           /// number of moments
#define NBITMAX 1000
#define EPSILON 1E-12 /// desired precision on the gradient's infinite norm
#define X_DEP 1E-3
#define START_DEV 0.5 /// max deviation from the starting point (0,1,0,0)
#define IDENT1 0
#define IDENT2 1
#define IDENT3 0
#define IDENT4 0

/// the main routine that finds the coefficients and puts them to xk
DLL_PUBLIC double cubic_solve(double *xk);

// input data
DLL_PUBLIC double InMom[13]; /// Input moments: 12 moments + InMom[0]:=1
DLL_PUBLIC double TgMom[4];  /// Target moments

#endif
