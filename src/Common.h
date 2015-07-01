/*
  File:             Common.h
  Created by:       Pavlo Mozharovskyi
  First published:  17.05.2013
  Last revised:     17.05.2013
  
  Commonly used functions.
*/

#pragma once

const double eps_pivot = 1e-10;

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

unsigned long long choose(unsigned long long n, unsigned long long k);
unsigned long long fact(unsigned long long n);
bool solveUnique(TDMatrix A, double* b, double* x, int d);

double determinant(bMatrix& m);
TDMatrix cov(TDMatrix X, int n, int d);

void GetDirections(TMatrix *directions, unsigned int k, unsigned int d);
void GetProjections(TMatrix& points, TMatrix& directions, TMatrix *projections);
