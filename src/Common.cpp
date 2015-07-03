/*
  File:             Common.cpp
  Created by:       Pavlo Mozharovskyi
  First published:  17.05.2013
  Last revised:     17.05.2013
  
  Commonly used functions.
*/

#include "stdafx.h"

TDMatrix asMatrix(double* arr, int n, int d){
	TDMatrix mat = new double*[n];
	for (int i = 0; i < n; i++)
		mat[i] = arr + i*d;
	return mat;
}

double** newM(int n, int d){
	double* a = new double[n*d];
	return asMatrix(a, n, d);
}

void deleteM(TDMatrix X){
	delete[] X[0];
	delete[] X;
}

void printMatrix(TDMatrix mat, int n, int d){
#ifdef _MSC_VER
	for (int i = 0; i < n; i++){
		for (int j = 0; j < d; j++)
			cout << mat[i][j] << "\t";
		cout << "\n";
	}
	cout << "\n";
#endif
}


unsigned long long choose(unsigned long long n, unsigned long long k){
	unsigned long long r = n--; unsigned long long d = 2;
	while (d <= k){ r *= n--; r /= d++; }
	return r;
}

unsigned long long fact(unsigned long long n){
	unsigned long long r = 1; unsigned long long i = 2;
	while (i <= n){ r *= i++; }
	return r;
}

/* -------------------------------------------------------------------------- */
/* By Rainer Dyckerhoff, modified by Pavlo Mozharovskyi                       */
/* Solves a uniquely solvable system of linear equations                      */
/* -------------------------------------------------------------------------- */
bool solveUnique(TDMatrix A, double* b, double* x, int d){
	int imax, jmax;
	int* colp = new int[d];
	double amax;
	for (int k = 0; k < d - 1; k++) {
		imax = k;
		amax = abs(A[k][k]);
		colp[k] = k;
		// Spaltenmaximum finden
		for (int i = k + 1; i < d; i++) {
			if (abs(A[i][k]) > amax) {
				amax = abs(A[i][k]);
				imax = i;
			}
		}
		// Spaltenmaximum gleich null => complete pivoting
		if (amax < eps_pivot) {
			for (int j = k + 1; j < d; j++) {
				for (int i = k; i < d; i++) {
					if (abs(A[i][j]) > amax) {
						amax = abs(A[i][j]);
						imax = i;
						jmax = j;
					}
				}
			}
			if (amax < eps_pivot) {
				delete[] colp;
				return false;
			}
			// Spaltentausch
			for (int i = 0; i < d; i++) {
				double tmp = A[i][k];
				A[i][k] = A[i][jmax];
				A[i][jmax] = tmp;
			}
			colp[k] = jmax;
		}
		// Zeilentausch
		if (imax != k) {
			for (int j = k; j < d; j++) {
				double tmp = A[k][j];
				A[k][j] = A[imax][j];
				A[imax][j] = tmp;
			}
			double tmp = b[k];
			b[k] = b[imax];
			b[imax] = tmp;
		}
		// Elimination
		for (int i = k + 1; i < d; i++) {
			double factor = A[i][k] / A[k][k];
			for (int j = k + 1; j < d; j++){
				A[i][j] -= factor * A[k][j];
			}
			b[i] -= factor * b[k];
		}
	}
	// Rücksubstituition
	colp[d - 1] = d - 1;
	for (int k = d - 1; k >= 0; k--) {
		x[k] = b[k] / A[k][k];
		for (int i = k - 1; i >= 0; i--) b[i] -= x[k] * A[i][k];
	}
	// Spaltenvertauschungen rückgängig machen
	for (int k = d - 1; k >= 0; k--) {
		if (colp[k] != k) {
			double temp = x[k];
			x[k] = x[colp[k]];
			x[colp[k]] = temp;
		}
	}
	delete[] colp;
	return true;
}

double determinant(bMatrix& m)
{
	bMatrix mLu(m);
	bPM pivots(m.size1());

	if (bnu::lu_factorize(mLu, pivots))
		return 0;

	double det = 1;
	for (size_t i = 0; i < pivots.size(); ++i)
	{
		if (pivots(i) != i)
			det *= -1;
		det *= mLu(i, i);
	}
	return det;
}

TDMatrix cov(TDMatrix X, int n, int d) {
	double* means = new double[d];
	double* dev = new double[d];
	// zeroing TDMatrix
	TDMatrix covX = newM(d, d);
	for (unsigned k = 0; k < d; k++)
		for (unsigned j = 0; j < d; j++)
			covX[k][j] = 0;
	// means
	for (unsigned i = 0; i < d; i++) {
		means[i] = 0.0;
		for (unsigned j = 0; j < n; j++)
			means[i] += X[j][i];
		means[i] /= n;
	}
	for (unsigned i = 0; i < n; i++) {
		// deviations
		for (unsigned k = 0; k < d; k++) {
			dev[k] = X[i][k] - means[k];
		}
		// add to cov
		for (unsigned k = 0; k < d; k++) {
			for (unsigned j = 0; j < d; j++) {
				covX[k][j] += dev[k] * dev[j];
			}
		}
	}
	//scale
	for (unsigned i = 0; i < d; i++) {
		for (unsigned j = 0; j < d; j++) {
			covX[i][j] /= n - 1;
		}
	}
	delete[] means;
	delete[] dev;
	return covX;
}

void GetDirections(TMatrix *directions, unsigned int k, unsigned int d){
	directions->resize(k);
	for (unsigned int i = 0; i < k; i++){
		TPoint direction(d);
		double sqrSum = 0;
		for (unsigned int j = 0; j < d; j++){
			direction[j] = normDist(rEngine);
			sqrSum += direction[j]*direction[j];
		}
		sqrSum = sqrt(sqrSum);
		for (unsigned int j = 0; j < d; j++){
			direction[j] = direction[j]/sqrSum;
		}
		(*directions)[i] = direction;
	}
}

void GetProjections(TMatrix& points, TMatrix& directions, TMatrix *projections){
	int d = points[0].size();
	int n = points.size();
	int k = directions.size();

	projections->resize(k);
	for (int i = 0; i < k; i++){
		TPoint projection(n);
		for (int j = 0; j < n; j++){
			double sum = 0;
			for (int l = 0; l < d; l++){
				sum += points[j][l]*directions[i][l];
			}
			projection[j] = sum;
		}
		(*projections)[i] = projection;
	}
}
