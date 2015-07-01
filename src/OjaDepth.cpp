#include "stdafx.h"

void OjaDepthsEx(TDMatrix X, TDMatrix x, int d, int n, int nx,
	double *depths){

	int* counters = new int[d + 1];
	bMatrix A (d + 1, d + 1);
	unsigned long long div0 = choose(n, d);

	TDMatrix covXtemp = cov(X,n,d);
	bMatrix covX(d, d); 
	for (unsigned k = 0; k < d; k++)
	for (unsigned j = 0; j < d; j++)
		covX(k,j) = covXtemp[k][j];
	deleteM(covXtemp);
	double S = pow(abs(determinant(covX)),-0.5);

	for (int obs = 0; obs < nx; obs++){
		long double sumVolume = 0;
		unsigned long long numSimplicesChecked = 0;

		int p = d - 1;
		for (int i = 0; i < p; i++){ counters[i] = i; }counters[p] = p - 1;
		while (counters[0] != n - (p + 1)){
			int i = p;
			while (i > 0 && counters[i] == n - (p + 1) + i){ i--; }
			counters[i]++; int j = i + 1;
			while (j < p + 1){ counters[j] = counters[j - 1] + 1; j++; }
			for (int j = 0; j < d; j++){
				for (int k = 0; k < d; k++){
					A(j + 1, k) = X[counters[k]][j];
				}
			}
			for (int j = 0; j < d; j++){
				A(j + 1, d) = x[obs][j];
			}
			for (int k = 0; k < d + 1; k++){
				A(0,k) = 1;
			}
			double volume = abs(determinant(A));
			sumVolume += volume;
			numSimplicesChecked ++;
		}
		bool sc = numSimplicesChecked == div0;
		double O = sumVolume / fact(d) / div0;

		double depth = 1/(1+O*S);
		depths[obs] = depth;
	}
	
	delete[] counters;
}

void OjaDepthsApx(TDMatrix X, TDMatrix x, int d, int n, int nx,
	unsigned long long k, double *depths){

	int* counters = new int[d + 1];
	bMatrix A(d + 1, d + 1);
	TDMatrix covXtemp = cov(X, n, d);
	bMatrix covX(d, d);
	for (unsigned k = 0; k < d; k++)
	for (unsigned j = 0; j < d; j++)
		covX(k, j) = covXtemp[k][j];
	deleteM(covXtemp);
	double S = pow(abs(determinant(covX)), -0.5);

	for (int obs = 0; obs < nx; obs++){

		long double sumVolume = 0;

		for (unsigned long long i = 0; i < k; i++){
			// Generate a combination of indices
			for (int j = 0; j < d; j++){
				bool _new = false;
				do{
					_new = true;
					counters[j] = random(n);
					for (int l = 0; l < j; l++){
						if (counters[l] == counters[j]){
							_new = false;
							break;
						}
					}
				} while (!_new);
			}
			// Construct the simplex out of it
			for (int j = 0; j < d; j++){
				for (int k = 0; k < d; k++){
					A(j + 1, k) = X[counters[k]][j];
				}
			}
			for (int j = 0; j < d; j++){
				A(j + 1, d) = x[obs][j];
			}
			for (int k = 0; k < d + 1; k++){
				A(0, k) = 1;
			}
			double volume = abs(determinant(A));
			sumVolume += volume;
		}
		double O = sumVolume / fact(d) / k;

		double depth = 1 / (1 + O*S);
		depths[obs] = depth;
	}
	delete[] counters;
}
