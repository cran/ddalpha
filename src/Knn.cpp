/*
  File:             Knn.cpp
  Created by:       Pavlo Mozharovskyi
  First published:  28.02.2013
  Last revised:     28.02.2013

  The realization of the fast binary affine-invariante KNN classifier.
*/

#include "stdafx.h"

#define DISTTYPE_EUCLIDEAN    1
#define DISTTYPE_MAXIMUM      2
#define DISTTYPE_STANDARDIZE 64

static TMatrix Sigma;

static int CompareValue(UPoint a, UPoint b)
/* This routine is passed to the sort routine. */
{
  return (a.value < b.value);
}

static int GetMean(TMatrix x, TPoint *mean){
	unsigned int n = x.size();if (n <= 0){return -1;}
	unsigned int d = x[0].size();if (d <= 0){return -1;}
	mean->resize(d);
	for (unsigned int i = 0; i < n; i++){
		for (unsigned int j = 0; j < d; j++){
			(*mean)[j] += x[i][j];
		}
	}
	for (unsigned int j = 0; j < d; j++){
		(*mean)[j] /= (double)n;
	}
	return 0;
}

static int GetCov(TMatrix x, TMatrix *cov){
	unsigned int n = x.size();if (n <= 0){return -1;}
	unsigned int d = x[0].size();if (d <= 0){return -1;}
	TPoint mean;GetMean(x, &mean);
	cov->resize(d);
	for (unsigned int i = 0; i < d; i++){(*cov)[i].resize(d);}
	for (unsigned int i = 0; i < n; i++){
		for (unsigned int j = 0; j < d; j++){
			for (unsigned int k = 0; k < d; k++){
				(*cov)[j][k] += (x[i][j] - mean[j])*(x[i][k] - mean[k]);
			}
		}
	}
	for (unsigned int j = 0; j < d; j++){
		for (unsigned int k = 0; k < d; k++){
			(*cov)[j][k] /= (double)(n - 1);
		}
	}
	return 0;
}

static int GetInverted(TMatrix x, TMatrix *inv){
	unsigned int d = x.size();if (d <= 0){return -1;}
	unsigned int _d = x[0].size();if (_d != d){return -1;}
 
#ifndef _MSC_VER

//  NumericMatrix A(d,d);
  
	boost::numeric::ublas::matrix<double> A(d, d);A.clear();
	for (unsigned int i = 0; i < d; i++){
		for (unsigned int j = 0; j < d; j++){
			A(i,j) = x[i][j];
		}
	}
	typedef boost::numeric::ublas::permutation_matrix<std::size_t> pmatrix;
	boost::numeric::ublas::matrix<double> Inv(d, d);Inv.clear();
	pmatrix pm(A.size1());
	int res = lu_factorize(A, pm);
	if (res != 0){return -1;}
	Inv.assign(boost::numeric::ublas::identity_matrix<double> (A.size1()));
	boost::numeric::ublas::lu_substitute(A, pm, Inv);
	
 /* Function solve("solve");
  NumericMatrix Inv = solve(A);*/
  
  inv->resize(d);for(unsigned int i = 0; i < d; i++){(*inv)[i].resize(d);}
  for (unsigned int i = 0; i < d; i++){
		for (unsigned int j = 0; j < d; j++){
			 (*inv)[i][j] = Inv(i,j);
		}
	}
	return 0;
#else
	return 1;
#endif // _MSC_VER
}

static double GetNormalized(TPoint dif){
	unsigned int d = dif.size();
	TPoint tmp(d);
	for (unsigned int i = 0; i < d; i++){
		for (unsigned int j = 0; j < d; j++){
			tmp[i] += dif[j]*Sigma[j][i];
		}
	}
	double res = 0;
	for (unsigned int i = 0; i < d; i++){
		res += tmp[i]*dif[i];
	}
	return res;
}

static double GetDistance(TPoint x, TPoint y, int d, int distType){
	double dist = 0.;
	if ((distType & DISTTYPE_EUCLIDEAN) == DISTTYPE_EUCLIDEAN){
		TPoint distVector(d);
		for (int i = 0; i < d; i++){distVector[i] = x[i] - y[i];}
		if ((distType & DISTTYPE_STANDARDIZE) == DISTTYPE_STANDARDIZE){
			dist = GetNormalized(distVector);
		}else{
			for (int i = 0; i < d; i++){dist += std::pow(x[i] - y[i],2);}
		}
	}
	if ((distType & DISTTYPE_MAXIMUM) == DISTTYPE_MAXIMUM){
		for (int i = 0; i < d; i++){
			if (abs(x[i] - y[i]) > dist){dist = abs(x[i] - y[i]);}
		}
	}
	return dist;
}

static int GetDistances(TMatrix x, TMatrix *dist, int distanceType){
	unsigned int n = x.size();if (n <= 0){return -1;}
	unsigned int d = x[0].size();if (d <= 0){return -1;}
	TMatrix cov;GetCov(x, &cov);GetInverted(cov, &Sigma);
	dist->resize(n);
	for (unsigned int i = 0; i < n; i++){(*dist)[i].resize(n);}
	for (unsigned int i = 0; i < n - 1; i++){
		for (unsigned int j = i + 1; j < n; j++){
			(*dist)[i][j] = (*dist)[j][i] = 
				GetDistance(x[i], x[j], d, distanceType);
		}
	}
	return 0;
}

static int GetMaxIndex(TVariables v){
	int index = 0;int maxValue = v[0];int d = v.size();
	for (int i = 1; i < d; i++){
		if (v[i] > maxValue){index = i; maxValue = v[i];}
	}
	return index;
}

int KnnCv(TMatrix points, TVariables labels, unsigned int kMax, int distType, unsigned int numFolds){
	// Collect basic statistics (and check it)
	unsigned int d = points[0].size();unsigned int n = points.size();
	int q = labels[GetMaxIndex(labels)] + 1;
	if (labels.size() != n){return -1;}
	// Prepare indicator array for Jack-Knifing
	TMatrix dist;GetDistances(points, &dist, distType);
	vector<vector<UPoint> > indicators;indicators.resize(n);
	for (unsigned i = 0; i < n; i++){
		indicators[i].resize(n);
		for (unsigned j = 0; j < n; j++){
			indicators[i][j] = UPoint(labels[j], dist[i][j]);
		}
	}
	for (unsigned i = 0; i < n; i++){indicators[i][i].value = -1;}
	for (unsigned i = 0; i < n; i++){
		sort(indicators[i].begin(), indicators[i].end(), CompareValue);
	}
	// Jack-knifing
	vector<TVariables> decisions(kMax + 1);
	for (unsigned int i = 0; i < kMax + 1; i++){decisions[i].resize(n);}
	for (unsigned int i = 0; i < n; i++){decisions[0][i] = labels[i];}
	for (unsigned int i = 0; i < n; i++){
		TVariables locVotes(q);
		for (unsigned j = 1; j < kMax + 1; j++){
			locVotes[indicators[i][j].pattern]++;
			decisions[j][i] = GetMaxIndex(locVotes);
		}
	}
	TVariables guessed(kMax + 1);
	for (unsigned i = 1; i < kMax + 1; i++){
		for (unsigned j = 0; j < n; j++){
			if (decisions[i][j] == decisions[0][j]){guessed[i]++;}
		}
	}
	return GetMaxIndex(guessed);
}

int Knn(TMatrix objects, TMatrix points, TVariables labels, unsigned int k, 
		int distType, TVariables *output){
	unsigned int n = points.size();if (n <= 0){return -1;}
	unsigned int d = points[0].size();if (d <= 0){return -1;}
	int q = labels[GetMaxIndex(labels)] + 1;
	unsigned int nobjects = objects.size();if (nobjects <= 0){return -1;}
	if (labels.size() != n){return -1;}if (objects[0].size() != d){return -1;}
	output->resize(nobjects);
	if ((distType & DISTTYPE_STANDARDIZE) == DISTTYPE_STANDARDIZE){
		TMatrix cov;GetCov(points, &cov);GetInverted(cov, &Sigma);
	}
	for (unsigned int i = 0; i < nobjects; i++){
		vector<UPoint> indicators(n);
		for (unsigned int j = 0; j < n; j++){
			indicators[j] = UPoint(labels[j], 
				GetDistance(objects[i], points[j], d, distType));
		}
		sort(indicators.begin(), indicators.end(), CompareValue);
		TVariables locVotes(q);
		for (unsigned j = 0; j < k; j++){
			locVotes[indicators[j].pattern]++;
		}
		(*output)[i] = GetMaxIndex(locVotes);
	}
	return 0;
}

int GetK_JK_Binary(TMatrix points, TVariables cardinalities, unsigned int maxk){
	// Collect basic statistics (and check it)
	unsigned int d = points[0].size();
	unsigned int n = points.size();
	unsigned int q = cardinalities.size();    if (q != 2){return -1;}
	// Prepare indicator array for Jack-Knife
	TMatrix dist;GetDistances(points, &dist, DISTTYPE_EUCLIDEAN | DISTTYPE_STANDARDIZE);
	vector<vector<UPoint> > indicators;
	indicators.resize(n);
	for (unsigned int i = 0; i < n; i++){
		indicators[i].resize(n);
		for (int j = 0; j < cardinalities[0]; j++){indicators[i][j] = UPoint(0, dist[i][j]);}
		for (unsigned j = cardinalities[0]; j < n; j++){indicators[i][j] = UPoint(1, dist[i][j]);}
	}
	for (unsigned int i = 0; i < n; i++){indicators[i][i].value = -1;}
	for (unsigned int i = 0; i < n; i++){sort(indicators[i].begin(), indicators[i].end(), CompareValue);}
	// Jack-knifing
	vector<TVariables> decisions(maxk);
	decisions[0].resize(n);for (unsigned int j = 0; j < n; j++){decisions[0][j] = indicators[j][1].pattern;}
	for (unsigned int i = 1; i < maxk; i++){
		decisions[i].resize(n);
		for (unsigned int j = 0; j < n; j++){
			decisions[i][j] = decisions[i - 1][j] + indicators[j][i + 1].pattern;
		}
	}
	for (unsigned i = 0; i < maxk; i++){
		for (unsigned j = 0; j < n; j++){
			decisions[i][j] = (decisions[i][j] > (i + 1)/2 ? 1 : 0);
		}
	}
	TVariables errors(maxk);
	for (unsigned int i = 0; i < maxk; i++){
		for (int j = 0; j < cardinalities[0]; j++){errors[i] += decisions[i][j];}
		for (unsigned j = cardinalities[0]; j < n; j++){errors[i] += 1 - decisions[i][j];}
	}
	int k = -1; int minErr = n + 1;
	for (unsigned int i = 0; i < maxk; i++){if (errors[i] < minErr){k = i + 1; minErr = errors[i];}}
	return k;
}

int Knn_Classify_Binary(TMatrix objects, TMatrix points, TVariables cardinalities, unsigned int k, TVariables *output){
	unsigned int n = points.size();if (n <= 0){return -1;}
	unsigned int d = points[0].size();if (d <= 0){return -1;}
	unsigned int nobjects = objects.size();if (nobjects <= 0){return -1;}
	if (objects[0].size() != d){return -1;}
	output->resize(nobjects);
	TMatrix cov;GetCov(points, &cov);
	GetInverted(cov, &Sigma);
	for (unsigned int i = 0; i < nobjects; i++){
		TPoint point = objects[i];
		TPoint tmp(d);
		TPoint dist(n);
		for (unsigned int j = 0; j < n; j++){
			for (unsigned int l = 0; l < d; l++){tmp[l] = point[l] - points[j][l];}
			dist[j] = GetNormalized(tmp);
		}
		vector<UPoint> indicators(n);
		for (int j = 0; j < cardinalities[0]; j++){indicators[j] = UPoint(0, dist[j]);}
		for (unsigned j = cardinalities[0]; j < n; j++){indicators[j] = UPoint(1, dist[j]);}
		sort(indicators.begin(), indicators.end(), CompareValue);
		unsigned decision = 0;
		for(unsigned int j = 0; j < k; j++){decision += indicators[j].pattern;}
		(*output)[i] = (decision > k/2 ? 1 : 0);
	}
	return 0;
}
