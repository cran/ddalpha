/*
  File:             ProjectionDepth.cpp
  Created by:       Pavlo Mozharovskyi
  First published:  17.05.2013
  Last revised:     17.05.2013
  
  Computation of the projection depth using random sampling.

  For a description of the method, see:
    Zuo, Y.J. and Serfling, R. (2000). General notions of statistical depth
	  function. Annals of Statistics 28, 461-482.
*/

#include "stdafx.h"

static int CompareAsc(OrderRec x, OrderRec y)
{
	return (x.value < y.value);
}

static int CompareDec(OrderRec x, OrderRec y)
{
	return (x.value > y.value);
}

static void GetMedMad(TPoint &points, double &median, double &mad){
	/* First, determine median */
	int n = points.size();
	sort(points.begin(), points.end());
	median = (points[(n + 1)/2 - 1] + points[(n + 2)/2 - 1])/2.;
	/* Obtain median absolute deviation (from median) (MAD) */
	TPoint deviations(n);
	for (int i = 0; i < n; i++){deviations[i] = abs(points[i] - median);}
	sort(deviations.begin(), deviations.end());
	mad = (deviations[(n + 1)/2 - 1] + deviations[(n + 2)/2 - 1])/2.;
}

void GetPtsPrjDepths(TPoint projection, TPoint objectsProjection,
					 TVariables cardinalities, TMatrix *ptsPrjDepths){
	/* Collect basic statistics */
	int q = cardinalities.size();
	int n = projection.size();
	int m = objectsProjection.size();
	ptsPrjDepths->resize(q);
	for (int i = 0; i < q; i++){
		/* Prepare data and obtain median and mad*/
		int beginIndex = 0;
		for (int j = 0; j < q; j++){
			if (j >= i){break;}
			beginIndex += cardinalities[j];
		}
		int endIndex = beginIndex + cardinalities[i];
		TPoint curClassProjection(projection.begin() + beginIndex,
			projection.begin() + endIndex);
		double median, mad;GetMedMad(curClassProjection, median, mad);
		/* Calculate i-class projectional univariate depths */
		(*ptsPrjDepths)[i] = TPoint(m);
		for (int j = 0; j < m; j++){
			(*ptsPrjDepths)[i][j] = (objectsProjection[j] - median)/mad;
		}
	}
}

int GetDepthsPrj(TMatrix points, TMatrix objects, TVariables cardinalities,
				  int k, bool newDirs, TMatrix *depths, TMatrix *directions,
				  TMatrix *projections){
	/* 1. Collecting basic statistics */
	int d = points[0].size();
	int n = points.size();
	int m = objects.size();
	int q = cardinalities.size();
	TMatrix objectsProjections;
	if (newDirs){
		GetDirections(directions, k, d);
		GetProjections(points, (*directions), projections);
	}
	GetProjections(objects, (*directions), &objectsProjections);
	/* 2. Calculate projection depths */
	vector<vector<vector<double> > > prjDepths(k);
	for (int i = 0; i < k; i++){
		prjDepths[i] = TMatrix(0);
		GetPtsPrjDepths((*projections)[i], objectsProjections[i], cardinalities,
			&prjDepths[i]);
	}
	/* 3. Merge depths */
	depths->resize(m);
	for (int i = 0; i < m; i++){
		(*depths)[i].resize(q);
		for (int j = 0; j < q; j++){
			(*depths)[i][j] = DBL_MIN;
		}
	}
	for (int i = 0; i < k; i++){
		for (int j = 0; j < q; j++){
			for (int l = 0; l < m; l++){
				if (prjDepths[i][j][l] > (*depths)[l][j]){
					(*depths)[l][j] = prjDepths[i][j][l];
				}
			}
		}
	}
	for (int i = 0; i < m; i++){
		for (int j = 0; j < q; j++){		
			(*depths)[i][j] = 1/(1 + (*depths)[i][j]);
		}
	}
	return 0;
}
