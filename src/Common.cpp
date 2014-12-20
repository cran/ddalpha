/*
  File:             Common.cpp
  Created by:       Pavlo Mozharovskyi
  First published:  17.05.2013
  Last revised:     17.05.2013
  
  Commonly used functions.
*/

#include "stdafx.h"

void GetDirections(TMatrix *directions, unsigned int k, unsigned int d){
	boost::random::mt19937_64 rEngine;
	rEngine.seed(time(NULL));
	boost::random::normal_distribution<double> normDist;
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
