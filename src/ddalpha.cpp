/*
  File:             ddalpha.cpp
  Created by:       Pavlo Mozharovskyi
  First published:  28.02.2013
  Last revised:     15.05.2013

  Defines the exported functions for the 'ddalpha'-package.

  For a description of the algorithm, see:
    Lange, T., Mosler, K. and Mozharovskyi, P. (2012). Fast nonparametric classification based on data depth. Statistical Papers.
    Mozharovskyi, P., Mosler, K. and Lange, T. (2013). Classifying real-world data with the DDalpha-procedure. Mimeo.
*/

#include "stdafx.h"

#define EOF (-1)

#ifdef __cplusplus
extern "C" {
#endif

void Sum(double *a, double *b, double *res){
	res[0] = a[0] + b[0];
}

void IsInConvexes(double *points, int *dimension, int *cardinalities, int *numClasses, double *objects, int *numObjects, int *isInConvexes){
	int numPoints = 0;for (int i = 0; i < numClasses[0]; i++){numPoints += cardinalities[i];}
	TMatrix x(numPoints);
	for (int i = 0; i < numPoints; i++){x[i] = TPoint(dimension[0]);}
	for (int i = 0; i < numPoints; i++){
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
	}
	TMatrix o(numObjects[0]);
	for (int i = 0; i < numObjects[0]; i++){o[i] = TPoint(dimension[0]);}
	for (int i = 0; i < numObjects[0]; i++){
		for (int j = 0; j < dimension[0]; j++){
			o[i][j] = objects[i * dimension[0] + j];
		}
	}
	TVariables cars(numClasses[0]);
	for (int i = 0; i < numClasses[0]; i++){
		cars[i] = cardinalities[i];
	}
	TVariables answers;
	int error = 0;
	InConvexes(x, cars, o, error, &answers);
	for (int i = 0; i < numObjects[0]; i++){
		isInConvexes[i] = answers[i];
	}
}

void ZDepth(double *points, double *objects, int *numPoints, int *numObjects, int *dimension, double *depths){
	TMatrix x(numPoints[0]);
	for (int i = 0; i < numPoints[0]; i++){x[i] = TPoint(dimension[0]);}
	for (int i = 0; i < numPoints[0]; i++){
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
	}
	TPoint means;TPoint sds;
	GetMeansSds(x, &means, &sds);
	Standardize(x, means, sds);
	TMatrix z(numObjects[0]);
	for (int i = 0; i < numObjects[0]; i++){z[i] = TPoint(dimension[0]);}
	for (int i = 0; i < numObjects[0]; i++){
		for (int j = 0; j < dimension[0]; j++){
			z[i][j] = objects[i * dimension[0] + j];
		}		
		Standardize(z[i], means, sds);
		int error;
		depths[i] = ZonoidDepth(x, z[i], error);
	}
}

void HDepth(double *points, double *objects, int *numObjects, int *dimension, int *cardinalities, int *numClasses, double *directions, double *projections, int *k, int *sameDirs, double *depths){
	int numPoints = 0;for (int i = 0; i < numClasses[0]; i++){numPoints += cardinalities[i];}
	TMatrix x(numPoints);
	for (int i = 0; i < numPoints; i++){x[i] = TPoint(dimension[0]);}
	for (int i = 0; i < numPoints; i++){
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
	}
	TMatrix z(numObjects[0]);
	for (int i = 0; i < numObjects[0]; i++){z[i] = TPoint(dimension[0]);}
	for (int i = 0; i < numObjects[0]; i++){
		for (int j = 0; j < dimension[0]; j++){
			z[i][j] = objects[i * dimension[0] + j];
		}
	}
	TVariables cars(numClasses[0]);
	for (int i = 0; i < numClasses[0]; i++){
		cars[i] = cardinalities[i];
	}
	TMatrix dirs(0);
	TMatrix prjs(0);
	if (sameDirs[0]){
		dirs.resize(k[0]);
		prjs.resize(k[0]);
		for (int i = 0; i < k[0]; i++){
			dirs[i].resize(dimension[0]);
			for (int j = 0; j < dimension[0]; j++){
				dirs[i][j] = directions[i * dimension[0] + j];
			}
			prjs[i].resize(numPoints);
			for (int j = 0; j < numPoints; j++){
				prjs[i][j] = projections[i * numPoints + j];
			}
		}
	}
	for (int i = 0; i < numObjects[0]; i++){
		TPoint dps;
		GetDepths(z[i], x, cars, k[0], sameDirs[0], dirs, prjs, &dps);
		for (int j = 0; j < numClasses[0]; j++){
			depths[i * numClasses[0] + j] = dps[j];
		}
	}
}

void HDSpace(double *points, int *dimension, int *cardinalities, int *numClasses, int *k, int *sameDirs, double *dSpace, double *directions, double *projections){
	int numPoints = 0;for (int i = 0; i < numClasses[0]; i++){numPoints += cardinalities[i];}
	TMatrix x(numPoints);
	for (int i = 0; i < numPoints; i++){x[i] = TPoint(dimension[0]);}
	for (int i = 0; i < numPoints; i++){
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
	}
	TVariables cars(numClasses[0]);
	for (int i = 0; i < numClasses[0]; i++){
		cars[i] = cardinalities[i];
	}
	TMatrix dsps;
	TMatrix dirs;
	TMatrix prjs;
	GetDSpace(x, cars, k[0], sameDirs[0], &dsps, &dirs, &prjs);
	for (int i = 0; i < numPoints*numClasses[0]; i++){
		dSpace[i] = dsps[i/numClasses[0]][i%numClasses[0]];
	}
	if (sameDirs[0]){
		for (int i = 0; i < k[0]*dimension[0]; i++){
			directions[i] = dirs[i/dimension[0]][i%dimension[0]];
		}
		for (int i = 0; i < k[0]*numPoints; i++){
			projections[i] = prjs[i/numPoints][i%numPoints];
		}
	}
}

void AlphaLearn(double *points, int *numPoints, int *dimension, int *cardinalities, int *degree, int *minFeatures, double *ray){
	TMatrix x(numPoints[0]);
	for (int i = 0; i < numPoints[0]; i++){x[i] = TPoint(dimension[0]);}
	for (int i = 0; i < numPoints[0]; i++){
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
	}
	TVariables y(numPoints[0]);
	for (int i = 0; i < cardinalities[0]; i++){y[i] = 1;}
	for (int i = cardinalities[0]; i < numPoints[0]; i++){y[i] = -1;}
	TMatrix _x;
	ExtendWithProducts(x, degree[0], &_x);
	TPoint direction;
	Learn(_x, y, minFeatures[0], &direction);
	ray[0] = degree[0];
	for (int i = 0; i < direction.size(); i++){
		ray[i + 1] = direction[i];
	}
}

void AlphaLearnCV(double *points, int *numPoints, int *dimension, int *cardinalities, int *upToPower, int *numFolds, int *minFeatures, double *ray){
	TMatrix x(numPoints[0]);
	for (int i = 0; i < numPoints[0]; i++){x[i] = TPoint(dimension[0]);}
	for (int i = 0; i < numPoints[0]; i++){
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
	}
	TVariables y(numPoints[0]);
	for (int i = 0; i < cardinalities[0]; i++){y[i] = 1;}
	for (int i = cardinalities[0]; i < numPoints[0]; i++){y[i] = -1;}
	TPoint direction;int power;
	LearnCV(x, y, minFeatures[0], upToPower[0], numFolds[0], &direction, &power);
	ray[0] = power;
	for (int i = 0; i < direction.size(); i++){
		ray[i + 1] = direction[i];
	}
}

void AlphaClassify(double *points, int *numPoints, int *dimension, int *degree, double *ray, int *output){
	TMatrix x(numPoints[0]);
	for (int i = 0; i < numPoints[0]; i++){x[i] = TPoint(dimension[0]);}
	for (int i = 0; i < numPoints[0]; i++){
		TPoint curPoint();
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
	}
	TMatrix _x;
	ExtendWithProducts(x, degree[0], &_x);
	TPoint direction(_x[0].size());
	for (int i = 0; i < _x[0].size(); i++){
		direction[i] = ray[i + 1];
	}
	TVariables y;
	Classify(_x, direction, &y);
	for (int i = 0; i < numPoints[0]; i++){
		output[i] = y[i];
	}
}

void KnnAffInvLearnJK(double *points, int *dimension, int *cardinalities, int *maxk, int *k){
	int numPoints = cardinalities[0] + cardinalities[1];
	TMatrix x(numPoints);
	for (int i = 0; i < numPoints; i++){x[i] = TPoint(dimension[0]);}
	for (int i = 0; i < numPoints; i++){
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
	}
	TVariables cars(2);cars[0] = cardinalities[0];cars[1] = cardinalities[1];
	k[0] = GetK_JK_Binary(x, cars, maxk[0]);
}

void KnnAffInvClassify(double *objects, int *numObjects, double *points, int *dimension, int *cardinalities, int *k, int *output){
	int numPoints = cardinalities[0] + cardinalities[1];
	TMatrix x(numPoints);
	for (int i = 0; i < numPoints; i++){x[i] = TPoint(dimension[0]);}
	for (int i = 0; i < numPoints; i++){
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
	}
	TVariables cars(2);cars[0] = cardinalities[0];cars[1] = cardinalities[1];
	TMatrix y(numObjects[0]);
	for (int i = 0; i < numObjects[0]; i++){y[i] = TPoint(dimension[0]);}
	for (int i = 0; i < numObjects[0]; i++){
		TPoint curPoint();
		for (int j = 0; j < dimension[0]; j++){
			y[i][j] = objects[i * dimension[0] + j];
		}
	}
	TVariables result;
	Knn_Classify_Binary(y, x, cars, k[0], &result);
	for (int i = 0; i < numObjects[0]; i++){
		output[i] = result[i];
	}
}

#ifdef __cplusplus
}
#endif
