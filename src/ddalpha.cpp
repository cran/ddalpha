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

void setSeed(unsigned int random_seed){
	if (random_seed != 0) {
		setseed(random_seed);
		rEngine.seed(random_seed);
	}
	else {
		setseed(time(NULL));
		rEngine.seed(time(NULL));
	}
}

void IsInConvexes(double *points, int *dimension, int *cardinalities, int *numClasses, double *objects, int *numObjects, int *seed, int *isInConvexes){
	setSeed(*seed);
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
	TIntMatrix answers(o.size());
	int error = 0;
	InConvexes(x, cars, o, error, &answers);
	for (int i = 0; i < numObjects[0]; i++)
    for (int j = 0; j < numClasses[0]; j++){
  		isInConvexes[numClasses[0]*i+j] = answers[i][j];
  	}
}

void ZDepth(double *points, double *objects, int *numPoints, int *numObjects, int *dimension, int *seed, double *depths){
	setSeed(*seed);
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

void HDepth(double *points, double *objects, int *numObjects, int *dimension, int *cardinalities, int *numClasses, double *directions, double *projections, int *k, int *sameDirs, int *seed, double *depths){
	setSeed(*seed);
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
	TMatrix dirs;
	TMatrix prjs;
	for (int i = 0; i < numObjects[0]; i++){
		TPoint dps;
		GetDepths(z[i], x, cars, k[0], i == 0 ? 0 : sameDirs[0] /*at the first step fill the matrices*/, dirs, prjs, &dps);
		for (int j = 0; j < numClasses[0]; j++){
			depths[i * numClasses[0] + j] = dps[j];
		}
	}
}

void HDSpace(double *points, int *dimension, int *cardinalities, int *numClasses, int *k, int *sameDirs, int *seed, double *dSpace, double *directions, double *projections){
	setSeed(*seed);
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
	for (unsigned i = 0; i < direction.size(); i++){
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
	TPoint direction; unsigned int power;
	LearnCV(x, y, minFeatures[0], upToPower[0], numFolds[0], &direction, &power);
	ray[0] = power;
	for (unsigned i = 0; i < direction.size(); i++){
		ray[i + 1] = direction[i];
	}
}

void AlphaClassify(double *points, int *numPoints, int *dimension, int *degree, double *ray, int *output){
	TMatrix x(numPoints[0]);
	for (int i = 0; i < numPoints[0]; i++){x[i] = TPoint(dimension[0]);}
	for (int i = 0; i < numPoints[0]; i++){
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
	}
	TMatrix _x;
	ExtendWithProducts(x, degree[0], &_x);
	TPoint direction(_x[0].size());
	for (unsigned i = 0; i < _x[0].size(); i++){
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

void KnnLearnJK(double *points, int *labels, int *numPoints, int *dimension, int *kmax, int *distType, int *k){
	TMatrix x(numPoints[0]);TVariables y(numPoints[0]);
	for (int i = 0; i < numPoints[0]; i++){
		x[i] = TPoint(dimension[0]);
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
		y[i] = labels[i];
	}
	k[0] = KnnCv(x, y, kmax[0], distType[0], 0);
}

void KnnClassify(double *objects, int *numObjects, double *points, int *labels, int *numPoints, int *dimension, int *k, int *distType, int *output){
	TMatrix x(numPoints[0]);TVariables y(numPoints[0]);
	for (int i = 0; i < numPoints[0]; i++){
		x[i] = TPoint(dimension[0]);
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
		y[i] = labels[i];
	}
	TMatrix z(numObjects[0]);
	for (int i = 0; i < numObjects[0]; i++){
		z[i] = TPoint(dimension[0]);
		for (int j = 0; j < dimension[0]; j++){
			z[i][j] = objects[i * dimension[0] + j];
		}
	}
	TVariables result;
	Knn(z, x, y, k[0], distType[0], &result);
	for (int i = 0; i < numObjects[0]; i++){
		output[i] = result[i];
	}
}

void PolynomialLearnCV(double *points, int *numPoints, int *dimension, int *cardinalities, int *maxDegree, int *chunkNumber, int *seed, /*OUT*/ int *degree, /*OUT*/ int *axis, /*OUT*/ double *polynomial){
	setSeed(*seed);
	TMatrix x(numPoints[0]);
	for (int i = 0; i < numPoints[0]; i++){ x[i] = TPoint(dimension[0]); }
	for (int i = 0; i < numPoints[0]; i++){
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
	}
	TVariables y(numPoints[0]);
	for (int i = 0; i < cardinalities[0]; i++){ y[i] = 1; }
	for (int i = cardinalities[0]; i < numPoints[0]; i++){ y[i] = -1; }
	
	TPoint pol = PolynomialLearnCV(x, cardinalities[0], *maxDegree, *chunkNumber, degree, axis);

	for (unsigned i = 0; i < pol.size(); i++){
		polynomial[i] = pol[i];
	}
}
/* everything implemented in R
void PolynomialClassify(double *points, int *numPoints, int *dimension, int *degree, double *ray, int *output){
	TMatrix x(numPoints[0]);
	for (int i = 0; i < numPoints[0]; i++){ x[i] = TPoint(dimension[0]); }
	for (int i = 0; i < numPoints[0]; i++){
		for (int j = 0; j < dimension[0]; j++){
			x[i][j] = points[i * dimension[0] + j];
		}
	}
	TMatrix _x;
	ExtendWithProducts(x, degree[0], &_x);
	TPoint direction(_x[0].size());
	for (unsigned i = 0; i < _x[0].size(); i++){
		direction[i] = ray[i + 1];
	}
	TVariables y;
	Classify(_x, direction, &y);
	for (int i = 0; i < numPoints[0]; i++){
		output[i] = y[i];
	}
}
*/

void ProjectionDepth(double *points, double *objects, int *numObjects,
					 int *dimension, int *cardinalities, int *numClasses,
					 double *directions, double *projections, int *k,
					 int *newDirs, int *seed, double *depths){
	setSeed(*seed);
	int numPoints = 0;
	for (int i = 0; i < numClasses[0]; i++){
		numPoints += cardinalities[i];
	}
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
	if (!newDirs[0]){
		dirs.resize(k[0]);prjs.resize(k[0]);
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
	TMatrix _depths;
	GetDepthsPrj(x, z, cars, k[0], newDirs[0], &_depths, &dirs, &prjs);
	for (int i = 0; i < numObjects[0]; i++){
		for (int j = 0; j < numClasses[0]; j++){
			depths[i * numClasses[0] + j] = _depths[i][j];
		}
	}
	if (newDirs[0]){
		for (int i = 0; i < k[0]*dimension[0]; i++){
			directions[i] = dirs[i/dimension[0]][i%dimension[0]];
		}
		for (int i = 0; i < k[0]*numPoints; i++){
			projections[i] = prjs[i/numPoints][i%numPoints];
		}
	}
}

#ifdef __cplusplus
}
#endif
