#include "stdafx.h"

static int CompareAsc(OrderRec x, OrderRec y)
{
	return (x.value < y.value);
}

static int CompareDec(OrderRec x, OrderRec y)
{
	return (x.value > y.value);
}

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

void GetProjections(TMatrix points, TMatrix directions, TMatrix *projections){
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

void GetPrjDepths(TPoint projection, TVariables cardinalities, int curClass, TVariables *prjDepths){
	//Collecting basic statistics
	int n = projection.size();
	int beginIndex = 0;
	for (int i = 0; i < cardinalities.size(); i++){
		if (i >= curClass){break;}
		beginIndex += cardinalities[i];
	}
	int endIndex = beginIndex + cardinalities[curClass] - 1;

	//Preparing structures
	vector<OrderRec> prjSort(n);
	for (int i = 0; i < n; i++){
		prjSort[i].order = i;
		prjSort[i].value = projection[i];
	}
	
	//Calculating projection depths
	TVariables depthsForwards(n);
	TVariables depthsBackwards(n);
	//Forwards
	sort(prjSort.begin(), prjSort.end(), CompareAsc);
	int curDepth = 0;
	for (int i = 0; i < n; i++){
		if ((prjSort[i].order >= beginIndex) && (prjSort[i].order <= endIndex)){curDepth++;}
		depthsForwards[prjSort[i].order] = curDepth;
	}
	//Backwards
	sort(prjSort.begin(), prjSort.end(), CompareDec);
	curDepth = 0;
	for (int i = 0; i < n; i++){
		if ((prjSort[i].order >= beginIndex) && (prjSort[i].order <= endIndex)){curDepth++;}
		depthsBackwards[prjSort[i].order] = curDepth;
	}
	//Merge
	prjDepths->resize(n);
	for (int i = 0; i < n; i++){
		if (depthsForwards[i] < depthsBackwards[i]){
			(*prjDepths)[i] = depthsForwards[i];
		}else{
			(*prjDepths)[i] = depthsBackwards[i];
		}
	}
}

void GetPtPrjDepths(TPoint projection, double point, TVariables cardinalities, TPoint *ptPrjDepths){
	int q = cardinalities.size();
	int n = projection.size();
	ptPrjDepths->resize(q);
	for (int i = 0; i < q; i++){
		int beginIndex = 0;
		for (int j = 0; j < q; j++){
			if (j >= i){break;}
			beginIndex += cardinalities[j];
		}
		int endIndex = beginIndex + cardinalities[i];
		int nPtsBelow = 0;
		int nPtsAbove = 0;
		for (int j = beginIndex; j < endIndex; j++){
			if (projection[j] <= point){nPtsBelow++;}
			if (projection[j] >= point){nPtsAbove++;}
		}
		(*ptPrjDepths)[i] = (nPtsBelow <= nPtsAbove)?(double)nPtsBelow:(double)nPtsAbove;
	}
}

//Indexing from zero
void GetDSpace(TMatrix points, TVariables cardinalities, int k, bool atOnce, TMatrix *dSpace, TMatrix *directions, TMatrix *projections){
	//1. Collecting basic statistics
	int d = points[0].size();
	int n = points.size();
	int q = cardinalities.size();
	if (!atOnce){
		dSpace->resize(n);
		for (int i = 0; i < n; i++){
			GetDepths(points[i], points, cardinalities, k, false, TMatrix(0), TMatrix(0), &(*dSpace)[i]);
		}
		return;
	}
	GetDirections(directions, k, d);
	GetProjections(points, (*directions), projections);
	//2. Calculate projection depths
	vector<vector<TVariables> > prjDepths(k);
	for (int i = 0; i < k; i++){
		vector<TVariables> onePrjDepths(q);
		prjDepths[i] = onePrjDepths;
		for (int j = 0; j < q; j++){
			GetPrjDepths((*projections)[i], cardinalities, j, &prjDepths[i][j]);
		}
	}
	//3. Merge depths
	dSpace->resize(n);
	for (int i = 0; i < n; i++){
		(*dSpace)[i].resize(q);
		for (int j = 0; j < q; j++){
			(*dSpace)[i][j] = cardinalities[j] + 1;
		}
	}
	for (int i = 0; i < k; i++){
		for (int j = 0; j < q; j++){
			for (int l = 0; l < n; l++){
				if (prjDepths[i][j][l] < (*dSpace)[l][j]){
					(*dSpace)[l][j] = prjDepths[i][j][l];
				}
			}
		}
	}
	for (int i = 0; i < q; i++){
		for (int j = 0; j < n; j++){
			(*dSpace)[j][i] /= cardinalities[i];
		}
	}
}

void GetDepths(TPoint point, TMatrix points, TVariables cardinalities, int k, bool atOnce, TMatrix directions, TMatrix projections, TPoint *depths){
	//1. Collecting basic statistics
	int d = points[0].size();
	int n = points.size();
	int q = cardinalities.size();
	int _k = 0;
	TMatrix _directions;
	TMatrix _projections;
	if (!atOnce){
		_k = k;
		GetDirections(&_directions, _k, d);
		GetProjections(points, _directions, &_projections);
	}else{
		_k = directions.size();
		_directions = directions;
		_projections = projections;
	}	
	//2. Calculate projection depths
	TPoint pointProjections(_k);
	for (int i = 0; i < _k; i++){
		double curPrj = 0;
		for (int j = 0; j < d; j++){
			curPrj += point[j]*_directions[i][j];
		}
		pointProjections[i] = curPrj;
	}
	TMatrix ptPrjDepths(_k);
	for (int i = 0; i < _k; i++){
//		TPoint onePrjDepths(q);
//		ptPrjDepths[i] = onePrjDepths;
		GetPtPrjDepths(_projections[i], pointProjections[i], cardinalities, &ptPrjDepths[i]);
	}
	//3. Merge depths
	depths->resize(q);
	for (int i = 0; i < q; i++){
		(*depths)[i] = cardinalities[i] + 1;
	}
	for (int i = 0; i < _k; i++){
		for (int j = 0; j < q; j++){
			if (ptPrjDepths[i][j] < (*depths)[j]){
				(*depths)[j] = ptPrjDepths[i][j];
			}
		}
	}
	for (int i = 0; i < q; i++){
		(*depths)[i] /= (double)cardinalities[i];
	}
}
