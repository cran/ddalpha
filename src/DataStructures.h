/*
  File:             DataStructures.cpp
  Created by:       Pavlo Mozharovskyi
  First published:  28.02.2013
  Last revised:     28.02.2013

  Defines the data structures used in the package.
*/

typedef vector<double> TPoint;
typedef vector<vector<double> > TMatrix;
typedef vector<int> TVariables;

struct UPoint{
	 int pattern;
	 double value;
	 UPoint(int pattern = 0, double value = 0){
		 this->pattern = pattern;
		 this->value = value;
	 }
 };

 struct Feature{
	 unsigned int order;
	 int number;
	 double angle;
	 unsigned int error;
	 Feature(unsigned int order = 0, int number = 0, double angle = 0, unsigned int error = INT_MAX){
		 this->order = order;
		 this->number = number;
		 this->angle = angle;
		 this->error = error;
	 }
 };

struct OrderRec {
	int order;
	double value;
	OrderRec(int order = -1, double value = 0) {
		this->order = order;
		this->value = value;
	}
};

struct SortRec { 
	double v; 
	TPoint* p; 
	SortRec(double v = 0, TPoint* p = NULL) {
		this->v = v;
		this->p = p;
	}
};

struct SortRecClassified {
	int c; 
	double v; 
	SortRecClassified(int c = 0, double v = 0, int p = 0) {
		this->c = c;
		this->v = v;
	}
};

typedef vector<Feature> Features;
typedef vector<UPoint> UPoints;
