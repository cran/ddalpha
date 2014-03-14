/*
  File:             Knn.h
  Created by:       Pavlo Mozharovskyi
  First published:  28.02.2013
  Last revised:     28.02.2013
  
  The realization of the KNN classifier.
*/

int GetK_JK_Binary(TMatrix points, TVariables cardinalities, unsigned int maxk);
int Knn_ClassifyOne_Binary(TPoint point, TMatrix points, 
						   TVariables cardinalities, unsigned int k);
int Knn_Classify_Binary(TMatrix objects, TMatrix points, 
						TVariables cardinalities, unsigned int k, TVariables *output);
int GetDDK_JK_Binary(TMatrix points, TVariables cardinalities, unsigned int maxk);
int DDKnn_ClassifyOne(TPoint point, TMatrix points, TVariables cardinalities, 
					  int k);
int KnnCv(TMatrix points, TVariables labels, unsigned int kMax, int distType, 
		  unsigned int numFolds);
int Knn(TMatrix objects, TMatrix points, TVariables labels, unsigned int k, 
		int distType, TVariables *output);
