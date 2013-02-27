int GetK_JK_Binary(TMatrix points, TVariables cardinalities, int maxk);
int Knn_ClassifyOne_Binary(TPoint point, TMatrix points, TVariables cardinalities, int k);
int Knn_Classify_Binary(TMatrix objects, TMatrix points, TVariables cardinalities, int k, TVariables *output);
int GetDDK_JK_Binary(TMatrix points, TVariables cardinalities, int maxk);
int DDKnn_ClassifyOne(TPoint point, TMatrix points, TVariables cardinalities, int k);
