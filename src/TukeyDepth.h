/*
  File:             TukeyDepth.h
  Created by:       Pavlo Mozharovskyi
  First published:  28.02.2013
  Last revised:     28.02.2013

  Computation of the random Tukey data depth.

  For a description of the algorithm, see:
    Cuesta-Albertos, J. A. and Nieto-Reyes, A. (2008). The random Tukey depth. Computational Statistics & Data Analysis 52, 11 (July 2008), 4979-4988.
    Mozharovskyi, P., Mosler, K. and Lange, T. (2013). Classifying real-world data with the DDalpha-procedure. Mimeo.
*/

void GetDSpace(vector<TPoint> points, TVariables cardinalities, int k, bool atOnce, vector<TPoint> *dSpace, vector<TPoint> *directions, vector<TPoint> *projections);
void GetDepths(TPoint point, vector<TPoint> points, TVariables cardinalities, int k, bool atOnce, vector<TPoint> directions, vector<TPoint> projections, TPoint *depths);
