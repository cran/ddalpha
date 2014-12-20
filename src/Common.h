/*
  File:             Common.h
  Created by:       Pavlo Mozharovskyi
  First published:  17.05.2013
  Last revised:     17.05.2013
  
  Commonly used functions.
*/

void GetDirections(TMatrix *directions, unsigned int k, unsigned int d);
void GetProjections(TMatrix& points, TMatrix& directions, TMatrix *projections);
