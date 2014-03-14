/*
  File:             ProjectionDepth.h
  Created by:       Pavlo Mozharovskyi
  First published:  17.05.2013
  Last revised:     17.05.2013
  
  Computation of the projection depth using random sampling.

  For a description of the method, see:
    Zuo, Y.J. and Serfling, R. (2000). General notions of statistical depth
	  function. Annals of Statistics 28, 461-482.
*/

int GetDepthsPrj(TMatrix points, TMatrix objects, TVariables cardinalities,
				 int k, bool newDirs, TMatrix *depths, TMatrix *directions,
				 TMatrix *projections);
