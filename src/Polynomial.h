/*
  File:             Polynomial.h
  Created by:       Oleksii Pokotylo
  First published:  07.05.2014
  Last revised:     07.05.2014

  Contains the polynomial classifier the DD-plot classification.

  For a description of the algorithm, see:
    Li, J., Cuesta-Albertos, J. A. and Liu, R. Y. (2012). DD-classifier: Nonparametric classification procedure based on
DD-plot, Journal of the American Statistical Association 107(498): 737 - 753.
*/

TPoint GetPolynomial(unsigned degree, TMatrix& points);

TPoint PolynomialLearnCV(TMatrix& input, unsigned numClass1, unsigned int maxDegree, unsigned int chunkNumber, int *degree, int *axis);
