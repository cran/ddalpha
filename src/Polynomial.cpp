/*
  File:             Polynomial.cpp
  Created by:       Oleksii Pokotylo
  First published:  07.05.2014
  Last revised:     07.05.2014

  Contains the polynomial classifier the DD-plot classification.

  For a description of the algorithm, see:
    Li, J., Cuesta-Albertos, J. A. and Liu, R. Y. (2012). DD-classifier: Nonparametric classification procedure based on
DD-plot, Journal of the American Statistical Association 107(498): 737 - 753.
*/

#include "stdafx.h"

#include <limits>
#include <boost/math/special_functions/binomial.hpp>

/**
Calculates the empirical risk for two classes on the basis of given depths
and approximates it to get continuous derivative

@param polynomial : Polynomial as a vector of coefficients starting with the first degree(a0 = 0 always)
@param points : nx2 matrix of depths, where each column contains the depths against the corresponding class
@param numClass1:  Number of points belonging to the first class
@return polynomial as a vector of coefficients starting with the first degree (a0 = 0 always)

@throws Smoothed empirical risk
*/
double GetEmpiricalRiskSmoothed(TPoint& polynomial, TMatrix& points, unsigned numClass1){
	const float smoothingConstant = 100;
	unsigned degree = polynomial.size();

	double risk = 0;
	int sign = 1;
	for (unsigned i = 0; i < points.size(); i++){
		if (i >= numClass1)
			sign = -1;

		double val = points[i][0];
		double res = 0;
		for (unsigned j = 0; j < degree; j++){
			res += polynomial[j] * std::pow(val, j);
		}
		risk += 1 / (1 + exp(-smoothingConstant*(points[i][1] - res)*sign));
	}

	return risk / points.size();
}

/**
Calculates the empirical risk for two classes on the basis of given depths

@param polynomial : Polynomial as a vector of coefficients starting with the first degree(a0 = 0 always)
@param points : nx2 matrix of depths, where each column contains the depths against the corresponding class
@param numClass1:  Number of points belonging to the first class
@return polynomial as a vector of coefficients starting with the first degree (a0 = 0 always)

@throws empirical risk
*/
double GetEmpiricalRisk(TPoint& polynomial, TMatrix& points, unsigned numClass1){
	unsigned degree = polynomial.size();

	double risk = 0;
	int sign = 1;
	for (unsigned i = 0; i < points.size(); i++){
		if (i >= numClass1)
			sign = -1;

		double val = points[i][0];
		double res = 0;
		for (unsigned j = 0; j<degree; j++){
			res += polynomial[j] * std::pow(val, j);
		}
		if ((points[i][1] - res)*sign > 0){    // for class1 depths[i,2] > res, for class 2 <
			risk++;
		}
	}

	return risk / points.size();
}

/**
Calculates the coefficients of the polynomial of a given degree
pathing through given points an the origin

@param degree : degree of the polynomial, should be equal the number of points
@param points : degreex2 points for the polynomial to path through

@return polynomial as a vector of coefficients starting with the first degree (a0 = 0 always)

@throws runtime_error in case of singularity
*/
TPoint GetPolynomial(unsigned degree, TMatrix& points) {

	bMatrix A(degree, degree);
	for (unsigned i = 0; i < degree; i++){
		for (unsigned j = 0; j < degree; j++){
			A(i, j) = (std::pow(points[i][0], j + 1));
		}
	}

	bVector b(degree);
	for (unsigned i = 0; i < degree; i++){
		b[i] = points[i][1];
	}

	bPM pm(A.size1());
	boost::numeric::ublas::lu_factorize(A, pm);
	boost::numeric::ublas::lu_substitute(A, pm, b);

	TPoint polynomial(degree);
	for (unsigned i = 0; i < degree; i++){
    if (!(b[i] < std::numeric_limits<double>::max()) || b[i] < -std::numeric_limits<double>::max()){  throw "not a value or inf"; }
		polynomial[i] = b[i];
	}
  
	return polynomial;
}

/**
Chooses the best in classification sense polynomial among
"cardinality" randomly chosen polynomials, passing through
"degree" arbitrary points

@param points:       nx2 matrix of points where first column is an absciss,	n = numClass1 + numClass2
@param numClass1:    Number of points belonging to the first class
@param degree:       Degree of the polynomial
@param n_polynomials:  Number of randomly chosen polynomials

@return polynomial as a vector of coefficients starting with the first degree (a0 = 0 always)
*/
TPoint GetRandomMinPolynomial(TMatrix& points, unsigned numClass1, unsigned degree, unsigned n_polynomials){

	vector<int> usedIndexesX(points.size());
	vector<int> usedIndexesY(points.size());
	int nx = 0, ny = 0;

	for (unsigned i = 0; i<points.size(); i++){
		if (points[i][0] != 0){
			usedIndexesX[nx++] = i;
			if (points[i][1] != 0)
				usedIndexesY[ny++] = i;
		}
	}

	int numOfCombinations = boost::math::binomial_coefficient<double>(nx - 1, degree - 1) * ny * 0.3; // 1/3 of all combination
	int numCandidates = (numOfCombinations > n_polynomials
		? n_polynomials
		: numOfCombinations);

	TPoint minPolynomial(degree);
	double minEmpRisk = 1;
	for (int i = 0; i < numCandidates; i++){
		// generate sample
		set<int> smp;
		smp.insert(usedIndexesY[random(ny)]);
		while (smp.size() < degree){
			smp.insert(usedIndexesX[random(nx)]);
		}

		TMatrix sample(degree);
		set <int>::const_iterator s = smp.begin();
		for (unsigned j = 0; j < degree; j++, s++) {
			sample[j] = points[*s];
		}

		try{
			TPoint pol = GetPolynomial(degree, sample);
			double rsk = GetEmpiricalRisk(pol, points, numClass1);
			if (rsk < minEmpRisk) {
  			minPolynomial = pol;
  			minEmpRisk = rsk;        
			}
		}
		catch (runtime_error &e){ /* singular matrix*/ }
		catch (...){ /* NA or inf */ }
	}
	return minPolynomial;
}

#ifndef _MSC_VER

#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
double CGetEmpiricalRiskSmoothed(NumericVector& polynomial, NumericMatrix& points, int numClass1, int numClass2){
	const float smoothingConstant = 100;
	unsigned degree = polynomial.size();

	double risk = 0;
	int sign = 1;
	for (unsigned i = 0; i < numClass1 + numClass2; i++){
		if (i >= numClass1)
			sign = -1;

		double val = points(i, 0);
		double res = 0;
		for (unsigned j = 0; j < degree; j++){
			res += polynomial[j] * std::pow(val, j + 1);
		}
		risk += 1 / (1 + exp(-smoothingConstant*(points(i, 1) - res)*sign));
	}

	return risk / (numClass1 + numClass2);
}

TPoint nlm_optimize(TMatrix& points, TPoint& minCandidate, int numClass1, int numClass2){

	NumericVector r_minCandidate(minCandidate.begin(), minCandidate.end());
	NumericMatrix r_points(points.size(), 2);
	for (int i = 0; i < points.size(); i++){ r_points(i, 0) = points[i][0]; r_points(i, 1) = points[i][1]; }

  Environment env("package:ddalpha");

  Function nlm("nlm"); 
 // Function getRisk = env["CGetEmpiricalRiskSmoothed"];
  Function optimize = env["nlm_optimize_r"];

	NumericVector r_pol = optimize(r_minCandidate, r_points, numClass1, numClass2);

	return as<TPoint>(r_pol);
}

#endif

/**
Chooses the best in classification sense polynomial

@param points:       nx2 matrix of points where first column is an absciss,	n = numClass1 + numClass2
@param numClass1:    Number of points belonging to the first class
@param degree:       Degree of the polynomial
@param presize:		 if true - run evaluation 5 times

@return polynomial as a vector of coefficients starting with the first degree (a0 = 0 always)
*/
TPoint GetOptPolynomial(TMatrix& points, unsigned numClass1, unsigned degree, bool presize /*default = FALSE*/){

	double minError = 100.1;
	TPoint minPol;
	
	for (int i = 0; i < (presize ? 5 : 1); i++){
		TPoint minCandidate = GetRandomMinPolynomial(points, numClass1, degree, 10 ^ degree);

#ifdef VStudio
		// no optimization
		minPol = minCandidate;
#endif
#ifndef VStudio

#ifdef DEBUG
Rcpp::Rcout << "candminPol: ";
for (int i = 0; i< minCandidate.size(); i++){
  Rcpp::Rcout << minCandidate[i] << " ";
}
Rcpp::Rcout <<  " \n";
#endif


		TPoint optPolynomial = nlm_optimize(points, minCandidate, numClass1, points.size() - numClass1);
		int err = GetEmpiricalRisk(optPolynomial, points, numClass1);
		if (err < minError){
			minPol = optPolynomial;
			minError = err;
		}
  
#ifdef DEBUG
Rcpp::Rcout << "minPol: ";
for (int i = 0; i< minPol.size(); i++){
  Rcpp::Rcout << minPol[i] << " ";
}
Rcpp::Rcout << " ; error = "<< minError << " \n";
#endif

#endif	
		
	}  

	return(minPol);
}


/**
Calculates classification error of "degree" - degree polynomial using cross - validation approach

@param points:       nx2 matrix of points where first column is an absciss,	n = numClass1 + numClass2
@param numClass1:    Number of points belonging to the first class
@param degree:       Degree of the polynomial
@param chunkNumber:  Number of chunks in which data set should be splitted, chunkNumber should be a divider of n(n%%chunkNumber = 0)

@return Number of errors
*/
int GetCvError(TMatrix& points, unsigned numClass1, unsigned degree, unsigned chunkNumber){

	unsigned n = points.size();
	unsigned chunkSize = ceil((double)n / chunkNumber);

	TMatrix learnpoints(n - chunkSize); TMatrix checkpoints(chunkSize);

	int chunk = 0;
	int n1 = 0; // number of Class1 points in checkpoints
	for (unsigned j = 0, l = 0, c = 0; j < n; j++){
		if (j%chunkNumber)
			learnpoints[l++] = points[j];
		else
			checkpoints[c++] = points[j];
		if (j < numClass1 && (j%chunkNumber == 0)) n1++;
	}

	double err = 0;
	bool bigch = true;
	for (; chunk < chunkNumber; chunk++){

		if (chunk > 0){
			if (bigch && (chunkNumber)*(chunkSize - 1) + chunk == n){
				bigch = false;
				chunkSize--;
				checkpoints.resize(chunkSize);
				learnpoints.resize(n - chunkSize);
				learnpoints[n - chunkSize - 1] = points[n - 1];
			}

			for (int i = 0; i < chunkSize; i++){
				checkpoints[i] = learnpoints[(chunkNumber - 1)*i + chunk - 1];
				learnpoints[(chunkNumber - 1)*i + chunk - 1] = points[chunkNumber*i + chunk - 1];
				if (chunkNumber*i + chunk == numClass1)
					n1--;
			}
		}

		TPoint minPolynomial = GetOptPolynomial(learnpoints, numClass1 - n1, degree, false);
		double curErr = GetEmpiricalRisk(minPolynomial, checkpoints, n1);
		err += curErr;//  chunkSize;
	}

	return err/n;
}

TPoint PolynomialLearnCV(TMatrix& input, unsigned numClass1, unsigned int maxDegree, unsigned int chunkNumber, int *degree, int *axis){
	unsigned numPoints = input.size();

	unsigned polOptDegree = 0;
	unsigned polOptError = numPoints;
	unsigned polOptAxis = 0;


	TMatrix input2(numPoints); // copy
	for (int i = 0, tmp; i < numPoints; i++){ input2[i] = TPoint(2); input2[i][0] = input[i][1]; input2[i][1] = input[i][0]; } // swap columns

	for (int degree = 1; degree <= maxDegree; degree++){
		double polError = GetCvError(input, numClass1, degree, chunkNumber);
		if (polError < polOptError){
			polOptAxis = 0;
			polOptDegree = degree;
			polOptError = polError;
		}
		polError = GetCvError(input, numClass1, degree, chunkNumber);
		if (polError < polOptError){
			polOptAxis = 1;
			polOptDegree = degree;
			polOptError = polError;
		}
	}
	TPoint polynomial = polOptAxis == 0
		? GetOptPolynomial(input, numClass1, polOptDegree, true) 
		: GetOptPolynomial(input2, numClass1, polOptDegree, true);

	*axis = polOptAxis;
	*degree = polOptDegree;

	return polynomial;
}
