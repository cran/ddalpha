// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// CGetEmpiricalRiskSmoothed
double CGetEmpiricalRiskSmoothed(NumericVector& polynomial, NumericMatrix& points, int numClass1, int numClass2);
RcppExport SEXP ddalpha_CGetEmpiricalRiskSmoothed(SEXP polynomialSEXP, SEXP pointsSEXP, SEXP numClass1SEXP, SEXP numClass2SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector& >::type polynomial(polynomialSEXP );
        Rcpp::traits::input_parameter< NumericMatrix& >::type points(pointsSEXP );
        Rcpp::traits::input_parameter< int >::type numClass1(numClass1SEXP );
        Rcpp::traits::input_parameter< int >::type numClass2(numClass2SEXP );
        double __result = CGetEmpiricalRiskSmoothed(polynomial, points, numClass1, numClass2);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}