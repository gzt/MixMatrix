// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// posdefsqrt
arma::mat posdefsqrt(arma::mat& x);
RcppExport SEXP _matrixdist_posdefsqrt(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(posdefsqrt(x));
    return rcpp_result_gen;
END_RCPP
}
// posdefinvsqrt
arma::mat posdefinvsqrt(arma::mat& x);
RcppExport SEXP _matrixdist_posdefinvsqrt(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(posdefinvsqrt(x));
    return rcpp_result_gen;
END_RCPP
}
// testsymmetric
bool testsymmetric(arma::mat x, double tol);
RcppExport SEXP _matrixdist_testsymmetric(SEXP xSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(testsymmetric(x, tol));
    return rcpp_result_gen;
END_RCPP
}
// dmatnorm_calc
arma::colvec dmatnorm_calc(arma::cube& x, arma::mat& mean, arma::mat& U, arma::mat& V);
RcppExport SEXP _matrixdist_dmatnorm_calc(SEXP xSEXP, SEXP meanSEXP, SEXP USEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type U(USEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(dmatnorm_calc(x, mean, U, V));
    return rcpp_result_gen;
END_RCPP
}
// dmat_t_calc
arma::colvec dmat_t_calc(arma::cube& x, double df, arma::mat& mean, arma::mat& U, arma::mat& V);
RcppExport SEXP _matrixdist_dmat_t_calc(SEXP xSEXP, SEXP dfSEXP, SEXP meanSEXP, SEXP USEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type U(USEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(dmat_t_calc(x, df, mean, U, V));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_matrixdist_posdefsqrt", (DL_FUNC) &_matrixdist_posdefsqrt, 1},
    {"_matrixdist_posdefinvsqrt", (DL_FUNC) &_matrixdist_posdefinvsqrt, 1},
    {"_matrixdist_testsymmetric", (DL_FUNC) &_matrixdist_testsymmetric, 2},
    {"_matrixdist_dmatnorm_calc", (DL_FUNC) &_matrixdist_dmatnorm_calc, 4},
    {"_matrixdist_dmat_t_calc", (DL_FUNC) &_matrixdist_dmat_t_calc, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_matrixdist(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
