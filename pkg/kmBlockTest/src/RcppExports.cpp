// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "kmBlock_types.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// meanByBlocks
Rcpp::List meanByBlocks(const Array& M, const IVector& clu, const int dimensions, const IVector& n, const std::string diagonal, const std::string& sBorders, const Rcpp::Nullable<Array>& bordersMatLower, const Rcpp::Nullable<Array>& bordersMatUpper, const Rcpp::Nullable<DMatrix>& bordersSeperateLower, const Rcpp::Nullable<DMatrix>& bordersSeperateUpper);
RcppExport SEXP _kmBlockTest_meanByBlocks(SEXP MSEXP, SEXP cluSEXP, SEXP dimensionsSEXP, SEXP nSEXP, SEXP diagonalSEXP, SEXP sBordersSEXP, SEXP bordersMatLowerSEXP, SEXP bordersMatUpperSEXP, SEXP bordersSeperateLowerSEXP, SEXP bordersSeperateUpperSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Array& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const IVector& >::type clu(cluSEXP);
    Rcpp::traits::input_parameter< const int >::type dimensions(dimensionsSEXP);
    Rcpp::traits::input_parameter< const IVector& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const std::string >::type diagonal(diagonalSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type sBorders(sBordersSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Array>& >::type bordersMatLower(bordersMatLowerSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Array>& >::type bordersMatUpper(bordersMatUpperSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<DMatrix>& >::type bordersSeperateLower(bordersSeperateLowerSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<DMatrix>& >::type bordersSeperateUpper(bordersSeperateUpperSEXP);
    rcpp_result_gen = Rcpp::wrap(meanByBlocks(M, clu, dimensions, n, diagonal, sBorders, bordersMatLower, bordersMatUpper, bordersSeperateLower, bordersSeperateUpper));
    return rcpp_result_gen;
END_RCPP
}
// kmBlock
Rcpp::List kmBlock(const Array& M, const IVector& clu, const Array& weights, const IVector& n, const IVector& nClu, const std::string& diagonal, const std::string& sBorders, const Rcpp::Nullable<Array>& bordersMatLower, const Rcpp::Nullable<Array>& bordersMatUpper, const Rcpp::Nullable<DMatrix>& bordersSeperateLower, const Rcpp::Nullable<DMatrix>& bordersSeperateUpper);
RcppExport SEXP _kmBlockTest_kmBlock(SEXP MSEXP, SEXP cluSEXP, SEXP weightsSEXP, SEXP nSEXP, SEXP nCluSEXP, SEXP diagonalSEXP, SEXP sBordersSEXP, SEXP bordersMatLowerSEXP, SEXP bordersMatUpperSEXP, SEXP bordersSeperateLowerSEXP, SEXP bordersSeperateUpperSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Array& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const IVector& >::type clu(cluSEXP);
    Rcpp::traits::input_parameter< const Array& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const IVector& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const IVector& >::type nClu(nCluSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type diagonal(diagonalSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type sBorders(sBordersSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Array>& >::type bordersMatLower(bordersMatLowerSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Array>& >::type bordersMatUpper(bordersMatUpperSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<DMatrix>& >::type bordersSeperateLower(bordersSeperateLowerSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<DMatrix>& >::type bordersSeperateUpper(bordersSeperateUpperSEXP);
    rcpp_result_gen = Rcpp::wrap(kmBlock(M, clu, weights, n, nClu, diagonal, sBorders, bordersMatLower, bordersMatUpper, bordersSeperateLower, bordersSeperateUpper));
    return rcpp_result_gen;
END_RCPP
}
// critFunction
double critFunction(const Array& M, const IVector& clu, const Array& weights, const int dimensions, const IVector& n, const std::string& diagonal, const std::string& sBorders, const Rcpp::Nullable<Array>& bordersMatLower, const Rcpp::Nullable<Array>& bordersMatUpper, const Rcpp::Nullable<DMatrix>& bordersSeperateLower, const Rcpp::Nullable<DMatrix>& bordersSeperateUpper);
RcppExport SEXP _kmBlockTest_critFunction(SEXP MSEXP, SEXP cluSEXP, SEXP weightsSEXP, SEXP dimensionsSEXP, SEXP nSEXP, SEXP diagonalSEXP, SEXP sBordersSEXP, SEXP bordersMatLowerSEXP, SEXP bordersMatUpperSEXP, SEXP bordersSeperateLowerSEXP, SEXP bordersSeperateUpperSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Array& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const IVector& >::type clu(cluSEXP);
    Rcpp::traits::input_parameter< const Array& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const int >::type dimensions(dimensionsSEXP);
    Rcpp::traits::input_parameter< const IVector& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type diagonal(diagonalSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type sBorders(sBordersSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Array>& >::type bordersMatLower(bordersMatLowerSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<Array>& >::type bordersMatUpper(bordersMatUpperSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<DMatrix>& >::type bordersSeperateLower(bordersSeperateLowerSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<DMatrix>& >::type bordersSeperateUpper(bordersSeperateUpperSEXP);
    rcpp_result_gen = Rcpp::wrap(critFunction(M, clu, weights, dimensions, n, diagonal, sBorders, bordersMatLower, bordersMatUpper, bordersSeperateLower, bordersSeperateUpper));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_kmBlockTest_meanByBlocks", (DL_FUNC) &_kmBlockTest_meanByBlocks, 10},
    {"_kmBlockTest_kmBlock", (DL_FUNC) &_kmBlockTest_kmBlock, 11},
    {"_kmBlockTest_critFunction", (DL_FUNC) &_kmBlockTest_critFunction, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_kmBlockTest(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
