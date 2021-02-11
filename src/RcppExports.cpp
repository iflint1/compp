// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// draw_points_cpp
Rcpp::NumericMatrix draw_points_cpp(int total_points, Rcpp::Function draw_points_lambda, double lambda_max);
RcppExport SEXP _compp_draw_points_cpp(SEXP total_pointsSEXP, SEXP draw_points_lambdaSEXP, SEXP lambda_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type total_points(total_pointsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type draw_points_lambda(draw_points_lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_max(lambda_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(draw_points_cpp(total_points, draw_points_lambda, lambda_max));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_compp_draw_points_cpp", (DL_FUNC) &_compp_draw_points_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_compp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}