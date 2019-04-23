// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// b_updateC
arma::colvec b_updateC(arma::mat X, arma::colvec y, arma::colvec u, arma::colvec z, double rho, unsigned int maxiter, double toler, double b, double alpha);
RcppExport SEXP _aeffp_b_updateC(SEXP XSEXP, SEXP ySEXP, SEXP uSEXP, SEXP zSEXP, SEXP rhoSEXP, SEXP maxiterSEXP, SEXP tolerSEXP, SEXP bSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type toler(tolerSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(b_updateC(X, y, u, z, rho, maxiter, toler, b, alpha));
    return rcpp_result_gen;
END_RCPP
}
// admmlasso_logC
arma::colvec admmlasso_logC(arma::mat X, arma::colvec y, double lam, double rho, unsigned int maxit, double tol);
RcppExport SEXP _aeffp_admmlasso_logC(SEXP XSEXP, SEXP ySEXP, SEXP lamSEXP, SEXP rhoSEXP, SEXP maxitSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(admmlasso_logC(X, y, lam, rho, maxit, tol));
    return rcpp_result_gen;
END_RCPP
}
// b_updateC_tabled
arma::colvec b_updateC_tabled(arma::mat X, arma::mat y, arma::colvec u, arma::colvec z, double rho, unsigned int maxiter, double toler, double b, double alpha);
RcppExport SEXP _aeffp_b_updateC_tabled(SEXP XSEXP, SEXP ySEXP, SEXP uSEXP, SEXP zSEXP, SEXP rhoSEXP, SEXP maxiterSEXP, SEXP tolerSEXP, SEXP bSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type toler(tolerSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(b_updateC_tabled(X, y, u, z, rho, maxiter, toler, b, alpha));
    return rcpp_result_gen;
END_RCPP
}
// admmlasso_logC_tabled
arma::colvec admmlasso_logC_tabled(arma::mat X, arma::mat y, double lam, double rho, unsigned int maxit, double tol);
RcppExport SEXP _aeffp_admmlasso_logC_tabled(SEXP XSEXP, SEXP ySEXP, SEXP lamSEXP, SEXP rhoSEXP, SEXP maxitSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type lam(lamSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(admmlasso_logC_tabled(X, y, lam, rho, maxit, tol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_aeffp_b_updateC", (DL_FUNC) &_aeffp_b_updateC, 9},
    {"_aeffp_admmlasso_logC", (DL_FUNC) &_aeffp_admmlasso_logC, 6},
    {"_aeffp_b_updateC_tabled", (DL_FUNC) &_aeffp_b_updateC_tabled, 9},
    {"_aeffp_admmlasso_logC_tabled", (DL_FUNC) &_aeffp_admmlasso_logC_tabled, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_aeffp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}