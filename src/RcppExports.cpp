// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// scale
arma::mat scale(const arma::mat& X, int dim);
RcppExport SEXP _X_ING_scale(SEXP XSEXP, SEXP dimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type dim(dimSEXP);
    rcpp_result_gen = Rcpp::wrap(scale(X, dim));
    return rcpp_result_gen;
END_RCPP
}
// unique
std::vector<int> unique(std::vector<int>& vec);
RcppExport SEXP _X_ING_unique(SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int>& >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(unique(vec));
    return rcpp_result_gen;
END_RCPP
}
// Lq_func
double Lq_func(arma::mat sjk_sqr, arma::mat mujk, arma::mat alphajk, arma::vec vk_sqr, arma::mat pik, arma::mat Lambda, arma::mat betahat);
RcppExport SEXP _X_ING_Lq_func(SEXP sjk_sqrSEXP, SEXP mujkSEXP, SEXP alphajkSEXP, SEXP vk_sqrSEXP, SEXP pikSEXP, SEXP LambdaSEXP, SEXP betahatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sjk_sqr(sjk_sqrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mujk(mujkSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type alphajk(alphajkSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type vk_sqr(vk_sqrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type pik(pikSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type betahat(betahatSEXP);
    rcpp_result_gen = Rcpp::wrap(Lq_func(sjk_sqr, mujk, alphajk, vk_sqr, pik, Lambda, betahat));
    return rcpp_result_gen;
END_RCPP
}
// XING_starting
List XING_starting(arma::mat betahat, arma::mat Lambda, int iterT, double eps_thres, double vk_init, double bound, double pi_init);
RcppExport SEXP _X_ING_XING_starting(SEXP betahatSEXP, SEXP LambdaSEXP, SEXP iterTSEXP, SEXP eps_thresSEXP, SEXP vk_initSEXP, SEXP boundSEXP, SEXP pi_initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< int >::type iterT(iterTSEXP);
    Rcpp::traits::input_parameter< double >::type eps_thres(eps_thresSEXP);
    Rcpp::traits::input_parameter< double >::type vk_init(vk_initSEXP);
    Rcpp::traits::input_parameter< double >::type bound(boundSEXP);
    Rcpp::traits::input_parameter< double >::type pi_init(pi_initSEXP);
    rcpp_result_gen = Rcpp::wrap(XING_starting(betahat, Lambda, iterT, eps_thres, vk_init, bound, pi_init));
    return rcpp_result_gen;
END_RCPP
}
// XING_single_data
List XING_single_data(const arma::mat& betahat, const arma::mat& Lambda, int PC, List results_alg1, double eps_thresh, int iterT, bool use_true_X, double bound);
RcppExport SEXP _X_ING_XING_single_data(SEXP betahatSEXP, SEXP LambdaSEXP, SEXP PCSEXP, SEXP results_alg1SEXP, SEXP eps_threshSEXP, SEXP iterTSEXP, SEXP use_true_XSEXP, SEXP boundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< int >::type PC(PCSEXP);
    Rcpp::traits::input_parameter< List >::type results_alg1(results_alg1SEXP);
    Rcpp::traits::input_parameter< double >::type eps_thresh(eps_threshSEXP);
    Rcpp::traits::input_parameter< int >::type iterT(iterTSEXP);
    Rcpp::traits::input_parameter< bool >::type use_true_X(use_true_XSEXP);
    Rcpp::traits::input_parameter< double >::type bound(boundSEXP);
    rcpp_result_gen = Rcpp::wrap(XING_single_data(betahat, Lambda, PC, results_alg1, eps_thresh, iterT, use_true_X, bound));
    return rcpp_result_gen;
END_RCPP
}
// XING
List XING(const arma::mat& betahat1, const arma::mat& betahat2, const arma::mat& Lambda1, const arma::mat& Lambda2, int CC, int PC1, int PC2, const List& results_alg1_dat1, const List& results_alg1_dat2, const List& results_alg2_dat1, const List& results_alg2_dat2, double eps_thresh, int iterT, double bound);
RcppExport SEXP _X_ING_XING(SEXP betahat1SEXP, SEXP betahat2SEXP, SEXP Lambda1SEXP, SEXP Lambda2SEXP, SEXP CCSEXP, SEXP PC1SEXP, SEXP PC2SEXP, SEXP results_alg1_dat1SEXP, SEXP results_alg1_dat2SEXP, SEXP results_alg2_dat1SEXP, SEXP results_alg2_dat2SEXP, SEXP eps_threshSEXP, SEXP iterTSEXP, SEXP boundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type betahat1(betahat1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type betahat2(betahat2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Lambda1(Lambda1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Lambda2(Lambda2SEXP);
    Rcpp::traits::input_parameter< int >::type CC(CCSEXP);
    Rcpp::traits::input_parameter< int >::type PC1(PC1SEXP);
    Rcpp::traits::input_parameter< int >::type PC2(PC2SEXP);
    Rcpp::traits::input_parameter< const List& >::type results_alg1_dat1(results_alg1_dat1SEXP);
    Rcpp::traits::input_parameter< const List& >::type results_alg1_dat2(results_alg1_dat2SEXP);
    Rcpp::traits::input_parameter< const List& >::type results_alg2_dat1(results_alg2_dat1SEXP);
    Rcpp::traits::input_parameter< const List& >::type results_alg2_dat2(results_alg2_dat2SEXP);
    Rcpp::traits::input_parameter< double >::type eps_thresh(eps_threshSEXP);
    Rcpp::traits::input_parameter< int >::type iterT(iterTSEXP);
    Rcpp::traits::input_parameter< double >::type bound(boundSEXP);
    rcpp_result_gen = Rcpp::wrap(XING(betahat1, betahat2, Lambda1, Lambda2, CC, PC1, PC2, results_alg1_dat1, results_alg1_dat2, results_alg2_dat1, results_alg2_dat2, eps_thresh, iterT, bound));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_X_ING_scale", (DL_FUNC) &_X_ING_scale, 2},
    {"_X_ING_unique", (DL_FUNC) &_X_ING_unique, 1},
    {"_X_ING_Lq_func", (DL_FUNC) &_X_ING_Lq_func, 7},
    {"_X_ING_XING_starting", (DL_FUNC) &_X_ING_XING_starting, 7},
    {"_X_ING_XING_single_data", (DL_FUNC) &_X_ING_XING_single_data, 8},
    {"_X_ING_XING", (DL_FUNC) &_X_ING_XING, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_X_ING(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}