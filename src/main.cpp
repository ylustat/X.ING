#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]

// Custom implementation of scale function
// [[Rcpp::export]]
arma::mat scale(const arma::mat& X, int dim = 0) {
  arma::mat scaledX = X;

  if (dim == 0) {
    // Scale along columns
    for (int j = 0; j < X.n_cols; ++j) {
      double meanVal = arma::mean(X.col(j));
      double sdVal = arma::stddev(X.col(j), 1);

      if (sdVal != 0.0) {
        scaledX.col(j) = (X.col(j) - meanVal) / sdVal;
      } else {
        scaledX.col(j).zeros();
      }
    }
  } else if (dim == 1) {
    // Scale along rows
    for (int i = 0; i < X.n_rows; ++i) {
      double meanVal = arma::mean(X.row(i));
      double sdVal = arma::stddev(X.row(i), 1);

      if (sdVal != 0.0) {
        scaledX.row(i) = (X.row(i) - meanVal) / sdVal;
      } else {
        scaledX.row(i).zeros();
      }
    }
  }

  return scaledX;
}
// [[Rcpp::export]]
std::vector<int> unique(std::vector<int>& vec) {
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
  return vec;
}

// perform CCA
Environment pkg = Environment::namespace_env("CCA");
Function cc = pkg["cc"];

Environment pkg2 = Environment::namespace_env("RGCCA");
Function rgcca = pkg2["rgcca"];

// Lower bound Lq
// [[Rcpp::export]]
double Lq_func(arma::mat sjk_sqr, arma::mat mujk, arma::mat alphajk, arma::vec vk_sqr, arma::mat pik, arma::mat Lambda, arma::mat betahat){
  int M = betahat.n_rows;
  int K = betahat.n_cols;
  arma::vec inner1(M);
  arma::vec inner2(M);
  arma::vec inner3(M);
  arma::vec inner4(M);
  arma::mat B = betahat - alphajk % mujk;
  arma::mat one_mat = arma::ones<arma::mat>(M,K);
  inner1 = sum((B * Lambda) % B, 1);
  inner2 = sum((alphajk % (arma::square(mujk) + sjk_sqr) - arma::square(alphajk) % arma::square(mujk)) * Lambda.diag(), 1);
  inner3 = sum(alphajk % (1 + (log(sjk_sqr)-one_mat * arma::diagmat(log(vk_sqr))) - (sjk_sqr + arma::square(mujk)) * arma::diagmat(1/vk_sqr)), 1);
  inner4 = sum(alphajk % (log(alphajk) - log(pik)) + (1 - alphajk) % (log(1 - alphajk) - log(1 - pik)), 1);

  double Lq = -0.5 * accu(inner1) - 0.5 * accu(inner2) + 0.5 * accu(inner3) - accu(inner4);
  return(Lq);
}

// [[Rcpp::export]]
List XING_starting(arma::mat betahat, arma::mat Lambda, int iterT = 10, double eps_thres = 1e-3,
               double vk_init = 0.1, double bound = 1e-4, double pi_init = 0.1) {

  int M = betahat.n_rows;
  int K = betahat.n_cols;

  arma::mat sjk_sqrm1(M, K, arma::fill::ones);
  arma::mat mujkm1(M, K, arma::fill::zeros);

  arma::vec vk_sqrm1(K);
  vk_sqrm1.fill(vk_init);

  arma::mat alphajkm1(M, K);
  alphajkm1.fill(pi_init);

  arma::mat pikm1(M, K);
  pikm1.fill(pi_init);

  arma::mat ujk(M, K, arma::fill::zeros);

  arma::vec Lq_iter(iterT, arma::fill::zeros);
  Lq_iter[0] = -INFINITY;

  double upbound = 1 - bound;
  double lowbound = 0 + bound;

  for(int t = 1; t < iterT; ++t) {

    // E step
    for (int k = 0; k < K; ++k) {
      sjk_sqrm1.col(k).fill(1 / (Lambda.diag()[k] + (1 / vk_sqrm1[k])));
    }
    mujkm1 = (betahat * Lambda - ((alphajkm1 % mujkm1) * Lambda - (alphajkm1 % mujkm1) * arma::diagmat(Lambda.diag()))) * arma::diagmat(1 / (Lambda.diag() + arma::pow(vk_sqrm1, -1)));
    ujk = log(pikm1 / (1 - pikm1)) + 0.5 * ((log(sjk_sqrm1) - arma::ones<arma::mat>(M, K) * arma::diagmat(log(vk_sqrm1))) + arma::pow(mujkm1, 2) / sjk_sqrm1);

    alphajkm1 = 1 / (1 + arma::exp(-ujk));
    alphajkm1 = arma::clamp(alphajkm1, lowbound, upbound);
    // M step
    for(int k = 0; k < K; ++k) {
      vk_sqrm1[k] = sum(alphajkm1.col(k) % (sjk_sqrm1.col(k) + arma::pow(mujkm1.col(k), 2))) / sum(alphajkm1.col(k));
    }
    for(int k = 0; k < K; ++k) {
      pikm1.col(k).fill(mean(alphajkm1.col(k)));
      pikm1.col(k).fill(clamp(pikm1(0, k), lowbound, upbound));
    }
    Lq_iter[t] = Lq_func(sjk_sqrm1, mujkm1, alphajkm1, vk_sqrm1, pikm1, Lambda, betahat); // Lq_func is assumed to be already implemented

    // if(Lq_iter[t] < Lq_iter[t - 1]){
    //   Rcout << "Lq value decreases!" << "\n";
    // } else if(Lq_iter[t] - Lq_iter[t - 1] < eps_thres){
    //   Rcout << "Algorithm converges in iterations = " << t - 1 << "\n";
    //   break;
    // }
  }

  return List::create(Named("sjk_sqrm1") = sjk_sqrm1,
                      Named("mujkm1") = mujkm1,
                      Named("betahat") = betahat,
                      Named("Lambda") = Lambda,
                      Named("vk_sqrm1") = vk_sqrm1,
                      Named("pikm1") = pikm1,
                      Named("alphajkm1") = alphajkm1,
                      Named("Lq_iter") = Lq_iter);
}


// [[Rcpp::export]]
List XING_single_data(const arma::mat& betahat, const arma::mat& Lambda, int PC = 2,
               List results_alg1 = R_NilValue, double eps_thresh = 1e-2,
               int iterT = 5, bool use_true_X = false, double bound = 1e-4) {

  int M = betahat.n_rows;
  int K = betahat.n_cols;
  double upbound = 1 - bound;
  double lowbound = bound;

  // Define and initialize variables
  arma::mat sjk_sqrm2 = results_alg1["sjk_sqrm1"];
  arma::mat mujkm2 = results_alg1["mujkm1"];
  arma::mat alphajkm2 = results_alg1["alphajkm1"];
  arma::vec vk_sqrm2 = results_alg1["vk_sqrm1"];
  arma::mat pikm2 = results_alg1["alphajkm1"];

  arma::mat pi = results_alg1["pikm1"];
  arma::mat x0k = log(pi / (1 - pi));

  arma::mat X_full(M, K);
  X_full.fill(0.0);
  arma::mat X_svd(M, K);
  X_svd.fill(0.0);
  arma::mat X_red(M, K);
  X_red.fill(0.0);
  arma::mat ujk2(M, K);
  ujk2.fill(0.0);

  std::vector<arma::mat> X_red_list(iterT);
  arma::vec Lq_iter(iterT);
  Lq_iter.fill(0.0);
  Lq_iter[0] = -arma::datum::inf;

  for (int t = 1; t < iterT; t++) {
    // E step
    for (int k = 0; k < K; k++) {
      sjk_sqrm2.col(k).fill(1 / (Lambda.diag()[k] + 1 / vk_sqrm2[k]));
    }

    mujkm2 = (betahat * Lambda -
      ((alphajkm2 % mujkm2) * Lambda -
      (alphajkm2 % mujkm2) * diagmat(Lambda.diag()))) *
      diagmat(1.0 / (Lambda.diag() + 1/vk_sqrm2));

    ujk2 = log(pikm2 / (1 - pikm2)) +
      0.5 * ((log(sjk_sqrm2) - repmat(log(vk_sqrm2).t(), M, 1)) +
      pow(mujkm2, 2) / sjk_sqrm2);

    alphajkm2 = 1.0 / (1.0 + exp(-ujk2));
    alphajkm2.elem(find(alphajkm2 > upbound)).fill(upbound);
    alphajkm2.elem(find(alphajkm2 < lowbound)).fill(lowbound);

    // M step
    X_full = log(alphajkm2)-log(1-alphajkm2) - x0k;

    arma::mat U, V;
    arma::vec S;
    svd_econ(U, S, V, scale(X_full));  // equivalent to svd function in R

    X_red = U.cols(0, PC-1) * diagmat(S.head(PC)) * trans(V.cols(0, PC-1));  // low-rank approximation
    arma::mat sd1 = stddev(X_full, 0, 0);
    arma::mat mean1 = mean(X_full, 0);

    X_red = X_red * diagmat(sd1) + repmat(mean1,  X_red.n_rows, 1);


    for(int k = 0; k < K; k++){
      vk_sqrm2[k] = sum(alphajkm2.col(k) % (sjk_sqrm2.col(k) + pow(mujkm2.col(k), 2))) / sum(alphajkm2.col(k));
    }
    pikm2 = 1.0 / (1.0 + exp(-X_red - x0k));
    pikm2.elem(find(pikm2 > upbound)).fill(upbound);
    pikm2.elem(find(pikm2 < lowbound)).fill(lowbound);

    //Lq compute
    Lq_iter[t] = Lq_func(sjk_sqrm2, mujkm2, alphajkm2, vk_sqrm2, pikm2, Lambda, betahat);  // check this line. Function implementation should be provided

    // if (Lq_iter[t] < Lq_iter[t - 1]) {
    //   Rcpp::Rcout << "Lq value decreases!" << std::endl;
    // } else if (Lq_iter[t] - Lq_iter[t - 1] < eps_thresh) {
    //   Rcpp::Rcout << "Algorithm converges in iterations = " << t - 1 << std::endl;
    //   break;
    // }
  }

  return List::create(_["sjk_sqrm1"] = sjk_sqrm2, _["mujkm1"] = mujkm2, _["betahat"] = betahat,
                      _["Lambda"] = Lambda, _["vk_sqrm1"] = vk_sqrm2, _["pikm1"] = pikm2, _["alphajkm1"] = alphajkm2,
                        _["X_est"] = X_red, _["Lq_iter"] = Lq_iter, _["x0k"] = x0k, _["X_full"] = X_full, _["X_red_list"] = X_red_list);
}

// [[Rcpp::export]]
List XING(const arma::mat& betahat1, const arma::mat& betahat2, const arma::mat& Lambda1, const arma::mat& Lambda2,
          int CC, int PC1, int PC2, const List& results_alg1_dat1, const List& results_alg1_dat2,
          const List& results_alg2_dat1, const List& results_alg2_dat2, double eps_thresh = 1e-2,
          int iterT=20,  double bound = 1e-4) {

  double upbound = 1 - bound;
  double lowbound = bound;

  int M = betahat1.n_rows;
  int K1 = betahat1.n_cols;
  int K2 = betahat2.n_cols;
  // Defining and initializing new variables using results from Alg 1 and 2
  arma::mat sjk_sqrm2_dat1 = results_alg2_dat1["sjk_sqrm1"]; // posterior variance for beta
  arma::mat mujkm2_dat1 = results_alg2_dat1["mujkm1"]; // posterior mean for beta
  arma::mat sjk_sqrm2_dat2 = results_alg2_dat2["sjk_sqrm1"]; // posterior variance for beta
  arma::mat mujkm2_dat2 = results_alg2_dat2["mujkm1"]; // posterior mean for beta
  // posterior for association probability pi_jk
  arma::mat alphajkm2_dat1 = results_alg2_dat1["alphajkm1"];
  arma::mat alphajkm2_dat2 = results_alg2_dat2["alphajkm1"];
  // Parameters
  arma::vec vk_sqrm2_dat1 = results_alg2_dat1["vk_sqrm1"];
  arma::vec vk_sqrm2_dat2 = results_alg2_dat2["vk_sqrm1"];
  arma::mat pikm2_dat1 = results_alg2_dat1["alphajkm1"];
  arma::mat pikm2_dat2 = results_alg2_dat2["alphajkm1"];
  arma::mat p1 = results_alg1_dat1["pikm1"];
  arma::mat p2 = results_alg1_dat2["pikm1"];
  arma::mat x0k_dat1 = log(p1 / (1 - p1));
  arma::mat x0k_dat2 = log(p2 / (1 - p2));
  arma::mat X_full_dat1 = arma::zeros<arma::mat>(M, K1);
  arma::mat X_svd_res1 = arma::zeros<arma::mat>(M, K1);
  arma::mat X_red_res1 = arma::zeros<arma::mat>(M, K1);
  arma::mat X_est1 = arma::zeros<arma::mat>(M, K1);
  arma::mat X_est2 = arma::zeros<arma::mat>(M, K2);
  arma::mat X_res1 = arma::zeros<arma::mat>(M, K1);
  arma::mat X_res2 = arma::zeros<arma::mat>(M, K2);
  arma::mat X_full_dat2 = arma::zeros<arma::mat>(M, K2);
  arma::mat X_red_res2 = arma::zeros<arma::mat>(M, K2);
  arma::mat X_svd_res2 = arma::zeros<arma::mat>(M, K2);
  arma::mat Lambdam2_dat1 = Lambda1;
  arma::mat betahatm2_dat1 = betahat1;
  arma::mat Lambdam2_dat2 = Lambda2;
  arma::mat betahatm2_dat2 = betahat2;

  // EM algorithm with variational approximation
  // Total number of iterations (upper bound)
  // Stopping rule threshold
  arma::mat ujk2_dat1 = arma::zeros<arma::mat>(M, K1);
  arma::mat ujk2_dat2 = arma::zeros<arma::mat>(M, K2);

  // Vector for storing lower bound to be maximized
  arma::vec Lq_iter_dat1 = arma::zeros<arma::vec>(iterT);
  arma::vec Lq_iter_dat2 = arma::zeros<arma::vec>(iterT);

  Lq_iter_dat1(0) = -std::numeric_limits<double>::infinity();
  Lq_iter_dat2(0) = -std::numeric_limits<double>::infinity();

  for (int t = 1; t < iterT; ++t) {
    // E step
    for (int k = 0; k < K1; ++k) {
      sjk_sqrm2_dat1.col(k).fill(1 / (Lambda1.diag()[k] + 1 / vk_sqrm2_dat1[k]));
    }

    mujkm2_dat1 = (betahatm2_dat1 * Lambda1 -
      ((alphajkm2_dat1 % mujkm2_dat1) * Lambda1 -
      (alphajkm2_dat1 % mujkm2_dat1) * diagmat(Lambda1.diag())))
      * diagmat(1/(Lambda1.diag() + (1 / vk_sqrm2_dat1)));


    ujk2_dat1 = log(pikm2_dat1 / (1 - pikm2_dat1)) +
    0.5 * ((log(sjk_sqrm2_dat1) - repmat(log(vk_sqrm2_dat1).t(),M,1)) +
    pow(mujkm2_dat1, 2) / sjk_sqrm2_dat1);
    alphajkm2_dat1 = 1 / (1 + exp(-ujk2_dat1));
    alphajkm2_dat1 = clamp(alphajkm2_dat1, lowbound, upbound);

    for (int k = 0; k < K2; ++k) {
      sjk_sqrm2_dat2.col(k).fill(1 / (Lambda2.diag()[k] + 1 / vk_sqrm2_dat2[k]));
    }
    mujkm2_dat2 = (betahatm2_dat2 * Lambda2 -
      ((alphajkm2_dat2 % mujkm2_dat2) * Lambda2 -
      (alphajkm2_dat2 % mujkm2_dat2) * diagmat(Lambda2.diag())))
      * diagmat(1/(Lambda2.diag() + (1 / vk_sqrm2_dat2)));

    ujk2_dat2 = log(pikm2_dat2 / (1 - pikm2_dat2)) +
    0.5 * ((log(sjk_sqrm2_dat2) - repmat(log(vk_sqrm2_dat2).t(),M,1)) +
    pow(mujkm2_dat2, 2) / sjk_sqrm2_dat2);
    alphajkm2_dat2 = 1 / (1 + exp(-ujk2_dat2));
    alphajkm2_dat2 = clamp(alphajkm2_dat2, lowbound, upbound);

    // M step
    // Update X
    X_full_dat1 = log(alphajkm2_dat1) - log(1 - alphajkm2_dat1) - x0k_dat1;
    X_full_dat2 = log(alphajkm2_dat2) - log(1 - alphajkm2_dat2) - x0k_dat2;

    arma::mat sd1 = stddev(X_full_dat1, 0, 0);
    arma::mat mean1 = mean(X_full_dat1, 0);
    X_full_dat1 = scale(X_full_dat1);

    arma::mat sd2 = stddev(X_full_dat2, 0, 0);
    arma::mat mean2 = mean(X_full_dat2, 0);

    X_full_dat2 = scale(X_full_dat2);
    for (int k = 0; k < K1; ++k) {
      vk_sqrm2_dat1(k) = sum(alphajkm2_dat1.col(k) % (sjk_sqrm2_dat1.col(k) + pow(mujkm2_dat1.col(k), 2))) / sum(alphajkm2_dat1.col(k));
    }

    for (int k = 0; k < K2; ++k) {
      vk_sqrm2_dat2(k) = sum(alphajkm2_dat2.col(k) % (sjk_sqrm2_dat2.col(k) + pow(mujkm2_dat2.col(k), 2))) / sum(alphajkm2_dat2.col(k));
    }
    List ccas = Rcpp::as<List>(cc(X_full_dat1, X_full_dat2));
    arma::mat Ahat = ccas["xcoef"];
    arma::mat Bhat = ccas["ycoef"];

    arma::mat XX = X_full_dat1 * Ahat;
    arma::mat YY = X_full_dat2 * Bhat;
    arma::mat X_est1 = XX.cols(0, CC - 1) * pinv(Ahat.cols(0, CC - 1));
    arma::mat X_est2 = YY.cols(0, CC - 1) * pinv(Bhat.cols(0, CC - 1));

    arma::mat X_res1 = X_full_dat1 - X_est1;
    arma::mat X_res2 = X_full_dat2 - X_est2;

    X_est1 = X_est1 * diagmat(sd1) + repmat(mean1,  X_est1.n_rows, 1);
    X_est2 = X_est2 * diagmat(sd2) + repmat(mean2,  X_est2.n_rows, 1);

    X_res1 = X_full_dat1 - X_est1;
    X_res2 = X_full_dat2 - X_est2;
    arma::mat U1, V1, U2, V2;
    arma::vec S1, S2;
    svd_econ(U1, S1, V1, scale(X_res1));  // equivalent to svd function in R
    X_red_res1 = U1.cols(0, PC1-1) * diagmat(S1.head(PC1)) * trans(V1.cols(0, PC1-1));  // low-rank approximation
    sd1 = stddev(X_res1, 0, 0);
    mean1 = mean(X_res1, 0);
    X_red_res1 = X_red_res1 * diagmat(sd1) + repmat(mean1,  X_red_res1.n_rows, 1);

    svd_econ(U2, S2, V2, scale(X_res2));  // equivalent to svd function in R
    X_red_res2 = U2.cols(0, PC2-1) * diagmat(S2.head(PC2)) * trans(V2.cols(0, PC2-1));  // low-rank approximation
    sd2 = stddev(X_res2, 0, 0);
    mean2 = mean(X_res2, 0);
    X_red_res2 = X_red_res2 * diagmat(sd2) + repmat(mean2,  X_red_res2.n_rows, 1);


    pikm2_dat1 = 1 / (1 + exp(-X_est1 - X_red_res1 - x0k_dat1));
    pikm2_dat1 = clamp(pikm2_dat1, lowbound, upbound);

    pikm2_dat2 = 1 / (1 + exp(-X_est2 - X_red_res2 - x0k_dat2));
    pikm2_dat2 = clamp(pikm2_dat2, lowbound, upbound);

    // Lq compute
    Lq_iter_dat1(t) = Lq_func(sjk_sqrm2_dat1, mujkm2_dat1, alphajkm2_dat1, vk_sqrm2_dat1,
                 pikm2_dat1, Lambdam2_dat1, betahatm2_dat1);

    Lq_iter_dat2(t) = Lq_func(sjk_sqrm2_dat2, mujkm2_dat2, alphajkm2_dat2, vk_sqrm2_dat2,
                 pikm2_dat2, Lambdam2_dat2, betahatm2_dat2);
  }

  return List::create(Named("alphajkm2_dat1") = alphajkm2_dat1,
                      Named("alphajkm2_dat2") = alphajkm2_dat2,
                      Named("Lq_iter_dat1") = Lq_iter_dat1,
                      Named("Lq_iter_dat2") = Lq_iter_dat2,
                      Named("mujkm2_dat1") = mujkm2_dat1,
                      Named("mujkm2_dat2") = mujkm2_dat2,
                      Named("X_full_dat1") = X_full_dat1,
                      Named("X_full_dat2") = X_full_dat2,
                      Named("sjk_sqrm2_dat1") = sjk_sqrm2_dat1,
                      Named("sjk_sqrm2_dat2") = sjk_sqrm2_dat2,
                      Named("pikm2_dat1") = pikm2_dat1,
                      Named("pikm2_dat2") = pikm2_dat2,
                      Named("vk_sqrm2_dat1") = vk_sqrm2_dat1,
                      Named("vk_sqrm2_dat2") = vk_sqrm2_dat2);
}
