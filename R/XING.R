#' @title Cross INtegrative Genomics method (X-ING)
#'
#' @description X-ING models and estimates the latent binary association status of each statistic, and captures the omic-shared, context-shared, omic-specific and context-specific major patterns in a hierarchical Bayesian model.
#'
#' @details
#' This function \code{XING}, jointly estimates the posterior means of input statistics, models the latent binary association probabilities based on the posterior means, and outputs the posterior mean and probability of association for each input statistic.
#'
#' @param z_list a list of association summary statistics, each with multivariate contexts.
#' @param Lambda_list a list of tissue-tissue correlation matrix among all tissues due to potential sample overlap.
#' @param CC number of canonical components to keep in canonical correlation analysis (CCA).
#' @param PC_list numbers of principal components to keep in principal correlation analysis (PCA).
#' @param iterT maximum iteration times.
#' @param method methods to extract shared patterns.
#' @param tolerance tolerance for iteration of generalized canonical correlation analysis.
#' @export
#' @return a list with following elements:
#' \item{mu}{posterior means of input statistics.}
#' \item{alpha}{posterior probability of association for each input statistic.}
#' @examples
#' library(MASS)
#' if (!require("CCA")) install.packages("CCA")
#' if (!require("RGCCA")) install.packages("RGCCA")
#' library(CCA)
#' library(RGCCA)
#' data(example)
#' z_list <- lapply(example,function(x) x[[1]])
#' res <- XING(z_list = z_list[1:2])
#' res <- XING(z_list = z_list)

XING <- function(z_list, Lambda_list = NULL, CC = 2,
                 PC_list = NULL, iterT = 20, method = "CCAPCA",
                 tolerance = 1e-8){
  num_data <- length(z_list)
  M <- unique(sapply(z_list,nrow))
  if(length(M) > 1) stop("All data should have same number of rows!")
  if(length(Lambda_list) == 0) Lambda_list <- lapply(z_list,function(x) solve(cov(x)))
  if(length(PC_list) == 0) PC_list <- rep(2,length(z_list))
  K <- sapply(z_list, ncol)

  # Fit the starting model first
  XING_single_result <- mapply(function(z, Lambda)
    XING_single(z = z,
                Lambda = Lambda,
                iterT = 20, eps_thres = 1e-3),z_list,Lambda_list,SIMPLIFY = F)

  # Defining and initializing new variables using results from starting model
  s2 <- lapply(XING_single_result, function(x) x$s2)
  mu <- lapply(XING_single_result, function(x) x$mu)

  # Posterior for association probability
  alpha <- lapply(XING_single_result, function(x) x$alpha)

  # Parameters
  v2 <- lapply(XING_single_result, function(x) x$v2)
  pi <- lapply(XING_single_result, function(x) x$alpha)
  x0k <- lapply(XING_single_result, function(x) log(x$pi/(1-x$pi)))

  ## EM algorithm with variational approximation:
  # Vector for storing lower bound to be maximized
  Lq_iter <- rep(list(c(-Inf,rep(NA, iterT-1))),num_data)
  uu2 <- X_est <- list(); bound = 1e-8
  upbound <- 1-bound;lowbound <- bound
  for(t in 2:iterT){
    print(paste0("Iteration times: ",t))
    ## E step
    for (i in 1:num_data) {
      for(k in 1:K[i]){
        s2[[i]][,k] <- 1/(diag(Lambda_list[[i]])[k] + (1/v2[[i]][k]))
      }
      mu[[i]] <- (z_list[[i]] %*% Lambda_list[[i]] -
                    ((alpha[[i]]*mu[[i]]) %*% Lambda_list[[i]] -
                       (alpha[[i]]*mu[[i]]) %*% diag(diag(Lambda_list[[i]]))))%*%
        (diag((diag(Lambda_list[[i]]) + (1/v2[[i]]))^(-1)))
      uu2[[i]] <- log(pi[[i]]/(1-pi[[i]])) +
        0.5*((log(s2[[i]])- matrix(1,M,K[i]) %*%
                diag(log(v2[[i]]))) + mu[[i]]^2/s2[[i]])
      alpha[[i]] <- 1/(1+exp(-uu2[[i]]))
      alpha[[i]] <- ifelse(alpha[[i]]>upbound,upbound, alpha[[i]])
      alpha[[i]] <- ifelse(alpha[[i]]<lowbound,lowbound, alpha[[i]])
    }
    ## M step
    ## Update X
    X_full <- mapply(function(alpha, x0) log(alpha) - log(1-alpha) - x0, alpha, x0k, SIMPLIFY = F)
    sd_list <- lapply(X_full,function(x) apply(x,2,sd))
    mean_list <- lapply(X_full,function(x) apply(x,2,mean))
    X_full <- lapply(X_full,scale)
    for (i in 1:num_data) {
      for(k in 1:K[i]){
        v2[[i]][k] <- sum(alpha[[i]][,k]*(s2[[i]][,k] + mu[[i]][,k]^2))/sum(alpha[[i]][,k])
      }
    }

    if(method %in% c("CCA","CCAPCA")) {
      # Use CCA first to extract shared patterns across data types
      if(num_data == 2) { # simple CCA case
        ccas <- cc(X_full[[1]], X_full[[2]])
        Ahat <- ccas$xcoef
        Bhat <- ccas$ycoef
        U <- X_full[[1]] %*% Ahat
        V <- X_full[[2]] %*% Bhat
        X_est[[1]] <- U[,1:CC] %*% ginv(Ahat[,1:CC])
        X_est[[2]] <- V[,1:CC] %*% ginv(Bhat[,1:CC])
        # predicted probability using both CCA component and individual component
        pi[[1]] <- 1/(1+exp(-X_est[[1]] - x0k[[1]]))
        pi[[2]] <- 1/(1+exp(-X_est[[2]] - x0k[[2]]))
      } else if (num_data > 2) {
        data_list <- c(X_full, list(do.call("cbind", X_full)))
        C <- matrix(0,length(data_list),length(data_list))
        C[nrow(C),1:(ncol(C)-1)] <- C[1:(nrow(C)-1),ncol(C)] <- 1
        options(warn=-1)
        cca.with.rgcca = rgcca(A = data_list,
                               C = C,
                               ncomp = rep(CC,length(data_list)),
                               tau = rep(0,length(data_list)),
                               verbose = F,
                               tol = tolerance)
        options(warn=0)
        U <- cca.with.rgcca$Y[1:num_data]
        Ahat <- mapply(function(x, u) ginv(x) %*% u, X_full, U, SIMPLIFY = F)
        X_est <- mapply(function(u, a) u[,1:CC] %*% ginv(a[,1:CC]), U, Ahat, SIMPLIFY = F)
      }

      # Use PCA to extract shared patterns across contexts
      if(method == "CCAPCA") {
        X_res <- mapply(function(xfull,xest) xfull - xest, X_full, X_est, SIMPLIFY = F)
        X_res <- lapply(X_res, function(x) {x[abs(x) > 100] <- sign(x[abs(x) > 100]) * 100; x}) # to avoid problem in SVD
        X_svd_res <- lapply(X_res, function(x) svd(scale(x)))
        X_red_res <- mapply(function(xsvd, pc) xsvd$u[,1:pc] %*% diag(xsvd$d[1:pc]) %*% t(xsvd$v[,1:pc]), X_svd_res, PC_list, SIMPLIFY = F)
        pi <- mapply(function(xest, xred, x0) 1/(1+exp(-xest - xred - x0)), X_est, X_red_res, x0k, SIMPLIFY = F)  # predicted probability using both CCA component and individual component
      } else {
        pi <- mapply(function(xest,x0) 1/(1+exp(-xest - x0)), X_est, x0k, SIMPLIFY = F)
      }
    }

    pi <- lapply(pi, function(x) {
      x <- ifelse(x > upbound, upbound, x)
      x <- ifelse(x < lowbound, lowbound, x)
      x
    })
    # Compute lower bound
    for (i in 1:num_data) {
      Lq_iter[[i]][t] <- Lq(s2 = s2[[i]], mu = mu[[i]],
                            alpha =alpha[[i]], v2 = v2[[i]],
                            pi = pi[[i]], Lambda = Lambda_list[[i]],
                            z = z_list[[i]])
    }
    if(t > 2 & prod(sapply(Lq_iter, function(x) (x[t] - x[t-1])/x[t-1] < 1e-4))==1) break

  }
  return(list(mu = mu,
              alpha = alpha))
}

