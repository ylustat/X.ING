#'
#' @title Starting model for X-ING
#'
#' @description \code{X-ING} starting model without considering shared patterns across contexts and data types.
#'
#' @details
#' This function \code{XING_single}, X-ING estimates the posterior means of input statistics, models the latent binary association probabilities based on the posterior means, and outputs the posterior mean and probability of association for each input statistic without considering shared patterns across contexts and data types.
#'
#' @param z a matrix of association summary statistics with multivariate contexts.
#' @param Lambda a matrix of tissue-tissue correlation matrix among all tissues due to potential sample overlap.
#' @param iterT maximum iteration times.
#' @param eps_thres: stopping rule threshold.
#' @param vk_init: initial value for the variance of latent genetic association distribution.
#' @param pi_init: initial probability for having latent genetic associations.
#' @export

XING_single <- function(z, Lambda, iterT = 20, eps_thres = 1e-3,
                 vk_init = 0.1,bound = 1e-4, pi_init = 1e-4){
  M <- nrow(z)
  K <- ncol(z)
  # Defining and initializing new variables;
  s2 <- matrix(1, nrow = M, ncol = K) # posterior variance for beta
  mu <- matrix(0, nrow = M, ncol = K) # posterior mean for beta

  # Posterior for association probability
  alpha <- matrix(pi_init, nrow = M, ncol = K)

  # Variational parameters
  v2 <- rep(vk_init, K)
  pi <- matrix(pi_init, nrow = M, ncol = K, byrow = T)

  upbound <- 1-bound
  lowbound <- bound

  # EM algorithm with variational approximation:
  uu <- matrix(NA, nrow = M, ncol = K)
  # Vector for storing lower bound to be maximized
  Lq_iter <- rep(NA, iterT)
  Lq_iter[1] <- -Inf
  for(t in 2:iterT){
    ## E step
    for (k in 1:K) {
      s2[,k] <- 1/(diag(Lambda)[k] + (1/v2[k]))
    }
    mu = (z%*%Lambda - ((alpha*mu)%*%Lambda - (alpha*mu)%*%diag(diag(Lambda))))%*%(diag((diag(Lambda) + (v2)^(-1))^(-1)))
    uu = log(pi/(1-pi)) + 0.5*((log(s2) - matrix(1,M,K)%*%diag(log(v2))) + mu^2/s2)
    alpha <- 1/(1+exp(-uu))
    alpha <- ifelse(alpha > upbound, upbound, alpha)
    alpha <- ifelse(alpha < lowbound, lowbound, alpha)

    ## M step
    for(k in 1:K){
      v2[k] <- sum(alpha[,k]*(s2[,k] + mu[,k]^2))/sum(alpha[,k])
    }
    for(k in 1:K){
      pi[,k] <- rep(colMeans(alpha)[k], M)
      pi[,k] <- rep(ifelse(pi[1,k] > upbound, upbound, pi[1,k]), M)
      pi[,k] <- rep(ifelse(pi[1,k] < lowbound, lowbound, pi[1,k]), M)
    }

    Lq_iter[t] <- Lq(s2 = s2,
                     mu=mu,
                     alpha=alpha,
                     v2 = v2,
                     pi = pi,
                     Lambda = Lambda,
                     z = z)

    if(Lq_iter[t] - Lq_iter[t-1] < eps_thres){
      break
    }
  }
  # Output posterior mean and probability
  return(list(alpha = alpha, mu = mu,
              gamma = gamma, s2 = s2,
              z = z, Lambda = Lambda,
              v2 = v2, pi = pi, Lq_iter=Lq_iter))
}

