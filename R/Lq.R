#' Q function lower bound in variational approximation that needs to be maximized in EM

Lq <- function(s2, mu, alpha, v2, pi, Lambda, z){
  M <- nrow(z)
  K <- ncol(z)
  inner1 <- rep(NA, M)
  inner2 <- rep(NA, M)
  inner3 <- rep(NA, M)
  inner4 <- rep(NA, M)

  B <- z - alpha*mu

  one_mat <- matrix(1,M,K)

  inner1 = rowSums((B%*%Lambda)*(B))
  inner2 = rowSums((alpha*(mu^2 + s2)- (alpha^2)*(mu^2))%*%diag(diag(Lambda)))
  inner3 = rowSums(alpha*(1 + (log(s2)-one_mat%*%diag(log(v2))) - ((s2 + mu^2)%*%diag(v2^(-1)))))
  inner4 = rowSums(alpha*(log(alpha)- log(pi)) + (1-alpha)*(log(1-alpha)- log(1-pi)))

  Lq <- -0.5*sum(inner1) - 0.5*sum(inner2) + 0.5*sum(inner3) - sum(inner4) + (M/2)*log(det(Lambda)) + 0.5*K*M
  return(Lq)
}

