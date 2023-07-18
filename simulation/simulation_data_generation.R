# functions for simulation_formal.R
genotype_gen <- function(sample_size, block, block_size, r) {
  LD_corr <- ar1_cor(block_size, r) # generate the LD correlation matrix within each block
  genotype_list <- list()
  for (j in 1:block) {
    genotype_cont <- mvtnorm::rmvnorm(sample_size, rep(0,block_size), LD_corr)
    genotype <- genotype_cont
    genotype_list[[j]] <- genotype
  }
  genotype_list <- do.call("cbind",genotype_list)
  return(genotype_list)
}

latent_X_gen <- function(M,K,R) {
  X <- mvtnorm::rmvnorm(M, mean = rep(0,sum(K)), R)
  column_index <- c(0,cumsum(K))
  X_list <- list()
  for (i in 1:length(K)) {
    X_list[[i]] <- X[,(column_index[i]+1):column_index[i+1]]
  }
  return(X_list)
}

latent_gamma_gen <- function(X,threshold = 0.9) {
  gamma <- lapply(X, function(x) {
    pi <- exp(x)/(1+exp(x))
    return(ifelse(pi > threshold, 1, 0))
  })
  return(gamma)
}


beta_gen <- function(M,K,h) {
  # beta <- mvtnorm::rmvnorm(M,mean = rep(0,K),sigma = diag(h^2))
  beta <- matrix(rnorm(M*K,0,h),nrow=M)
  return(beta)
}

structured_beta_gen <- function(M,K,h,R){
  covMat <- diag(h) %*% R %*% diag(h)
  beta <- mvtnorm::rmvnorm(M,mean = rep(0,K),sigma = covMat)
  return(beta)
}


ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, 
                         ncol = n, byrow = TRUE) - (1:n - 1))
  rho^exponent
}

cor_mat_gen <- function(r1,r2,cut,K,diag = F) {
  cut <- mapply(function(x,k) c(x,k-x),cut,K,SIMPLIFY = F) %>% unlist()
  cor_mat <- matrix(0,nrow = sum(cut),ncol = sum(cut))
  cut <- c(0,cut)
  cut <- cumsum(cut)
  for (i in 1:(length(cut)-1)) {
    cor_mat[(cut[i]+1):cut[i+1],(cut[i]+1):cut[i+1]] <- r1[i]
  }
  for (i in 1:(length(cut)-3)) {
    if(diag) {
      diag(cor_mat[(cut[i]+1):cut[i+1],(cut[i+2]+1):cut[i+3]]) <- r2
    } else {
      cor_mat[(cut[i]+1):cut[i+1],(cut[i+2]+1):cut[i+3]] <- r2
    }
  }
  cor_mat[lower.tri(cor_mat)] <- 0
  cor_mat <- cor_mat + t(cor_mat)
  diag(cor_mat) <- 1
  return(cor_mat)
}

# cor_mat_gen <- function(r1,r2,cut,K) {
#   cut <- mapply(function(x,k) c(x,k-x),cut,K,SIMPLIFY = F) %>% unlist()
#   cor_mat <- matrix(r2,nrow = sum(cut),ncol = sum(cut))
#   cut <- c(0,cut)
#   cut <- cumsum(cut)
#   for (i in 1:(length(cut)-1)) {
#     cor_mat[(cut[i]+1):cut[i+1],(cut[i]+1):cut[i+1]] <- (1 - r1[i]) * Imat(cut[i+1]-cut[i]) + r1[i] * One_mat(cut[i+1]-cut[i])
#   }
#   return(cor_mat)
# }

# heatmap.2(R,dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',key=F)
# 
# cor_mat_gen <- function(r1,r2,L,K) {
#   cut <- mapply(function(x,k) c(x,k-x),L,K,SIMPLIFY = F) %>% unlist()
#   cor_mat <- matrix(0,nrow = sum(cut),ncol = sum(cut))
#   cut <- c(0,cut)
#   cut <- cumsum(cut)
#   for (i in 1:(length(cut)-1)) {
#       cor_mat[(cut[i]+1):cut[i+1],(cut[i]+1):cut[i+1]] <- r1[i]
#   }
#   
#   for (i in seq(from = 1, to = length(cut)-3, by = 2)) {
#     for (j in seq(from = i+2, to = length(cut)-1, by = 2)) {
#       cor_mat[(cut[i]+1):cut[i+1],(cut[j]+1):cut[j+1]] <- r2
#       cor_mat[(cut[i+1]+1):cut[i+2],(cut[j+1]+1):cut[j+2]] <- r2
#     }
#   }
#   cor_mat[lower.tri(cor_mat)] <- 0
#   cor_mat <- cor_mat + t(cor_mat)
#   diag(cor_mat) <- 1
#   return(cor_mat)
# }


I_mat <- function(n) {
  diag(1,n)
}

One_mat <- function(n) {
  matrix(1,n,n)
}

simplelm <- function (x, y) {
  ## number of data
  n <- length(x)
  ## centring
  y0 <- sum(y) / length(y); yc <- y - y0
  x0 <- sum(x) / length(x); xc <- x - x0
  ## fitting an intercept-free model: yc ~ xc + 0
  xty <- c(crossprod(xc, yc))
  xtx <- c(crossprod(xc))
  slope <- xty / xtx
  rc <- yc - xc * slope
  ## Pearson estimate of residual standard error
  sigma2 <- c(crossprod(rc)) / (n - 2)
  ## standard error for slope
  slope_se <- sqrt(sigma2 / xtx)
  ## t-score and p-value for slope
  tscore <- slope / slope_se
  pvalue <- 2 * pt(abs(tscore), n - 2, lower.tail = FALSE)
  ## return estimation summary for slope
  c("Estimate" = slope, "Std. Error" = slope_se, "t value" = tscore, "Pr(>|t|)" = pvalue)
}


rmse <- function(yhat,y) {
  sqrt(mean((yhat - y)^2))
}


