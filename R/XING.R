XING <- function(z_list, Lambda_list = NULL, CC = 2,
                 PC_list = NULL, iterT = 20,
                 tolerance = 1e-4){
  L <- length(z_list)
  if(length(PC_list) == 0) PC_list <- rep(2,length(z_list))
  if(length(Lambda_list) == 0) Lambda_list <- lapply(z_list,function(x) solve(cov(x)))
  if(L == 2) {
    res <- XING_two_data(z_list[[1]],z_list[[2]],Lambda_list[[1]],Lambda_list[[2]],
                  CC, PC_list[1], PC_list[2],
                  iterT = iterT)
    res <- list(mu = list(res$mu1,res$mu2),
                alpha = list(res$alpha1,res$alpha2))
  } else {
    res <- XING_multidata(z_list = z_list, Lambda_list = Lambda_list, CC = CC,
                   PC_list = PC_list, iterT = iterT, 
                   tolerance = tolerance)
  }
  return(res)
}