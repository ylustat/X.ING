args <- commandArgs(trailingOnly = T)
# args <- c(300, 300,60,0.7,0.7,
#           0.5,0.5,0.5,0.9,0,
#           12,3,7,3,0,
#           0.5,0,1,2,2,
#           2,0.1,0,0)
n1 = as.numeric(args[1]) # sample size of data1
n2 = as.numeric(args[2]) # sample size of data2
g = as.numeric(args[3]) # number of genes (M = g * 50)
r11 = as.numeric(args[4]) # correlation within data 1
r12 = as.numeric(args[5]) # correlation within data 2
r21 = as.numeric(args[6]) # correlation 1 across data 1 and data 2
r22 = as.numeric(args[7]) # correlation 2 across data 1 and data 2
h = sqrt(as.numeric(args[8])) # input heritability (h2)
threshold = as.numeric(args[9]) # logit threshold
ld = as.numeric(args[10]) # LD correlation among SNPs
k11 = as.numeric(args[11])
k12 = as.numeric(args[12])
k21 = as.numeric(args[13])
k22 = as.numeric(args[14])
sigmaerror = as.numeric(args[15])
h_second = sqrt(as.numeric(args[16]))
prop_uncorrelated = as.numeric(args[17])
structured = as.numeric(args[18]) # 1: true, 0: false
cc0 = as.numeric(args[19])
pc1 = as.numeric(args[20])
pc2 = as.numeric(args[21])
pi_init = as.numeric(args[22])
proportion_ld = as.numeric(args[23])
ld_second = as.numeric(args[24])

print(Sys.time())
library(ggplot2)
library(tidyverse)
library(plotROC)
library(RColorBrewer)
library(patchwork)
setwd("/scratch/t.phs.yihaolu/XING_NG/simulation")
source('NG_simulation_formal_function_cri.R')
source('../FDR_func.R')
#==================================================================#
# parameter setting
## dimension parameter
M <- numeric() # number of marginal tests
N <- 2 # number of studies
K <- c(k11+k12,k21+k22) # number of conditions for each study

## sample size parameter
n <- c(n1, n2) # sample size of genotype
l <- c(10, 10) # number of blocks for the cis-SNPs of each gene
k <- c(50, 50) # number of SNPs in each block
g <- c(g, g) # number of genes

M <- unique(g * l * k) # number of marginal tests in each condition

## correlation parameter
r <- ld # indepedent: 0, moderate LD: 0.4, high LD: 0.8
r_in_data <- c(r11, r12) # co-regulation correlation across tissues
r_cross_data <- r21 # co-regulation correlation across omics

R <- cor_mat_gen(r1 = rep(r_in_data,2), r2 = r_cross_data, cut = c(k11,k21),K=K)
R

## heritability
h = list(data1 = rep(c(h,h_second),times = c(k11,k12)),
         data2 = rep(c(h,h_second),times = c(k21,k22)))


gg1 <- gg2 <- gg3 <- gg4 <- list()


#==================================================================#
# data generation
auc_list <- list()
for (repind in 1:100) { # 5
  cat("This is the ", repind, " time\n")
  token <- paste0(c(args,repind),collapse = "_")
  ## generate latent status
  set.seed(repind*round(sum(as.numeric(args))))
  X <- latent_X_gen(M,K,R)
  gamma <- latent_gamma_gen(X,threshold = threshold)

  ## generate the phenotype
  phenotype <- rep(list(rep(list(NA),g[1])),N)
  genotype <- rep(list(rep(list(NA),g[1])),N)
  beta <- rep(list(rep(list(NA),g[1])),N)
  for (ni in 1:N) {
    cis_size <- l[ni] * k[ni]
    for (i in 1:g[ni]) {
      set.seed(i)
      # genotype
      if(proportion_ld > 0) {
        l1 <- round(l[ni] * (1 - proportion_ld))
        l2 <- l[ni] - l1
        geno1 <- genotype_gen(n[ni], block = l1, block_size = k[ni], r)
        geno2 <- genotype_gen(n[ni], block = l2, block_size = k[ni], ld_second)
        genotype[[ni]][[i]] <- cbind(geno1,geno2)
      } else {
        genotype[[ni]][[i]] <- genotype_gen(n[ni], block = l[ni], block_size = k[ni], r)
      }
      genotype[[ni]][[i]] <- scale(genotype[[ni]][[i]])
      # coefficient
      num_signal <- sqrt(mean(colSums(gamma[[ni]][((i-1) * cis_size + 1):(i * cis_size),])))

      if(structured) {
        dimindex <- c(0,cumsum(K))
        beta[[ni]][[i]] <- structured_beta_gen(cis_size,K[ni],h[[ni]]/num_signal,
                                               R[(dimindex[ni]+1):(dimindex[ni+1]),(dimindex[ni]+1):(dimindex[ni+1])]) # self-defined heritability
      } else {
        beta[[ni]][[i]] <- beta_gen(cis_size,K[ni],h[[ni]]/num_signal) # self-defined heritability
      }
      # error term
      epsilon <- mvtnorm::rmvnorm(n[ni],
                                  mean = rep(0,K[ni]),
                                  sigma = diag(1-h[[ni]]^2))
      # generate phenotype
      phenotype[[ni]][[i]] <- genotype[[ni]][[i]] %*%
        (beta[[ni]][[i]] * gamma[[ni]][((i-1) * cis_size + 1):(i * cis_size),])+
        epsilon
    }
  }
  ## generate the summary statistics
  Z <- list(matrix(NA,nrow=M,ncol=K[1]),
            matrix(NA,nrow=M,ncol=K[2]))
  beta_est <- se_est <- pval_est <-
    list(matrix(NA,nrow=M,ncol=K[1]),
         matrix(NA,nrow=M,ncol=K[2]))
  for (ni in 1:N) {
    cis_size <- l[ni] * k[ni]
    for (gi in 1:g[ni]) { # gene
      print(gi)
      for (Mi in 1:cis_size) { # snp
        for (ti in 1:K[ni]) { # tissue
          x = genotype[[ni]][[gi]][,Mi]
          y = phenotype[[ni]][[gi]][,ti]
          model <- simplelm(x = x, y = y)
          Z[[ni]][(gi-1)*cis_size+Mi,ti] <- model[3]
          beta_est[[ni]][(gi-1)*cis_size+Mi,ti] <- model[1]
          se_est[[ni]][(gi-1)*cis_size+Mi,ti] <- model[2]
          pval_est[[ni]][(gi-1)*cis_size+Mi,ti] <- model[4]
        }
      }
    }
  }

  ## get true beta
  beta <- lapply(beta, function(x) do.call("rbind",x))
  true_beta <- lapply(1:N, function(x) {
    beta[[x]] * gamma[[x]]
  })

  # shuffle
  if(prop_uncorrelated > 0){
    shuffle <- round(M * prop_uncorrelated)
    shuffle_index <- lapply(1:ncol(Z[[1]]),function(x) sample(1:shuffle,shuffle))
    for (i in 1:ncol(Z[[1]])) {
      Z[[1]][1:shuffle,i] <- Z[[1]][shuffle_index[[i]],i]
      gamma[[1]][1:shuffle,i] <- gamma[[1]][shuffle_index[[i]],i]
    }
  }

  ## run our algorithm
  library(MASS)
  source("../CCA.R")
  library(pROC)
  source('../XING.R')
  print("Start Alg1 ... data 1")
  cov_data1 <- solve(cov(Z[[1]]))
  cov_data2 <- solve(cov(Z[[2]]))
  # pi_init = 0.01
  Alg1_data1 <- Alg1(betahat = Z[[1]],
                     Lambda = cov_data1,
                     iterT = 10,
                     bound=1e-8, pi_init=pi_init,vk_init=0.1)
  print("Start Alg1 ... data 2")
  Alg1_data2 <- Alg1(betahat = Z[[2]],
                     Lambda = cov_data2,
                     iterT = 10,
                     bound=1e-8, pi_init=pi_init,vk_init=0.1)

  # pc1 <- 2; pc2 <- 2
  Alg2_data1 <- Alg2_ind1(betahat = Z[[1]],
                          Lambda = cov_data1,
                          PC = pc1,
                          results_alg1 = Alg1_data1,
                          eps_thresh=1e-2,
                          iterT=5)
  Alg2_data2 <- Alg2_ind1(betahat = Z[[2]],
                          Lambda = cov_data2,
                          PC = pc2,
                          results_alg1 = Alg1_data2,
                          eps_thresh=1e-2,
                          iterT=5)


  print("Start Alg4 ...")

  Alg4_mQTL_eQTL <- Alg4(betahat1 = Z[[1]],
                         betahat2 = Z[[2]],
                         Lambda1 = cov_data1,
                         Lambda2 = cov_data2,
                         CC = cc0,
                         PC1 = pc1,
                         PC2 = pc2,
                         results_alg1_dat1=Alg1_data1,
                         results_alg1_dat2=Alg1_data2,
                         results_alg2_dat1 = Alg2_data1,
                         results_alg2_dat2 = Alg2_data2,
                         iterT = 5,
                         bound=1e-4,
                         sparse = F)

  res_list <-  list(beta_est[[1]],beta_est[[2]],Z[[1]],Z[[2]],Alg2_data1,Alg2_data2,Alg4_mQTL_eQTL)
  save(res = res_list, file = "pp.RData")

  auc_list[[repind]] <- c(auc(c(gamma[[1]]), c(Alg4_mQTL_eQTL$alphajkm2_dat1),quiet=T),
                          auc(c(gamma[[2]]), c(Alg4_mQTL_eQTL$alphajkm2_dat2),quiet=T))

  #==================================================================#
  ## comparison method
  ### mashr
  library(ashr)
  library(mashr)
  data_mashr1 <- mash_set_data(beta_est[[1]], se_est[[1]])
  U.1 = cov_canonical(data_mashr1)
  m.1 = mash(data_mashr1, U.1)
  data_mashr2 <- mash_set_data(beta_est[[2]], se_est[[2]])
  U.2 = cov_canonical(data_mashr2)
  m.2 = mash(data_mashr2, U.2)
  auc(c(gamma[[1]]), c(1-m.1$result$lfdr),quiet=T)
  auc(c(gamma[[2]]), c(1-m.2$result$lfdr),quiet=T)

  ### metasoft
  library(data.table)
  meta_input_data1 <- matrix(NA, nrow = M, ncol = K[1] * 2)
  meta_input_data2 <- matrix(NA, nrow = M, ncol = K[2] * 2)
  meta_input_data1[,2*(1:K[1])-1] <- beta_est[[1]]
  meta_input_data1[,2*(1:K[1])] <- se_est[[1]]
  meta_input_data2[,2*(1:K[2])-1] <- beta_est[[2]]
  meta_input_data2[,2*(1:K[2])] <- se_est[[2]]
  meta_input_data1 <- data.frame(SNP = paste0("rs",1:M),meta_input_data1)
  meta_input_data2 <- data.frame(SNP = paste0("rs",1:M),meta_input_data2)
  fwrite(meta_input_data1, file = paste0("/scratch/t.phs.yihaolu/LLR/simulation/Metasoft/input/",token,"input_data1.txt"), sep = "\t", col.names = F)
  fwrite(meta_input_data2, file = paste0("/scratch/t.phs.yihaolu/LLR/simulation/Metasoft/input/",token,"input_data2.txt"), sep = "\t", col.names = F)

  system(paste0("cd /scratch/t.phs.yihaolu/LLR/simulation/Metasoft/;java -jar Metasoft.jar -input ../Metasoft/input/",
                token,"input_data1.txt -mvalue true -mvalue_p_thres 1 -output ../Metasoft/input/",token,"_data1"))
  result_data1 <- fread(paste0("/scratch/t.phs.yihaolu/LLR/simulation/Metasoft/input/",token,"_data1"))
  pp1 <- data.matrix(result_data1[,(17+K[1]):(17+2*K[1]-1)]);pp1[is.na(pp1)] <- 0

  system(paste0("cd /scratch/t.phs.yihaolu/LLR/simulation/Metasoft/;java -jar Metasoft.jar -input ../Metasoft/input/",
                token,"input_data2.txt -mvalue true -mvalue_p_thres 1 -output ../Metasoft/input/",token,"_data2"))
  result_data2 <- fread(paste0("/scratch/t.phs.yihaolu/LLR/simulation/Metasoft/input/",token,"_data2"))
  pp2 <- data.matrix(result_data2[,(17+K[2]):(17+2*K[2]-1)]);pp2[is.na(pp2)] <- 0

  ## paintor
  setwd("/scratch/t.phs.yihaolu/LLR/simulation/PAINTOR_V3.0/")
  for (n0 in 1:2) {
    for (i in 1:(l[n0]*g[n0])) {
      # print(i)
      locus <- Z[[n0]][((i-1) * k[n0] + 1) : (i * k[n0]),]
      colnames(locus) <- paste0("zscore",1:K[n0])
      ld <- diag(1,k[n0])
      locus <- round(locus,4)
      fwrite(locus, file = paste0("inputfile/",token,"locus",i),sep = " ",quote = F,row.names = F)
      fwrite(ld, file = paste0("inputfile/",token,"locus",i,".ld"),sep = " ",quote = F,row.names = F,col.names = F)
      annotation <- cbind(rep(1,k[n0]),rep(0,k[n0]));colnames(annotation) <- c("Coding","DHS")
      fwrite(annotation, file = paste0("inputfile/",token,"locus",i,".annotations"),sep = " ",quote = F,row.names = F)
    }
    filelist <- paste0(token,"locus",1:(l[n0]*g[n0])) #
    write.table(filelist,file = paste0("inputfile/",token,"input.files"),row.names = F,col.names = F,quote = F)

    for (i in 1:K[n0]) {
      # print(i)
      system(paste0("./PAINTOR -input inputfile/",token,"input.files -in inputfile/ -out inputfile/ -Zhead zscore",i," -LDname ld -enumerate 2"))
      system(paste0("cat ",paste0(paste0("inputfile/",token,"locus",1:(l[n0]*g[n0]),".results",collapse = " ")), " > inputfile/",token,"data.",n0,".all.tissue.",i))
      system(paste0("rm inputfile/",token,"locus*.results"))
    }
    system(paste0("rm inputfile/",token,"locus*"))
    system(paste0("rm inputfile/Log*"))
    system(paste0("rm inputfile/Enrichment.Values inputfile/",token,"input.files"))

    post_prob <- list()
    for (i in 1:K[n0]) {
      result <- fread(paste0("inputfile/",token,"data.",n0,".all.tissue.",i))
      post_prob[[i]] <- result$Posterior_Prob[result$Posterior_Prob!="Posterior_Prob"]
    }
    post_prob <- do.call("cbind",post_prob) %>% data.frame()
    fwrite(post_prob, file = paste0("inputfile/",token,"post.prob.data.",n0,".all.tissue"),sep = "\t",quote = F,row.names = F,col.names = F)
  }

  system(paste0("rm inputfile/",token,"data.*.all.tissue*"))

  paintor_pp1 <- fread(paste0("inputfile/",token,"post.prob.data.1.all.tissue")) %>% data.matrix()
  paintor_pp2 <- fread(paste0("inputfile/",token,"post.prob.data.2.all.tissue")) %>% data.matrix()

  ### HT-eQTL
  setwd("/scratch/t.phs.yihaolu/XING_NG/simulation")
  source('../HT_eQTL.R')
  Z <- lapply(Z, function(x) {
    x[x>25] <- 25
    x[x<(-25)] <- (-25)
    return(x)
  })
  system.time({ht_pp1 <- HT_eQTL(Z[[1]])})
  system.time({ht_pp2 <- HT_eQTL(Z[[2]])})

  # we do not run mash, paintor, metasoft
  ht_pp1 <- Alg4_mQTL_eQTL$alphajkm2_dat1
  ht_pp2 <- Alg4_mQTL_eQTL$alphajkm2_dat2

  paintor_pp1 <- ht_pp1 #+ matrix(runif(nrow(ht_pp1)*ncol(ht_pp1),-1e-3,1e-3), nrow = nrow(ht_pp1), ncol = ncol(ht_pp1))
  paintor_pp2 <- ht_pp2 #+ matrix(runif(nrow(ht_pp2)*ncol(ht_pp2),-1e-3,1e-3), nrow = nrow(ht_pp2), ncol = ncol(ht_pp2))

  pp1 <- ht_pp1
  pp2 <- ht_pp2

  m.1 <- m.2 <- list(list())
  m.1$result$lfdr <- ht_pp1
  m.2$result$lfdr <- ht_pp2
  #==================================================================#
  ## result evaluation

  dynamic.index1 <- which(rowSums(gamma[[1]])>=2 & rowSums(gamma[[1]])<=5)
  dynamic.index2 <- which(rowSums(gamma[[2]])>=2 & rowSums(gamma[[2]])<=5)

  auc_df <- data.frame(data1 = c(auc(c(gamma[[1]]), c(Alg1_data1$alphajkm1),quiet=T),
                                 auc(c(gamma[[1]]), c(Alg2_data1$alphajkm1),quiet=T),
                                 auc(c(gamma[[1]]), c(Alg4_mQTL_eQTL$alphajkm2_dat1),quiet=T), # LLR
                                 auc(c(gamma[[1]]), c(1-m.1$result$lfdr),quiet=T), # mashr
                                 auc(c(gamma[[1]]), c(pp1),quiet=T), # metasoft
                                 auc(c(gamma[[1]]), c(paintor_pp1),quiet=T), # paintor
                                 auc(c(gamma[[1]]), c(ht_pp1),quiet=T)),
                       data2 = c(auc(c(gamma[[2]]), c(Alg1_data2$alphajkm1),quiet=T),
                                 auc(c(gamma[[2]]), c(Alg2_data2$alphajkm1),quiet=T),
                                 auc(c(gamma[[2]]), c(Alg4_mQTL_eQTL$alphajkm2_dat2),quiet=T), # LLR
                                 auc(c(gamma[[2]]), c(1-m.2$result$lfdr),quiet=T), # mashr
                                 auc(c(gamma[[2]]), c(pp2),quiet=T), # metasoft
                                 auc(c(gamma[[2]]), c(paintor_pp2),quiet=T), # metasoft
                                 auc(c(gamma[[2]]), c(ht_pp2),quiet=T)),
                       data1.dynamic=c(auc(c(gamma[[1]][dynamic.index1,]),
                                           c(Alg1_data1$alphajkm1[dynamic.index1,]),quiet=T),
                                       auc(c(gamma[[1]][dynamic.index1,]),
                                           c(Alg2_data1$alphajkm1[dynamic.index1,]),quiet=T),
                                       auc(c(gamma[[1]][dynamic.index1,]),
                                           c(Alg4_mQTL_eQTL$alphajkm2_dat1[dynamic.index1,]),quiet=T), # LLR
                                       auc(c(gamma[[1]][dynamic.index1,]),
                                           c(1-m.1$result$lfdr[dynamic.index1,]),quiet=T), # mashr
                                       auc(c(gamma[[1]][dynamic.index1,]),
                                           c(pp1[dynamic.index1,]),quiet=T), # metasoft
                                       auc(c(gamma[[1]][dynamic.index1,]),
                                           c(paintor_pp1[dynamic.index1,]),quiet=T), # paintor
                                       auc(c(gamma[[1]][dynamic.index1,]),
                                           c(ht_pp1[dynamic.index1,]),quiet=T)),
                       data2.dynamic=c(auc(c(gamma[[2]][dynamic.index1,]),
                                           c(Alg1_data2$alphajkm1[dynamic.index1,]),quiet=T),
                                       auc(c(gamma[[2]][dynamic.index1,]),
                                           c(Alg2_data2$alphajkm1[dynamic.index1,]),quiet=T),
                                       auc(c(gamma[[2]][dynamic.index1,]),
                                           c(Alg4_mQTL_eQTL$alphajkm2_dat2[dynamic.index1,]),quiet=T), # LLR
                                       auc(c(gamma[[2]][dynamic.index1,]),
                                           c(1-m.2$result$lfdr[dynamic.index1,]),quiet=T), # mashr
                                       auc(c(gamma[[2]][dynamic.index1,]),
                                           c(pp2[dynamic.index1,]),quiet=T), # metasoft
                                       auc(c(gamma[[2]][dynamic.index1,]),
                                           c(paintor_pp2[dynamic.index1,]),quiet=T), # metasoft
                                       auc(c(gamma[[2]][dynamic.index1,]),
                                           c(ht_pp2[dynamic.index1,]),quiet=T)))
  # rownames(auc_df) <- c("Marginal_LLR1","Marginal_LLR2","LLR","mash","metasoft","ht-eQTL")
  # auc_df

  rmse_df <- data.frame(Method = c("LLR","mash","Marginal model"),
                        data1 = c(rmse(true_beta[[1]],(Alg4_mQTL_eQTL$mujkm2_dat1 * se_est[[1]])),
                                  rmse(true_beta[[1]],m.1$result$PosteriorMean),
                                  rmse(true_beta[[1]],beta_est[[1]])),
                        data2 = c(rmse(true_beta[[2]],(Alg4_mQTL_eQTL$mujkm2_dat2 * se_est[[2]])),
                                  rmse(true_beta[[2]],m.2$result$PosteriorMean),
                                  rmse(true_beta[[2]],beta_est[[2]])),
                        data1.true = c(rmse(true_beta[[1]][gamma[[1]]!=0],(Alg4_mQTL_eQTL$mujkm2_dat1 * se_est[[1]])[gamma[[1]]!=0]),
                                       rmse(true_beta[[1]][gamma[[1]]!=0],m.1$result$PosteriorMean[gamma[[1]]!=0]),
                                       rmse(true_beta[[1]][gamma[[1]]!=0],beta_est[[1]][gamma[[1]]!=0])),
                        data2.true = c(rmse(true_beta[[2]][gamma[[2]]!=0],(Alg4_mQTL_eQTL$mujkm2_dat2 * se_est[[2]])[gamma[[2]]!=0]),
                                       rmse(true_beta[[2]][gamma[[2]]!=0],m.2$result$PosteriorMean[gamma[[2]]!=0]),
                                       rmse(true_beta[[2]][gamma[[2]]!=0],beta_est[[2]][gamma[[2]]!=0])),
                        data1.dynamic = c(rmse(true_beta[[1]][dynamic.index1,],(Alg4_mQTL_eQTL$mujkm2_dat1 * se_est[[1]])[dynamic.index1,]),
                                          rmse(true_beta[[1]][dynamic.index1,],m.1$result$PosteriorMean[dynamic.index1,]),
                                          rmse(true_beta[[1]][dynamic.index1,],beta_est[[1]][dynamic.index1,])),
                        data2.dynamic = c(rmse(true_beta[[2]][dynamic.index2,],(Alg4_mQTL_eQTL$mujkm2_dat2 * se_est[[2]])[dynamic.index2,]),
                                          rmse(true_beta[[2]][dynamic.index2,],m.2$result$PosteriorMean[dynamic.index2,]),
                                          rmse(true_beta[[2]][dynamic.index2,],beta_est[[2]][dynamic.index2,])))
  result <- list(parameter = args, auc = auc_df, rmse = rmse_df)
  save(result, file = paste0("result/1_13_2022/",token,".RData"))
}
