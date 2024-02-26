\name{XING}
\alias{XING}
\title{EM Algorithm with Variational Approximation for Data Analysis}
\description{
  Implements an Expectation-Maximization (EM) algorithm with variational approximation for statistical analysis in a multivariate context. The function is designed to analyze two datasets simultaneously, using variational methods to optimize computational efficiency and accuracy.
}
\usage{
XING(z_list, Lambda_list, CC = 2, PC_list, iterT = 20, tolerance = 1e-4)
}
\arguments{
  \item{z_list}{A list of estimated z-scores across multiple contexts or tissues, where each element corresponds to one omics data type.}
  \item{Lambda_list}{A list of inversed estimated tissue-tissue correlation matrix among all tissues due to potential sample overlap.}
  \item{CC}{Number of canonical components to keep in canonical correlation analysis (CCA).}
  \item{PC_list}{Numbers of principal components to keep in principal correlation analysis (PCA) for each omics data type.}
  \item{iterT}{Number of iterations for the EM algorithm.}
  \item{tolerance}{The stopping value for convergence for GCCA.}
}
\value{
  Returns a list containing:
  \itemize{
    \item \code{alpha}: Posterior probabilities of associations.
    \item \code{mu}: Posterior means of input statistics.
  }
}
\details{
This function \code{XING}, jointly estimates the posterior means of input statistics, models the latent binary association probabilities based on the posterior means, and outputs the posterior mean and probability of association for each input statistic.
}
\examples{
library(MASS)
library(CCA)
library(RGCCA)
data(example_data)
z <- example$z
res <- XING(z_list = z[1:2])
res <- XING(z_list = z)
}

