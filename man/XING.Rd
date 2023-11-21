\name{XING}
\alias{XING}
\title{EM Algorithm with Variational Approximation for Data Analysis}
\description{
  Implements an Expectation-Maximization (EM) algorithm with variational approximation for statistical analysis in a multivariate context. The function is designed to analyze two datasets simultaneously, using variational methods to optimize computational efficiency and accuracy.
}
\usage{
XING(z1, z2, Lambda1, Lambda2, CC, PC1, PC2, result_starting1, result_starting2, result_single1, result_single2, iterT = 20)
}
\arguments{
  \item{z1}{Matrix of association summary statistics with multivariate contexts for the first omics.}
  \item{z2}{Matrix of association summary statistics with multivariate contexts for the second omics.}
  \item{Lambda1}{Tissue-tissue correlation matrix among all tissues due to potential sample overlap for the first omics.}
  \item{Lambda2}{Tissue-tissue correlation matrix among all tissues due to potential sample overlap for the second omics.}
  \item{CC}{Number of canonical components to keep in canonical correlation analysis (CCA).}
  \item{PC1}{Numbers of principal components to keep in principal correlation analysis (PCA) for the first data.}
  \item{PC2}{Numbers of principal components to keep in PCA for the second data}
  \item{result_starting1}{List of starting model results for the first dataset.}
  \item{result_starting2}{List of starting model results for the second dataset.}
  \item{result_single1}{List of results from a single data version of X-ING for the first dataset.}
  \item{result_single2}{List of results from a single data version of X-ING for the second dataset.}
  \item{iterT}{Number of iterations for the EM algorithm.}
}
\details{
  The function conducts an EM algorithm, an iterative process comprising the Expectation and Maximization steps. It uses variational approximation techniques to handle large datasets efficiently. The algorithm processes two datasets concurrently, leveraging their respective precision matrices and initial results to refine the analysis iteratively.
}
\value{
  Returns a list containing:
  \itemize{
    \item \code{alpha1}: Posterior probabilities of associations for the first dataset.
    \item \code{alpha2}: Posterior probabilities of associations for the second dataset.
    \item \code{mu1}: Posterior means of input statistics for the first dataset.
    \item \code{mu2}: Posterior means of input statistics for the second dataset.
  }
}
