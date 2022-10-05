# X.ING: Cross Integrative Genomics

This package implements methods for cross-feature integration of effects from multiple studies each with multivariate contexts. The latent binary association status of each statistic was modeled to capture the omics-shared and context-shared major patterns in a hierarchical Bayesian model.

## Installation

Install `X.ING` from GitHub: 

```R
library(devtools)
install_github("ylustat/X.ING")
```

Note that `X.ING` depends on the 'CCA', "RGCCA", "MASS" package. The package was tested on R/4.0.3. The operating system is macOS Monterey (Version 12.6)

## Example

An example of applying X-ING is:

```R
library(MASS)
library(CCA)
library(RGCCA)
data(example)
z_list <- lapply(example,function(x) x[[1]])
# L = 2
res <- XING(z_list = z_list[1:2])
# L = 3
res <- XING(z_list = z_list)
```

## Development

This package is maintained by Yihao Lu (yihaolu@uchicago.edu).

