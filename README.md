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
```

The expected output includes the posterior mean estimation and posterior probability. The computation time was listed.

$L=2$:

```R
# L = 2
> start.time <- proc.time()
> res <- XING(z_list = z_list[1:2])
[1] "Iteration times: 2"
[1] "Iteration times: 3"
[1] "Iteration times: 4"
[1] "Iteration times: 5"
[1] "Iteration times: 6"
> proc.time() - start.time
   user  system elapsed 
  1.282   0.394   1.700 
> names(res)
[1] "mu"    "alpha"
```

$L=3$:

```R
# L = 3
> start.time <- proc.time()
> res <- XING(z_list = z_list)
[1] "Iteration times: 2"
[1] "Iteration times: 3"
[1] "Iteration times: 4"
[1] "Iteration times: 5"
> proc.time() - start.time
   user  system elapsed 
 71.025  23.524  96.468 
> names(res)
[1] "mu"    "alpha"
```

## Development

This package is maintained by Yihao Lu (yihaolu@uchicago.edu).

## Licenses

The code is released under [GNU General Public License (GPL)](https://opensource.org/licenses/gpl-license).
