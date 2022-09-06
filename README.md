
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/logo.png" align="right" style="float:right; height:140px;" />

# BHMSMAfMRI: Bayesian Hierarchical Multi-Subject Multiscale Analysis of fMRI Data

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/BHMSMAfMRI)](https://CRAN.R-project.org/package=BHMSMAfMRI)
<!-- badges: end -->

BHMSMAfMRI performs Bayesian hierarchical multi-subject multiscale
analysis of function MRI (fMRI) data as described in Sanyal & Ferreira
([2012](#ref-paper)), or other multiscale data, using wavelet based
prior that borrows strength across subjects and provides posterior
smooth estimates of the effect sizes and samples from their posterior
distribution.

## Installation

### Install from CRAN

``` r
install.packages("BHMSMAfMRI")
```

### Install from GitHub

``` r
# install.packages("devtools")
devtools::install_github("nilotpalsanyal/BHMSMAfMRI")
```

## The main function:

BHMSMA is the main function which accepts fMRI data as a 4D array (see
code below) and a design matrix. For the time-series of all voxels, a
general linear model (GLM) is fit with all the regressors in the design
matrix. After that, the standardized regression coefficient map of a
regressor of interest is subjected to further analysis. The function
BHMSMA returns the posterior smoothed image of the regression
coefficients. Below is a basic illustration of its use. For a detailed
manual, see the package <a
href="https://cran.r-project.org/web/packages/BHMSMAfMRI/vignettes/BHMSMAfMRIvignette.pdf"
target="_blank">vignette</a>.

``` r
library(BHMSMAfMRI)
#> Loading required package: Rcpp
n <- 3    #number of subjects
grid <- 8   #total number of voxels is grid^2
ntime <- 4  #number of timepoints
data <- array(rnorm(n*grid*grid*ntime), dim=c(n,grid,grid,ntime))
designmat <- cbind(c(1,1,1,1),c(1,0,1,0))
k <- 2  #regressor of interest
analysis <- "multi"     #perform multi-subject analysis
BHMSMAmulti <- BHMSMA(n, grid, data, designmat, k, analysis)
analysis <- "single"     #perform single subject analysis
BHMSMAsingle <- BHMSMA(n, grid, data, designmat, k, analysis)
```

## References:

<div id="refs" class="references">

<div id="ref-paper">

Sanyal, Nilotpal, and Ferreira, Marco A.R. (2012). Bayesian hierarchical
multi-subject multiscale analysis of functional MRI data. Neuroimage,
63, 3, 1519-1531. <span
target="_blank"><https://doi.org/10.1016/j.neuroimage.2012.08.041></span>
