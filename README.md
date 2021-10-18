
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GREAT <img src="man/figures/logo.png" align="right" width="120"/>

<!-- badges: start -->

[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Codecov test
coverage](https://codecov.io/gh/ruthkr/GREAT/branch/master/graph/badge.svg?token=L6TNLEPLLO)](https://codecov.io/gh/ruthkr/GREAT)
<!-- badges: end -->

The goal of GREAT (Gene Registration from Expression and Time-courses)
to register (align) gene expression profiles between two species
(reference data and data to transform). Non-reference gene expression
profiles will be stretched and shifted. The optimality of registration
parameters (shifts and stretches) will be estimated using least-squares
criterion. This package is also designed to compare a non-registration
model versus a non-registration model, as well as determine whether
registration model performed better than non-registration
transformation.

## Installation

<!--
You can install the released version of GREAT from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("GREAT")
```
-->

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ruthkr/GREAT")
```

## Quick start

This is a basic example which shows you how to register (align) gene
expression profiles over time:

``` r
# library(GREAT)
```
