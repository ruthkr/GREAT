
<!-- README.md is generated from README.Rmd. Please edit that file -->

# greatR <img src="man/figures/logo.png" align="right" width="120"/>

<!-- badges: start -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/greatR)](https://cran.r-project.org/package=greatR)
[![lifecycle](https://lifecycle.r-lib.org/articles/figures/lifecycle-stable.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![R-CMD-check](https://github.com/ruthkr/greatR/workflows/R-CMD-check/badge.svg)](https://github.com/ruthkr/greatR/actions)
[![codecov](https://codecov.io/gh/ruthkr/greatR/branch/master/graph/badge.svg?token=L6TNLEPLLO)](https://app.codecov.io/gh/ruthkr/greatR)
[![pkgdown](https://github.com/ruthkr/greatR/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/ruthkr/greatR/actions/workflows/pkgdown.yaml)
[![GitHub last
commit](https://img.shields.io/github/last-commit/ruthkr/greatR)](https://github.com/ruthkr/greatR/commits/master)
<!-- badges: end -->

The goal of `greatR` (**G**ene **R**egistration from **E**xpression
**a**nd **T**ime-courses in **R**) to register (align) gene expression
profiles between two species (reference data and data to transform).
Non-reference gene expression profiles will be stretched and shifted.
The optimality of registration parameters (shifts and stretches) will be
estimated using least-squares criterion. This package is also designed
to compare a registration model versus a non-registration model, as well
as determine whether registration model performed better than
non-registration transformation.

## Package workflow

Flowchart below illustrates the workflow of the package given an input
data:

<img src="man/figures/greatR_workflow.png" width="85%" />

More details on how to use this package are available on function
documentations and vignettes:

1.  [Input
    requirements](https://ruthkr.github.io/greatR/articles/prepare-data.html)
2.  [Register
    data](https://ruthkr.github.io/greatR/articles/register-data.html)
3.  [Visualise registration
    results](https://ruthkr.github.io/greatR/articles/visualise-results.html)

## Installation

You can install the released version of `greatR` from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("greatR")
```

And the development version of `greatR` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ruthkr/greatR")
```

## Usage - quick start

This is a basic example which shows you how to register (align) gene
expression profiles over time:

``` r
# Load the package
library(greatR)
```

``` r
# Load a data frame from the sample data
b_rapa_data <- system.file("extdata/brapa_arabidopsis_all_replicates.csv", package = "greatR") |>
  utils::read.csv()

# Running the registration
registration_results <- register(
  b_rapa_data,
  reference = "Ro18",
  query = "Col0"
)
#> ℹ Will process 10 genes.
#>
#> ── Starting registration with optimisation ──────────────────────────────────────
#> ℹ Using computed stretches and shifts search space limits.
#> ✔ Optimising registration parameters for genes (10/10) [5m 14.8s]
```

## Reference

Calderwood, A., Hepworth, J., Woodhouse, … Morris, R. (2021).
Comparative transcriptomics reveals desynchronisation of gene expression
during the floral transition between Arabidopsis and Brassica rapa
cultivars. *Quantitative Plant Biology, 2*, E4.
[doi:10.1017/qpb.2021.6](https://www.cambridge.org/core/journals/quantitative-plant-biology/article/comparative-transcriptomics-reveals-desynchronisation-of-gene-expression-during-the-floral-transition-between-arabidopsis-and-brassica-rapa-cultivars/811BFDFA14F4BCC9C7F0ECC7CE103BB6)
