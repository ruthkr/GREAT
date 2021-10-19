
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GREAT <img src="man/figures/logo.png" align="right" width="120"/>

<!-- badges: start -->

[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Codecov test
coverage](https://codecov.io/gh/ruthkr/GREAT/branch/master/graph/badge.svg?token=L6TNLEPLLO)](https://codecov.io/gh/ruthkr/GREAT)
<!-- [![](https://img.shields.io/github/last-commit/federicomarini/GeneTonic.svg)](https://github.com/federicomarini/GeneTonic/commits/master) -->
<!-- badges: end -->

## Overview

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

And the development version of `GREAT` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ruthkr/GREAT")
```

## Usage - quick start

This is a basic example which shows you how to register (align) gene
expression profiles over time:

``` r
# Load the package
devtools::load_all()
library(dplyr)
```

``` r
# Define the dataframe from the sample data
# Gene expression data with replicates
all_data_df <- system.file("extdata/brapa_arabidopsis_all_replicates.csv", package = "GREAT") %>%
  utils::read.csv()

# Mean gene expression 
mean_df <- system.file("extdata/brapa_arabidopsis_mean.csv", package = "GREAT") %>%
  utils::read.csv()

# Running the registration
registration_results <- scale_and_register_data(
    mean_df = mean_df,
    all_data_df = all_data_df,
    stretches = c(2, 1.5, 1),
    shift_extreme = 4,
    num_shifts = 27,
    min_num_overlapping_points = 4,
    initial_rescale = FALSE,
    do_rescale = TRUE,
    accession_data_to_transform = "Col0",
    accession_data_ref = "Ro18",
    data_to_transform_time_added = 11,
    data_ref_time_added = 11
  )
```

``` r
# Plot registration result
plot_registered_gene_of_interest(registration_results[["imputed_mean_df"]], ncol = 3)
```

<img src="man/figures/README-plot-results-1.png" width="100%" />

More examples are available on functions documentations and vignettes,
please refer to the documentation.

## Reference

[Calderwood, A., Hepworth, J., Woodhouse, . . . Morris, R. (2021).
Comparative transcriptomics reveals desynchronisation of gene expression
during the floral transition between Arabidopsis and Brassica rapa
cultivars. Quantitative Plant Biology, 2, E4.
doi:10.1017/qpb.2021.6](https://www.cambridge.org/core/journals/quantitative-plant-biology/article/comparative-transcriptomics-reveals-desynchronisation-of-gene-expression-during-the-floral-transition-between-arabidopsis-and-brassica-rapa-cultivars/811BFDFA14F4BCC9C7F0ECC7CE103BB6)
