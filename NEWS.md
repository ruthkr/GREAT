# greatR 1.0.0

* Rewrote registration pipeline from scratch, deprecating unnecessary, and redundant auxiliary functions.
* Added L-BFGS-B and Nelder-Mead (now default) optimisation methods to {greatR}.
* Switched to manual calculation of log likelihood via `calc_loglik()` instead of `stats::logLik()`.
* Reduced computation time up to 1000 times, (x30 speed-up from package rewrite, and x35 speed-up from switching default optimisation method).
* Removed {dplyr}, {magrittr}, {purrr}, {rlang}, and {stringr} as package dependencies.
* Added {neldermead} as a package depedency.
* Updated list of exported functions:
  * `register()`
  * `summarise_registration()`
  * `get_approximate_stretch()`
  * `plot_registration_results()`
  * `plot_heatmap()`
  * `calculate_distance()`

## Improvements

* Simplified parameters of main `register()` function, and added `scaling_method`.
* Simplified structure of output object of `register()`.
* Simplfied parameters of `summarise_registration()`, `plot_registration_results()`, `plot_heatmap()`, `calculate_distance()` to simply require `results` object from `register()`, vastly simplifing usage.
* Improved messages, errors, and progress indicators with {cli}.
* Added correct pluralisation in {cli} messages.
* Rewrote unit tests to use {data.table} exclusively for data manipulation.
* Added unit tests for `calc_loglik_H1()`, `calc_loglik_H2()`, `calc_overlapping_percent()`, `calculate_distance()`, `cross_join()`, `get_search_space_limits_from_params()`, `get_search_space_limits()`, `objective_fun()`, `optimise()`, `plot_heatmap()`, `plot_registration_results()`, `preprocess_data()`, `register_manually()`, `register()`, `summary_registration()`, `validate_params()`.

## Bug fixes

* Fixed `match_names()` call when validating accession names in `register()`
* Fixed use of deprecated `aes_string()` by parsing `timepoint_var` using `!!ggplot2::sym()` call.
* Fixed `preds` left join in `plot_registration_results()`.
* Fixed issue in `plot_registration_results()` not working when all genes are unregistered with `type = "registered"`.
* Fixed calculation of `time_delta` in `preprocess_data()` to ensure it's grouped by `gene_id` and `accession` (not just `accession`).

# greatR 0.2.0

* Added Alex Calderwood as package co-author.
* Added vignette for optimisation process.
* Refactored `num_shifts` and `shift_extreme` parameters by simplified `shifts` parameter.

## Improvements

* Improved default parameter values in exported functions.
* Added {optimization}, {purrr} as package dependencies.
* Removed {cowplot}, {ggpubr}, {ggrepel}, {Rtsne}, and {viridis} as package dependencies.
* Cleaned up {cli} messages.
* Removed legacy AIC references, as it is no longer used.
* Updated `calculate_between_sample_distance()` to use `registration_results` as primary parameter instead of `mean_df`, `mean_df_sc`, and `imputed_mean_df`.
* Added warning if there is no comparable time points found using users' pre-defined parameters.
* Refactored `optimise_shift_extreme` as `maintain_min_num_overlapping_points`, properly defined and corrected the boundary box if number overlapping points whether needed to be maintained or not.

## Bug fixes

* Check that input accessions exist in the input data in `get_approximate_stretch()`.
* Manually create time point sorting levels for `x_sample` and `y_sample` columns according in `plot_heatmap()`.
* Properly handle `-` character in accession names in `plot_heatmap()` so that time points are parsed correctly.

## New features

* Added optional parameter optimisation process using Simulated Annealing through `optimise_registration_params()`.

## New functions

* `preprocess_data()` to simplify `scale_and_register_data()` code and reuse logic elsewhere.
* `get_best_stretch_and_shift_simplified()`.
* `get_BIC_from_registering_data()`.
* `get_boundary_box()`.
* `optimise_registration_params_single_gene()`.
* `optimise_registration_params()` as wrapper of `optimise_registration_params_single_gene()` for multiple genes.
* `get_best_stretch_and_shift_after_optimisation()`.

# greatR 0.1.0

* Initial release.
* Added a `NEWS.md` file to track changes to the package.
