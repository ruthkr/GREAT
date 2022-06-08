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
