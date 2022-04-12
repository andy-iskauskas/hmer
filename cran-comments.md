## Resubmission

This is a resubmission to address comments from the initial submission from 07/04/2022. In this submission:

* I have removed `print` statements and replaced them with a `warning`, or a `cat` statement called only if `verbose` is `TRUE`. Functions that did not have a verbose argument have had one added (applied in functions `summary_diag`, `Emulator$implausibility`, `preflight`, `full_wave`, `generate_new_runs` and children `<xxx>_sample`, `idemc`, `<xxx>_emulator_from_data`, `eval_funcs`).

* In `effect_strength`, plotting is now optional and enabled via a `plt` argument. In a similar vein, calling `diagnostic_wrap` returns a list of ggplot objects which can be stored or printed as the user deems appropriate, rather than automatically printing to console.

* Where a function calls `par(op)`, an `on.exit` statement has been added to the start of the function to ensure no user options are irrevocably changed on failure (applied to functions `behaviour_plot`, `validation_diagnostics`).

* Examples that were wrapped in `dontrun` have been changed to `donttest` where the computation time would exceed 5s (functions `idemc`, `bimodal_emulator_from_data`).

* I have clarified with CRAN checkers that functions that can save to file (`generate_new_runs`, `importance_sample`, `diagnostic_wrap`) do not violate CRAN policy on the matter since they require the user to explicitly provide a location (the default argument being `NULL`, in which case no attempt to save is made). No examples or vignettes use the optional argument.

* I have removed `wave_variance`; upon inspection it is clear that the result of this function can be easily replicated in other functions. I did not deem deprecation necessary due to the 'new submission' nature of the package.

## R CMD check results

0 errors | 0 warnings | 1 note

There was 1 NOTE:

* This is a new release.

## Downstream dependencies

There are no downstream dependencies.

## Test environments

Tested with

* R 4.1.3 (Windows) `devtools::check_win_release()`
* R 4.2.0 alpha (Windows) `devtools::check_win_devel()`
* `rhub::check_for_cran()`
* `use_github_actions()`
