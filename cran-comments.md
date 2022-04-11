## Resubmission
This is a resubmission to address comments from initial submission. In this submission I have:

* Removed `print` statements and replaced them with a `warning`, or a `cat` statement called only if `verbose` is `TRUE`. Functions that did not have a verbose argument have had one added (applied in functions `summary_diag`, `Emulator$implausibility`, `preflight`, `full_wave`, `generate_new_runs` and children `<xxx>_sample`, `idemc`, `<xxx>_emulator_from_data`).

* In `effect_strength`, plotting is now optional and enabled via a `plt` argument.

* Where a function calls `par(op)`, an `on.exit` statement has been added to the start of the function to ensure no user options are irrevocably changed (applied to functions `behaviour_plot`, `validation_diagnostics`).

* Examples that were wrapped in `dontrun` have been changed to `donttest` where the computation time would exceed 5s (functions `idemc`, `bimodal_emulator_from_data`).

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
