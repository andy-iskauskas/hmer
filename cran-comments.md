## Resubmission
This is a resubmission, due to automatic test failure: long build time on r-devel-windows-x86_64 due to vignette computational time. In this version I have:

* Pre-compiled all vignettes to ensure acceptable checktime

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
