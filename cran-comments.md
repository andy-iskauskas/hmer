## Resubmission

This is a resubmission due to failed test cases in automated CRAN checking.
Automated checking failed on single `testhat` block (`test-emulator.R`, lines 245-257). Problem traced to poorly-conceived randomisation in test that caused failure in ~1% of test runs. Unit test has been modified to accurately reflect the required test case. No other package content required modification.

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

## Downstream Dependencies
There are no downstream dependencies.

## Test Environments
* Local MacOS install R 4.3.1
* Local Windows 10 install R 4.3.1
* `devtools::check_win_release()`
* `devtools::check_win_devel()`
* `use_github_actions()`

## Original Submission

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

## Downstream Dependencies
There are no downstream dependencies.

Test Environments
* Local MacOS install R 4.3.1
* Local Windows 10 install R 4.3.1
* `devtools::check_win_release()`
* `devtools::check_win_devel()`
* `use_github_actions()`
* `rhub::check_for_cran()`
