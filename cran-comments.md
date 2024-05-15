## Resubmission
This is a re-resubmission following manual check:
* Added DOIs for representative references to history matching and emulation in DESCRIPTION
* Tarball should now be <5Mb; local compilation had caused unnecessary bloat in figures (checked using `check_win_release()` and `check_win_devel()`).

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

## Downstream Dependencies
There are no downstream dependencies.

## Test Environments
* Local MacOS install R 4.4.0
* Local Windows 10 install R 4.4.0
* `devtools::check_win_release()`
* `devtools::check_win_devel()`
* `use_github_actions()`
* `rhub::rhub_check()`
