# hmer 1.4.0

## Bug fixes

* Multiple fixes to deal with singletons: ensuring 1d examples work with `generate_new_runs`; single output systems behave appropriately under clustering, diagnostics, and implausibility for `bimodal_emulator_from_data`; modifications to `simulator_plot` for single output systems; bug-fix in `generate_new_runs` for user-provided single-element plausible sets.

* Small fixes to guarantee compatibility with new functionality.

* Fixed scoping issue with variance emulators where emulator would not get the correct `s_diag` function when initialised within a collection of emulators in `variance_emulator_from_data`

## Enhancements

* Beeswarm plots added as alternative to `simulator_plot` (credit to T.J. McKinley)

* Custom multi-emulator implausibility now supported within `generate_new_runs`: multiple conditions can be supplied as part of the point-screening process via `accept_measure`. Structure of optional arguments for `generate_new_runs` has been modified (with backwards-compatibility for older code) - see help file for details.

* Added check to `emulator_from_data` to handle mismatched input names and ranges

## New Functionality

* `Proto_emulator` introduced: the `hmer` framework of diagnostics, visualisation, and point proposal can be used with entirely custom objects.

# hmer 1.2.0

## Bug fixes

* Fix to implausibility for variance emulators to take account of ensemble size

* Functions `full_wave`, `variance_emulator_from_data`, `bimodal_emulator_from_data` all accept an `na.rm` argument to handle missing data

* Fixed edge cases where `generate_new_runs` could get stuck at a particular implausibility cutoff, and increased stability of termination for points generation from variance/bimodal emulators

* Other small fixes, including modification to `ggplot` functions to handle deprecation of `size` aesthetic

## Enhancements

* Optimisation of emulator calculations, particularly within correlation matrices

* Optimised point generation to leverage Latin Hypercube Designs, where useful

* Modification to `standard_errors` to highlight points of interest

* Modifications made to faciliate support for custom emulators and implausibility measures (to come in a later update).

# hmer 1.0.1

## Bug Fixes

* Implemented more error catching/handling for correlation functions, including explicit stop calls if hyperparameters are not provided/ill-specified and automatic coercion to data.matrix to ensure compatibility with derivative functions

* Internal functions `collect_emulators`, `scale_input`, `convertRanges` and `multiply_function` modified to handle various edge case usages and issues with multiple waves of stochastic or bimodal emulators

* Modified `lhs_gen_cluster` so that if emulators cannot/need not be clustered, default `lhs_gen` behaviour is used

* Fixed `Emulator` code to address calculation issues with derivatives that meant that partial derivatives did not commute

* Various other small fixes

## Enhancements

* Unit testing framework implemented and tracked; unit tests for all relevant functions are now in place

* Functions `generate_new_runs` and `nth_implausible` will try to determine a sensible value of `n` or `nth` when calculating nth-maximum implausibility if none is supplied by the user.

* Badges and sticker added to readme

# hmer 1.0.0

* Added a `NEWS.md` file to track changes to the package.
