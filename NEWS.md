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
