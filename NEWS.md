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
