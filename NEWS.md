# hmer 1.5.5

## Major changes

* New functionality: `diagnostic_pass` for automated diagnostics and modifications of
emulators; `hit_by_wave` for visualisation of history matching progress.

## Bug fixes

* Fixed issue in `simulator_plot` et al. where some combinations of `wave_numbers`
and `zero_in` caused out-of-index errors.
* Edge case fixes for training individual variance emulators or one-dimensional
variance emulators.

## Enhancements

* Vastly improved efficiency of emulator training when `emulator_type = "variance"` or
`emulator_type = "covariance"`.
* Improvement to hyperparameter estimation when training covariance emulators.

# hmer 1.5.0

## Major changes

* Change to nomenclature: all `xx_emulator_from_data` functions are now a single function; `variance_emulator_from_data` and `bimodal_emulator_from_data` are now called using `emulator_from_data` with argument `emulator_type = 'variance'` and `emulator_type = 'multistate'` respectively. To avoid confusion about the output of the function, `generate_new_runs` has been renamed to `generate_new_design`. Older functions have been deprecated and will be removed in a subsequent version.

## Bug fixes

* Changes to `generate_new_design` to avoid implausibility asymptoting in edge cases.

* Explicit checks included to ensure that ordering of parameter ranges matches with data.frames provided for training

* Fixes to ensure that 1d systems behave as expected.

* Fix to avoid singular correlation matrices in some situations where `corr_type = 'exp_sq'`.

* Fixes for edge-case multistate emulator training problems, and for situations where `model.matrix` behaviour fails (due to `deparse` truncation issues in core R functions).

## Enhancements

* Modifications to emulator design via `emulator_from_data`: in particular hyperparameter estimation has been made more robust; emulator regression surfaces now support cubic terms; variance and multistate emulators now more robust to different numbers of repetitions at different input sites.

* Covariance emulation introduced (via `emulator_from_data(..., emulator_type = 'covariance')`), allowing a full covariance matrix to be robustly emulated. Results are presented using a prototype emulator matrix `EmulatedMatrix` to efficiently make predictions, including checks to ensure that predicted matrices are semi-positive definite. In the future, covariance matrices will be supported in implausibility measures and point proposals.

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
