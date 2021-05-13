# hmer
Bayes Linear Emulation Package (R)

## A package for Bayes Linear Emulation and History Matching ##

### Background ###
In many fields, it is necessary to build a complex model to represent real-life systems. Such simulations can encompass a large, or high-dimensional, input space, and the simulator can be slow to evaluate over points in the space. Examples include [simulation of galaxy formation](http://dro.dur.ac.uk/15070/), [HIV transmission models](http://dro.dur.ac.uk/22952/), or [hormonal crosstalk in biological systems](https://arxiv.org/abs/1801.01538). In these, and other, situations, it is useful to be able to generate large portions of the simulator space without explicit recourse to the simulator.

An emulator is a statistical model of the simulator, using a relatively small number of sample points taken from simulator runs. The broad, simplistic form of the emulator is <img src="/tex/4aadb878bfd05154ff967d5f5d539b9a.svg?invert_in_darkmode&sanitize=true" align=middle width=146.3736615pt height=24.65753399999998pt/>, where the <img src="/tex/ffcbbb391bc04da2d07f7aef493d3e2a.svg?invert_in_darkmode&sanitize=true" align=middle width=30.61077854999999pt height=24.65753399999998pt/> and <img src="/tex/8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode&sanitize=true" align=middle width=10.16555099999999pt height=22.831056599999986pt/> create a regression surface, and the <img src="/tex/320b3450fd8b780975b68c70115439b3.svg?invert_in_darkmode&sanitize=true" align=middle width=31.590693749999986pt height=24.65753399999998pt/> encodes correlation between points in the input space. In this way, the emulator represents the simulator results, but also has a built-in understanding of the uncertainties of both the emulation and of the simulator uncertainties. The value of the emulator at a point x is given by <img src="/tex/952617a3912bcad4e130e08f0389caf0.svg?invert_in_darkmode&sanitize=true" align=middle width=52.31752514999999pt height=24.65753399999998pt/>, with associated variance <img src="/tex/0986ffe08d5d83601df0378fabb990fa.svg?invert_in_darkmode&sanitize=true" align=middle width=66.74676194999999pt height=24.65753399999998pt/>.

Given such an emulator structure, outputs from the simulator can be used to 'train' the emulator using Bayesian methods. In this package, we forego a full Bayesian approach and instead focus on Bayes linear updates (for details, see e.g. [this book](https://onlinelibrary.wiley.com/doi/book/10.1002/9780470065662)). This has the advantage of being quick to evaluate, and does not require us to supply full probabilistic specifications for all parameters of the emulator.

Often, we don't care about the structure of the emulator (or the simulator), but instead the following question: **what input values would give us agreement with observations?** To answer this question, we can use History Matching. We define an *implausibility measure* I(x), which in some sense indicates how likely a set of input values are to give a particular output. By applying a suitable cut-off to this measure, we can generate a non-implausible space for a desired output, from which we can sample more points, from which we can train a better emulator. The process can be repeated until either:
- We have enough non-implausible points to make an informed statement about the structure of the model;
- The emulator uncertainty is smaller than the uncertainties inherent in the system (e.g. observation error, model discrepancy)

The method of Bayes Linear emulation and history matching is a powerful tool to efficiently and robustly explore large, complex, time-intensive models, and find informative regions of input space that could potentially be overlooked by recourse only to the simulator.

### Package Details ###
**Note: the package is a work in progress! Use with care.**

The package is composed of an object definition (as an R6 class), and functions that use such objects. For details of usage, see the associated help files and examples therein. The following is a summary.

 The R6 class is
- ``Emulator``: This is the key object in the package. Given the second-order specifications for <img src="/tex/8217ed3c32a785f0b5aad4055f432ad8.svg?invert_in_darkmode&sanitize=true" align=middle width=10.16555099999999pt height=22.831056599999986pt/> and <img src="/tex/320b3450fd8b780975b68c70115439b3.svg?invert_in_darkmode&sanitize=true" align=middle width=31.590693749999986pt height=24.65753399999998pt/> (namely their expectations, variances, and covariances), and the basis functions <img src="/tex/ffcbbb391bc04da2d07f7aef493d3e2a.svg?invert_in_darkmode&sanitize=true" align=middle width=30.61077854999999pt height=24.65753399999998pt/>, an emulator is generated. From this we can find expectation and variance of the emulator at input points; we can use data points from a simulator to train new emulators based on these prior specifications; and we can determine implausibility for a given input against a desired output. The specifications can be definde by hand, up to and including nugget terms; alternatively, emulators for known data can be derived as a 'first pass' using `emulator_from_data`.

Associated functions (not an exhaustive list; see ``lsf.str("package:emulatorr")`` and associated help files):
- ``emulator_from_data``: Generates basic univariate emulators from simulator data, along with some simple prior information
- Diagnostic functions (``standard_error``, ``comparison_diagnostics``, ...): Shows diagnostic tests on an emulator given a validation set of data
- ``nth_implausible``: For a collection of outputs, and univariate emulators trained to each, finds the <img src="/tex/55a049b8f161ae7cfeb0197d75aff967.svg?invert_in_darkmode&sanitize=true" align=middle width=9.86687624999999pt height=14.15524440000002pt/>-th maximum implausibility for an input point against a desired set of outputs.
- ``generate_new_runs``: Given a set of trained emulators, suggests a new set of points to be entered into the simulator.
- Plotting functions (``output_plot``, ``emulator_plot``, ...): Various plotting functions for visualising output.
- ``full_wave``: Rolls together the key components of a wave, generating emulators from data, checking diagnostics, finding implausibilities, and creating a new set of points on which to run the simulator.

### Installation ###
Using RStudio with ``devtools`` enabled, the command ``devtools::install_github("Tandethsquire/emulatorr")`` will install the package. Henceforth, entering ``library(emulatorr)`` will make it available for use. For generic help on the package, use ``help(package = 'emulatorr')`` and for a specific function use the ``?`` syntax: e.g. ``?generate_new_runs``.
An introductory vignette is available in the ``vignettes`` folder: compiled versions are available in the ``doc`` folder.

### Comments and Suggestions ###
Suggestions for improvement, and comments on functionality, are welcome. For preference these should be addressed via pull requests and/or issue logging, but if in doubt email andrew.iskauskas@durham.ac.uk.
