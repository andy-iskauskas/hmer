# --------------------
# Emulator Prototype Documentation
# --------------------
#' @title Prototype class for emulator-like objects
#'
#' @description Converts a prediction object into a form amenable to hmer.
#'
#'     The history matching process can be used for objects that are not
#'     created by the \code{hmer} package: most notably for Gaussian Process
#'     Emulators, but even for simple linear models. This R6 class converts
#'     an object into a form that can be called reliably by the methods of the
#'     package, including for visualisation and diagnostics.
#'
#' @name Proto_emulator
#'
#' @importFrom R6 R6Class
#'
#' @section Constructor: \code{Proto_emulator$new(ranges, output_name, predict_func, variance_func, ...)}
#'
#' @section Arguments:
#'
#'    Required:
#'
#'    \code{ranges} A list of ranges for the inputs to the model.
#'    \code{output_name} The name of the output modelled.
#'
#'
#'    \code{predict_func} The function that provides the predictions at a new
#'    point. This should be a function with first argument, \code{x}, which
#'    can be a single point or a data.frame of points.
#'
#'    \code{variance_func} The function that encodes the prediction error due to
#'    the model of choice. This, too, takes a first argument \code{x} of the same form
#'    as that of \code{predict_func}.
#'
#'    Optional:
#'
#'    \code{implausibility_func} A function that takes points \code{x} and a
#'    target \code{z} (and possibly a cutoff value \code{cutoff} and additional arguments) 
#'    and returns a measure of closeness of the predicted value to the target (or a boolean
#'    representing whether the prediction is within the specified amount).
#'    If not provided, then the standard implausibility is used: namely the
#'    absolute value of the prediction minus the observation, dividied by the
#'    square root of the sum in quadrature of the errors.
#'
#'    \code{print_func} If the prediction object has a suitable print function
#'    that one wishes to transfer to the R6 class, it is specified here.
#'
#'    \code{...} Additional objects to pass to emulators and/or implausibility
#'    measures.
#'
#' @section Constructor Details:
#'
#'    The constructor must take, as a minimum, the first four arguments (the
#'    input ranges, output name, and functions to provide prediction and
#'    prediction error). Default behaviour exists if an implausibility function
#'    is not specified. The output of the constructor is an R6 object with the
#'    classes "Emulator" and "EmProto".
#'
#' @section Accessor Methods:
#'
#'    Note that these have the same external structure as those in
#'    \code{\link{Emulator}}.
#'
#'    \code{get_exp(x)} Returns the prediction.
#'
#'    \code{get_cov(x)} Returns the prediction error.
#'
#'    \code{implausibility(x, z, cutoff = NULL)} Returns the 'implausibility'.
#'
#'    \code{print()} Prints relevant details, where appropriate.
#'
#' @export
#'
#' @examples
#'    # Use linear regression with an "error" on the SIR dataset
#'    ranges <- list(aSI = c(0.1, 0.8), aIR = c(0, 0.5), aSR = c(0, 0.05))
#'    targets <- SIREmulators$targets
#'    lms <- purrr::map(names(targets),
#'        ~step(lm(data = SIRSample$training[,c(names(ranges), .)],
#'        formula = as.formula(paste0(., "~(",
#'            paste0(names(ranges), collapse = "+"),
#'            ")^2"
#'        ))
#'    ), trace = 0))
#'    # Set up the proto emulators
#'    # Note that we name the list after creating it with the output names
#'    proto_ems <- purrr::map(seq_along(lms), function(l) {
#'        Proto_emulator$new(
#'            ranges,
#'            names(targets)[l],
#'            function(x) predict(lms[[l]], x),
#'            function(x) predict(lms[[l]], x, se.fit = TRUE)$se.fit^2 +
#'               predict(lms[[l]], x, se.fit = TRUE)$residual.scale^2,
#'            print_func = function() print(summary(lms[[l]]))
#'        )
#'    }) |> setNames(names(targets))
#'    ## Test it with hmer functions
#'    nth_implausible(proto_ems, SIRSample$validation, targets)
#'    emulator_plot(proto_ems)
#'    emulator_plot(proto_ems, 'imp', targets = targets)
#'    validation_diagnostics(proto_ems, targets, SIRSample$validation)
#'    new_points <- generate_new_runs(proto_ems, 100, targets)
#'
NULL
