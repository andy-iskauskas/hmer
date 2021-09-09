VarianceEmulator <- R6::R6Class(
  "VarianceEmulator",
  inherit = Emulator,
  public = list(
    v_em = NULL,
    kurt = NULL,
    initialize = function(v_em, kurt, ...) {
      super$initialize(...)
      self$v_em <- v_em
      self$kurt <- kurt
    },
    get_exp = function(x, n_vals = NULL) {
      if (!is.null(n_vals)) {
        if (length(n_vals) == 1) n_vals <- rep(n_vals, nrow(x))
        if (length(n_vals) != nrow(x)) stop("Number of sample sizes does not match number of points.")
        return(super$get_exp(x) + diag(self$v_em$get_exp(x)/n_vals))
      }
      return(super$get_exp(x))
    },
    get_cov = function(x, full = FALSE, n_vals = NULL) {
      if (!is.null(n_vals)) {
        if (length(n_vals) == 1) n_vals <- rep(n_vals, nrow(x))
        if (length(n_vals) != nrow(x)) stop("Number of sample sizes does not match number of points.")
        vT <- (self$v_em$get_exp(x)^2 + self$v_em$get_cov(x))/n_vals * (self$kurt-1+2)/(n_vals-1)
        if (full) return(super$get_cov(x, full = TRUE) + diag(vT))
        else return(super$get_cov(x) + vT)
      }
      return(super$get_cov(x, full = full))
    }
  )
)
