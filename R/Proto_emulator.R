Proto_emulator <- R6::R6Class(
  "EmProto",
  inherit = Emulator,
  public = list(
    ranges = NULL,
    output_name = NULL,
    predf = NULL,
    varf = NULL,
    impf = NULL,
    printf = NULL,
    active_vars = NULL,
    add_args = NULL,
    initialize = function(ranges, output_name, predict_func, variance_func,
                          implausibility_func = NULL, print_func = NULL, a_vars = NULL, ...) {
      self$ranges <- ranges
      self$add_args <- list(...)
      if (is.null(a_vars)) self$active_vars <- rep(TRUE, length(self$ranges))
      if (all(is.character(a_vars))) self$active_vars <- names(self$ranges) %in% a_vars
      if (length(self$active_vars) != length(self$ranges)) {
        if (length(self$active_vars) > length(self$ranges)) self$active_vars <- a_vars[1:length(self$ranges)]
        else self$active_vars <- c(a_vars, rep(FALSE, length(self$ranges) - length(a_vars)))
      }
      if (!all(is.logical(self$active_vars))) {
        tryCatch(as.logical(self$active_vars),
                 error = function(e) stop(paste("Could not parse active variables a_vars:", e)))
      }
      testPoint <- data.frame(matrix(purrr::map_dbl(ranges, mean), nrow = 1)) |> setNames(names(ranges))
      self$output_name <- output_name
      if (!is.function(predict_func)) stop("Prediction 'function' does not appear to be a function.")
      if (length(formals(predict_func)) < 1) stop("Prediction function requires at least one argument.")
      args <- formalArgs(predict_func)
      if (args[1] != "x") stop("First argument to 'predict_func' must be 'x'")
      if (length(args) > 1) {
        argsOK <- args[-1] %in% names(self$add_args)
        if (any(!argsOK)) stop(paste0("'predict_func' requires missing '",
                                      paste0(args[-1][!argsOK], collapse = ", "),
                                      "' object(s) to be passed to 'Proto_emulator'"))
        tryCatch(
          do.call(predict_func, c(list(x = testPoint), self$add_args[match(args[-1], names(self$add_args))])),
          error = function(e) {
            stop(paste("Predict function failing to evaluate at a point:", e))
          }
        )
      }
      else {
        tryCatch(
          predict_func(testPoint),
          error = function(e) {
            stop(paste("Predict function failing to evaluate at a point:", e))
          }
        )
      }
      self$predf <- predict_func
      if (!is.function(variance_func)) stop("Variance 'function' does not appear to be a function.")
      if (length(formals(variance_func)) < 1) stop("Variance function requires at least one argument.")
      args <- formalArgs(variance_func)
      if (args[1] != 'x') stop("First argument to 'variance_func' must be 'x'")
      if (length(args) > 1) {
        argsOK <- args[-1] %in% names(self$add_args)
        if (any(!argsOK)) stop(paste0("'variance_func' requires missing '",
                                      paste0(args[-1][!argsOK], collapse = ", "),
                                      "' object(s) to be passed to 'Proto_emulator'"))
        tryCatch(
          do.call(variance_func, c(list(x = testPoint),
                                   self$add_args[match(args[-1], names(self$add_args))])),
          error = function(e) {
            stop(paste("Variance function failing to evaluate at a point:", e))
          }
        )
      } else {
        tryCatch(
          variance_func(testPoint),
          error = function(e) {
            stop(paste("Variance function failing to evaluate at a point:", e))
          }
        )
      }
      self$varf <- variance_func
      if (!is.null(implausibility_func)) {
        if (!is.function(implausibility_func)) stop("Implausibility 'function' does not appear to be a function.")
        if (length(formals(implausibility_func)) < 2) stop("Implausibility function must have at least two arguments - input point(s) and target.")
        args <- formalArgs(implausibility_func)
        if (!identical(args[1:3], c('x', 'z', 'cutoff')))
          stop("First three arguments to 'implausibility_func' must be 'x', 'z', and 'cutoff")
        if (length(args) > 3) {
          argsOK <- args[-c(1:3)] %in% names(self$add_args)
          if (any(!argsOK)) stop(paste0("'implasuibilty_func' requires missing '",
                                        paste0(args[-c(1:3)][!argsOK], collapse = ", "),
                                        "' object(s) to be passed to 'Proto_emulator'"))
        }
      }
      self$impf <- implausibility_func
      if (!is.null(print_func)) {
        if (!is.function(print_func))
          stop("Print 'function' does not appear to be a function.")
        args <- formalArgs(print_func)
        if (!is.null(args)) {
          argsOK <- args %in% names(self$add_args)
          if (any(!argsOK)) stop(paste0("'print_func' requires missing '",
                                        paste0(args[!argsOK], collapse = ", "),
                                        "' object(s) to be passed to 'Proto_emulator'"))
        }
      }
      self$printf <- print_func
    },
    get_exp = function(x) {
      args <- formalArgs(self$predf)
      if (length(args) > 1) {
        return(do.call(self$predf,
                       c(list(x = x), self$add_args[match(args[-1], names(self$add_args))])))
      } else {
        return(self$predf(x))
      }
    },
    get_cov = function(x) {
      args <- formalArgs(self$varf)
      if (length(args) > 1) {
        return(do.call(self$varf,
                       c(list(x = x), self$add_args[match(args[-1], names(self$add_args))])))
      } else {
        return(self$varf(x))
      }
    },
    implausibility = function(x, z, cutoff = NULL) {
      if (is.null(self$impf)) {
        if (!is.numeric(z) && !is.null(z$val)) {
          if (is.null(z$sigma)) z$sigma <- 0
          imp_var <- self$get_cov(x) + z$sigma^2
          imp <- sqrt((z$val - self$get_exp(x))^2/imp_var)
        }
        else {
          pred <- self$get_exp(x)
          bound_check <- purrr::map_dbl(pred, function(y) {
            if (y <= z[2] && y >= z[1]) return(0)
            if (y < z[1]) return(-1)
            return(1)
          })
          which_compare <- purrr::map_dbl(bound_check, function(y) {
            if (y < 1) return(z[1])
            return(z[2])
          })
          uncerts <- self$get_cov(x)
          uncerts[uncerts <= 0] <- 1e-6
          imp <- bound_check * (pred - which_compare)/sqrt(uncerts)
        }
      } else {
        args <- formalArgs(self$impf)
        if (length(args) > 3) {
          return(do.call(self$impf,
                         c(list(x = x, z = z, cutoff = cutoff),
                           self$add_args[match(args[-c(1:3)], names(self$add_args))])))
        } else {
          return(self$impf(x, z, cutoff))
        }
      }
      if (!is.null(cutoff)) return(imp <= cutoff)
      return(imp)
    },
    print = function(...) {
      if (is.null(self$printf)) cat("Emulator prototype from custom object.\n")
      else {
        args <- formalArgs(self$printf)
        if(!is.null(args)) {
          return(do.call(self$printf, self$add_args[match(args, names(self$add_args))]))
        } else {
          return(self$printf)
        }
      }
    }
  )
)
