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
    add_args = NULL,
    initialize = function(ranges, output_name, predict_func, variance_func, implausibility_func = NULL, print_func = NULL, ...) {
      self$ranges <- ranges
      self$output_name <- output_name
      ## Want to add tests on the predict_func here
      self$predf <- predict_func
      ## Want to add tests on the variance_func here
      ## Also: question of whether we check for two arguments (eg for covariance matrix)
      self$varf <- variance_func
      ## Checks here (if not null)
      self$impf <- implausibility_func
      self$printf <- print_func
      self$add_args <- list(...)
    },
    get_exp = function(x) {
      args <- formalArgs(self$predf)
      if(args[1] != "x") stop("First argument to 'predf' must be 'x'")
      if(length(args) > 1) {
          argsOK <- args[-1] %in% names(self$add_args)
          if(any(!argsOK)) stop(paste0("You must pass '", paste0(args[-1], collapse = ", "), "' arguments for 'predf' when setting up emulator"))
          return(do.call(self$predf, c(list(x = x), self$add_args[match(args[-1], names(self$add_args))])))
      } else {
          return(self$predf(x))
      }
    },
    get_cov = function(x) {
      args <- formalArgs(self$varf)
      if(args[1] != "x") stop("First argument to 'varf' must be 'x'")
      if(length(args) > 1) {
          argsOK <- args[-1] %in% names(self$add_args)
          if(any(!argsOK)) stop(paste0("You must pass '", paste0(args[-1], collapse = ", "), "' arguments for 'varf' when setting up emulator"))
          return(do.call(self$varf, c(list(x = x), self$add_args[match(args[-1], names(self$add_args))])))
      } else {
          return(self$varf(x))
      }
    },
    implausibility = function(x, z, cutoff = NULL) {
      if (is.null(self$impf)) {
        if (!is.numeric(z) && !is.null(z$val)) {
          imp_var <- self$get_cov(x) + z$sigma^2
          imp <- sqrt((z$val - self$get_exp(x))^2/imp_var)
        }
        else {
          pred <- self$get_exp(x)
          bound_check <- purrr::map_dbl(pred, function(y) {
            if (y <= z[2] && y >= z[1]) return(0)
            if (y < z[1]) return(-1)
            if (y > z[2]) return(1)
          })
          which_compare <- purrr::map_dbl(bound_check, function(y) {
            if (y < 1) return(z[1])
            return(z[2])
          })
          uncerts <- self$get_cov(x)
          uncerts[uncerts <= 0] <- 0.0001
          imp <- bound_check * (pred - which_compare)/sqrt(uncerts)
        }
      } else {
        args <- formalArgs(self$impf)
        if(!identical(args[1:3], c("x", "z", "cutoff"))) stop("First three arguments to 'impf' must be 'x', 'z' and 'cutoff'")
        if(length(args) > 3) {
          argsOK <- args[-c(1:3)] %in% names(self$add_args)
          if(any(!argsOK)) stop(paste0("You must pass '", paste0(args[-c(1:3)], collapse = ", "), "' arguments for 'impf' when setting up emulator"))
          return(do.call(self$impf, c(list(x = x, z = z, cutoff = cutoff), self$add_args[match(args[-c(1:3)], names(self$add_args))])))
        } else {
          return(self$impf(x, z, cutoff))
        }
      }
    },
    print = function(...) {
      if (is.null(self$printf)) cat("Emulator protoype from custom object.\n")
      else {
        args <- formalArgs(self$printf)
        if(length(args) > 0) {
            argsOK <- args %in% names(self$add_args)
            if(any(!argsOK)) stop(paste0("You must pass '", paste0(args, collapse = ", "), "' arguments for 'printf' when setting up emulator"))
            return(do.call(self$printf, self$add_args[match(args, names(self$add_args))]))
        } else {
            return(self$printf())
        }
      }
    }
  )
)
