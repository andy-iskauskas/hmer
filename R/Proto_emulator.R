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
    initialize = function(ranges, output_name, predict_func, variance_func, implausibility_func = NULL, print_func = NULL) {
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
    },
    get_exp = function(x) {
      return(self$predf(x))
    },
    get_cov = function(x) {
      return(self$varf(x))
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
      }
      else
        imp <- self$impf(x, z)
      if (!is.null(cutoff)) return(imp <= cutoff)
      return(imp)
    },
    print = function(...) {
      if (is.null(self$printf)) cat("Emulator protoype from custom object.\n")
      else return(self$printf(...))
    }
  )
)
