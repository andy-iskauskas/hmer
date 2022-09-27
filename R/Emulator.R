Emulator <- R6Class(
  "Emulator",
  private = list(
    data_corrs = NULL,
    design_matrix = NULL,
    u_exp_modifier = NULL,
    u_var_modifier = NULL,
    beta_u_cov_modifier = NULL
  ),
  public = list(
    basis_f = NULL,
    beta_mu = NULL,
    beta_sigma = NULL,
    u_mu = NULL,
    u_sigma = NULL,
    corr = NULL,
    beta_u_cov = NULL,
    in_data = NULL,
    out_data = NULL,
    ranges = NULL,
    model = NULL,
    model_terms = NULL,
    active_vars = NULL,
    o_em = NULL,
    output_name = NULL,
    disc = list(internal = 0, external = 0),
    multiplier = 1,
    initialize = function(basis_f, beta, u, ranges, data = NULL, model = NULL,
                          original_em = NULL, out_name = NULL, a_vars = NULL,
                          discs = NULL, multiplier = 1) {
      self$model <- model
      self$model_terms <- tryCatch(
        c("1", labels(terms(self$model))),
        error = function(e) return(NULL)
      )
      self$o_em <- original_em
      self$basis_f <- basis_f
      self$beta_mu <- beta$mu
      self$beta_sigma <- beta$sigma
      self$u_mu <- function(x) 0
      self$multiplier <- multiplier
      if(is.numeric(u$sigma)) self$u_sigma <- self$multiplier * u$sigma
      else self$u_sigma <- function(x) self$multiplier * (u$sigma)(x)
      self$corr <- u$corr
      if (!is.null(out_name)) self$output_name <- out_name
      if (is.null(ranges)) stop("Ranges for the parameters must be specified.")
      self$ranges <- ranges
      if (is.null(a_vars)) {
        self$active_vars <- purrr::map_lgl(seq_along(ranges), function(x) {
          point_vec <- c(rep(0, x-1), 1, rep(0, length(ranges)-x))
          func_vals <- purrr::map_dbl(self$basis_f, purrr::exec, point_vec)
          sum(func_vals) > 1
        })
      }
      else self$active_vars <- a_vars
      if (all(self$active_vars == FALSE))
        self$active_vars <- rep(TRUE, length(ranges))
      if (!is.null(discs)) {
        self$disc$internal <- ifelse(!is.null(discs$internal),
                                     discs$internal, 0)
        self$disc$external <- ifelse(!is.null(discs$external),
                                     discs$external, 0)
      }
      self$beta_u_cov <- function(x) rep(0, length(self$beta_mu))
      if (!is.null(data)) {
        self$in_data <- data.matrix(eval_funcs(scale_input,
                                               data[,names(self$ranges)],
                                               self$ranges))
        self$out_data <- data[, !names(data) %in% names(self$ranges)]
      }
      if (!is.null(self$in_data)) {
        if (is.numeric(self$u_sigma))
          d_corr <- self$u_sigma^2 *
            self$corr$get_corr(self$in_data, actives = self$active_vars)
        else
          d_corr <- diag(
            apply(
              self$in_data, 1,
              self$u_sigma)^2, nrow = nrow(self$in_data)) %*%
            self$corr$get_corr(self$in_data, actives = self$active_vars)
        private$data_corrs <- tryCatch(
          private$data_corrs <- chol2inv(chol(d_corr)),
          error = function(e) {
            MASS::ginv(d_corr)
          }
        )
        private$design_matrix <- t(
          apply(
            self$in_data, 1, function(x)
              purrr::map_dbl(self$basis_f, purrr::exec, x)))
        if (nrow(private$design_matrix) == 1)
          private$design_matrix <- t(private$design_matrix)
        private$u_var_modifier <- private$data_corrs %*%
          private$design_matrix %*% self$beta_sigma %*%
          t(private$design_matrix) %*% private$data_corrs
        private$u_exp_modifier <- private$data_corrs %*%
          (self$out_data - private$design_matrix %*% self$beta_mu)
        private$beta_u_cov_modifier <- self$beta_sigma %*%
          t(private$design_matrix) %*% private$data_corrs
      }
    },
    get_exp = function(x, include_c = TRUE, c_data = NULL) {
      if (nrow(x) > 2000) {
        k <- ceiling(nrow(x)/2000)
        m <- ceiling(nrow(x)/k)
        s_df <- split(x, rep(1:k, each = m, length.out = nrow(x)))
        return(unlist(purrr::map(s_df, ~self$get_exp(., include_c, c_data)),
                      use.names = FALSE))
      }
      x <- x[, names(self$ranges)[names(self$ranges) %in% names(x)]]
      x <- eval_funcs(scale_input, x, self$ranges)
      if (!all(self$beta_sigma == 0) || is.null(self$model)) {
        g <- t(
          apply(
            x, 1, function(y) purrr::map_dbl(self$basis_f, purrr::exec, y)))
        if (length(self$beta_mu) == 1) beta_part <- g * self$beta_mu
        else beta_part <- g %*% self$beta_mu
      }
      else {
        beta_part <- predict(self$model, x)
      }
      x <- data.matrix(x)
      bu <- t(apply(x, 1, self$beta_u_cov))
      u_part <- apply(x, 1, self$u_mu)
      if (!is.null(self$in_data)) {
        if (is.null(c_data))
          c_data <- self$corr$get_corr(self$in_data, x, self$active_vars)
        if (is.numeric(self$u_sigma)) {
          if (length(self$beta_mu) == 1)
            u_part <- t(u_part + (t(bu) %*% t(private$design_matrix) +
                                    self$u_sigma^2 * c_data) %*%
                          private$u_exp_modifier)
          else
            u_part <- u_part + (bu %*% t(private$design_matrix) +
                                  self$u_sigma^2 * c_data) %*%
              private$u_exp_modifier
        }
        else {
          if (length(self$beta_mu) == 1)
            u_part <- t(u_part + (t(bu) %*% t(private$design_matrix) +
                                    sweep(
                                      sweep(
                                        c_data, 2,
                                        apply(
                                          self$in_data, 1,
                                          self$u_sigma), "*"), 1,
                                      apply(x, 1, self$u_sigma), "*")) %*%
                          private$u_exp_modifier)
          else
            u_part <- u_part + (bu %*% t(private$design_matrix) +
                                  sweep(
                                    sweep(
                                      c_data, 2,
                                      apply(self$in_data, 1,
                                            self$u_sigma), "*"), 1,
                                    apply(x, 1, self$u_sigma), "*")) %*%
              private$u_exp_modifier
        }
      }
      if (length(self$beta_mu) == 1) {
        if (include_c) return(c(beta_part + u_part))
        return(c(beta_part))
      }
      if (include_c)
        return(beta_part + u_part)
      return(beta_part)
    },
    get_cov = function(x, xp = NULL, full = FALSE, include_c = TRUE, c_x = NULL, c_xp = NULL) {
      beta_part <- 0
      if (nrow(x) > 2000 && !full) {
        k <- ceiling(nrow(x)/2000)
        m <- ceiling(nrow(x)/k)
        s_df <- split(x, rep(1:k, each = m, length.out = nrow(x)))
        return(unlist(purrr::map(s_df, ~self$get_cov(., xp, FALSE, include_c, c_x, c_xp)),
                      use.names = FALSE))
      }
      x <- eval_funcs(scale_input,
                      x[, names(self$ranges)[names(self$ranges) %in% names(x)]],
                      self$ranges)
      if (!all(self$beta_sigma == 0))
        g_x <- apply(x, 1, function(y)
          purrr::map_dbl(self$basis_f, purrr::exec, y))
      else g_x <- NULL
      x <- data.matrix(x)
      bupart_x <- apply(x, 1, self$beta_u_cov)
      null_flag <- FALSE
      if (is.null(xp)) {
        null_flag <- TRUE
        xp <- x
        g_xp <- g_x
        bupart_xp <- bupart_x
      }
      else {
        xp <- eval_funcs(
          scale_input,
          xp[, names(self$ranges)[names(self$ranges) %in% names(xp)]],
          self$ranges)
        if (!all(self$beta_sigma == 0))
          g_xp <- apply(xp, 1,
                      function(y) purrr::map_dbl(self$basis_f, purrr::exec, y))
        else g_xp <- NULL
        xp <- data.matrix(xp)
        bupart_xp <- apply(xp, 1, self$beta_u_cov)
      }
      if (full || nrow(x) != nrow(xp)) {
        x_xp_c <- self$corr$get_corr(xp, x, self$active_vars)
        if (!is.null(g_x)) {
          if (is.null(nrow(g_x))) beta_part <- g_x %*% self$beta_sigma %*% g_xp
          else beta_part <- t(g_x) %*% self$beta_sigma %*% g_xp
        }
        if (is.numeric(self$u_sigma))
          u_part <- self$u_sigma^2 * x_xp_c
        else
          u_part <- sweep(
            sweep(
              x_xp_c, 2, apply(xp, 1, self$u_sigma), "*"), 1,
            apply(x, 1, self$u_sigma), "*")
        if (!is.null(self$in_data)) {
          if (is.null(c_x))
            c_x <- self$corr$get_corr(self$in_data, x, self$active_vars)
          if (is.null(c_xp))
            c_xp <- if(null_flag) c_x else self$corr$get_corr(self$in_data,
                                                            xp,
                                                            self$active_vars)
          # if(nrow(x) == 1) {
          #   c_x <- t(c_x)
          #   c_xp <- t(c_xp)
          # }
          if (is.numeric(self$u_sigma)) {
            u_part <- u_part - self$u_sigma^4 * c_x %*%
              (private$data_corrs - private$u_var_modifier) %*% t(c_xp)
            if (!all(private$beta_u_cov_modifier == 0)) {
              bupart_x <- bupart_x - private$beta_u_cov_modifier %*%
                (private$design_matrix %*% bupart_x + self$u_sigma^2 * t(c_x))
              bupart_xp <- if (null_flag) bupart_x else bupart_xp -
                private$beta_u_cov_modifier %*%
                (private$design_matrix %*% bupart_xp + self$u_sigma^2 * t(c_xp))
            }
          }
          else {
            c_x <- sweep(
              sweep(
                matrix(
                  c_x, nrow = nrow(x)), 2,
                apply(self$in_data, 1, self$u_sigma), "*"), 1,
              apply(x, 1, self$u_sigma), "*")
            c_xp <- sweep(
              sweep(
                matrix(
                  c_xp, nrow = nrow(xp)), 2,
                apply(self$in_data, 1, self$u_sigma), "*"), 1,
              apply(xp, 1, self$u_sigma), "*")
            u_part <- u_part - c_x %*% (private$data_corrs -
                                          private$u_var_modifier) %*% t(c_xp)
            if (!all(private$beta_u_cov_modifier == 0)) {
              bupart_x <- bupart_x - private$beta_u_cov_modifier %*%
                (private$design_matrix %*% bupart_x + t(c_x))
              bupart_xp <- if(null_flag) bupart_x else bupart_xp -
                private$beta_u_cov_modifier %*%
                (private$design_matrix %*% bupart_xp + t(c_xp))
            }
          }
          if (is.null(nrow(g_x))) {
            bupart_x <- t(bupart_x)
            bupart_xp <- t(bupart_xp)
          }
        }
        if (!all(bupart_x == 0)) {
          if (is.null(nrow(g_x))) bupart <- outer(g_x, bupart_xp, "*") +
              outer(bupart_x, g_xp, "*")
          else bupart <- t(g_x) %*% bupart_xp + t(bupart_x) %*% g_xp
        }
        else {
          bupart <- 0
        }
      }
      else {
        point_seq <- seq_len(nrow(x))
        if (!is.null(g_x)) {
          if (is.null(nrow(g_x)))
            beta_part <- diag(diag(self$beta_sigma) * outer(g_x, g_xp))
          else
            beta_part <- purrr::map_dbl(point_seq,
                                        ~g_x[,.] %*% self$beta_sigma %*% g_xp[,.])
        }
        if (identical(x, xp))
        {
          if (is.numeric(self$u_sigma)) u_part <- rep(self$u_sigma^2, length = nrow(x))
          else u_part <- purrr::map_dbl(point_seq, ~self$u_sigma(x[.,])^2)
        }
        else {
          if (is.numeric(self$u_sigma))
            u_part <- self$u_sigma^2 * diag(
              self$corr$get_corr(x, xp, self$active_vars))
          else
            u_part <- purrr::map_dbl(
              point_seq, ~self$u_sigma(x[.,]) * self$u_sigma(xp[.,])) *
              diag(self$corr$get_corr(x, xp, self$active_vars))
        }
        if (!is.null(self$in_data)) {
          if (is.null(c_x))
            c_x <- self$corr$get_corr(self$in_data, x, self$active_vars)
          if (is.null(c_xp))
            c_xp <- if(null_flag) c_x else self$corr$get_corr(self$in_data,
                                                            xp,
                                                            self$active_vars)
          # if (nrow(x) == 1) {
          #   c_x <- t(c_x)
          #   c_xp <- t(c_xp)
          # }
          if (is.numeric(self$u_sigma)) {
            c_x <- self$u_sigma^2 * c_x
            c_xp <- self$u_sigma^2 * c_xp
          }
          else {
            c_x <- sweep(
              sweep(
                matrix(
                  c_x, nrow = nrow(x)), 2, apply(
                    self$in_data, 1, self$u_sigma), "*"), 1,
              apply(x, 1, self$u_sigma), "*")
            c_xp <- sweep(
              sweep(
                matrix(
                  c_xp, nrow = nrow(xp)), 2, apply(
                    self$in_data, 1, self$u_sigma), "*"), 1,
              apply(xp, 1, self$u_sigma), "*")
          }
          u_part <- u_part -
            rowSums((c_x %*%
                       (private$data_corrs - private$u_var_modifier)) * c_xp)
          if (!all(private$beta_u_cov_modifier == 0)) {
            bupart_x <- bupart_x -
              private$beta_u_cov_modifier %*%
              (private$design_matrix %*% bupart_x + t(c_x))
            bupart_xp <- if(null_flag) bupart_x else bupart_xp -
              private$beta_u_cov_modifier %*%
              (private$design_matrix %*% bupart_xp + t(c_xp))
          }
        }
        if (!all(bupart_x == 0)) {
          if(is.null(nrow(g_x)))
            bupart <- purrr::map_dbl(
              point_seq,
              ~t(purrr::map_dbl(
                self$basis_f,
                purrr::exec, x[.,])) *
                bupart_xp[.] +
                t(bupart_x[.] * purrr::map_dbl(
                  self$basis_f, purrr::exec, xp[.,])))
          else bupart <- purrr::map_dbl(
            point_seq,
            ~t(purrr::map_dbl(
              self$basis_f, purrr::exec, x[.,])) %*% bupart_xp[,.] +
              t(bupart_x[,.] %*%
                  purrr::map_dbl(self$basis_f, purrr::exec, xp[.,])))
        }
        else bupart <- 0
      }
      ## I don't like this, but there's a lot of rounding error going on
      if (include_c)
        result <- round(beta_part + u_part + bupart, 10)
      else
        result <- round(beta_part, 10)
      result[result < 0] <- 1e-6
      return(result)
    },
    ## This derivative form assumes u_sigma is constant.
    get_exp_d = function(x, p) {
      if (nrow(x) > 2000) {
        k <- ceiling(nrow(x)/2000)
        m <- ceiling(nrow(x)/k)
        s_df <- split(x, rep(1:k, each = m, length.out = nrow(x)))
        return(unlist(purrr::map(
          s_df, ~self$get_exp_d(., p)), use.names = FALSE))
      }
      deriv_func <- tryCatch(
        get(paste0(self$corr$corr_name, "_d")),
        error = function(e) return(NULL)
      )
      if(is.null(deriv_func))
        stop("No derivative expression available for correlation function.")
      range_scale <- diff(self$ranges[[p]])/2
      x <- x[, names(self$ranges)[names(self$ranges) %in% names(x)]]
      x <- eval_funcs(scale_input, x, self$ranges)
      if (!p %in% names(self$ranges)) return(0)
      p_ind <- which(names(self$ranges) == p)
      if (is.null(self$model_terms)) {
        model_terms <- purrr::map_chr(
          self$basis_f, function_to_names, names(self$ranges), FALSE)
        g_d <- purrr::map(model_terms, ~D(parse(text = .), p))
      }
      else {
        model_terms <- self$model_terms
        g_d <- purrr::map(
          self$model_terms,
          ~D(parse(text = sub("I\\((\\w*\\^\\d*)\\)", "\\1",
                              gsub(":", "*", .))), p))
      }
      g <- matrix(
        unlist(
          do.call(
            'rbind',
            purrr::map(
              seq_along(x[,1]),
              function(y) purrr::map(g_d, ~eval(., envir = x[y,]))))),
        nrow = nrow(x))
      g <- g/range_scale
      x <- data.matrix(x)
      bu <- t(apply(x, 1, self$beta_u_cov))
      if (length(self$beta_mu) == 1) beta_part <- g * self$beta_mu
      else beta_part <- g %*% self$beta_mu
      u_part <- apply(x, 1, self$u_mu)
      if (!is.null(self$in_data)) {
        c_data <- t(self$corr$get_corr_d(x, self$in_data,
                                         p_ind, actives = self$active_vars))
        c_data <- matrix(c_data/range_scale, nrow = nrow(x))
        if (is.numeric(self$u_sigma)) {
          if (length(self$beta_mu) == 1)
            u_part <- c(t(u_part + (t(bu) %*%
                                      t(private$design_matrix) +
                                      self$u_sigma^2 * c_data) %*%
                            private$u_exp_modifier))
          else
            u_part <- u_part + (bu %*% t(private$design_matrix) +
                                  self$u_sigma^2 * c_data) %*%
              private$u_exp_modifier
        }
        else {
          if (length(self$beta_mu) == 1)
            u_part <- c(t(u_part + (t(bu) %*%
                                      t(private$design_matrix) +
                                      sweep(
                                        sweep(
                                          c_data, 2,
                                          apply(
                                            self$in_data, 1,
                                            self$u_sigma), "*"), 1,
                                        apply(x, 1, self$u_sigma), "*")) %*%
                            private$u_exp_modifier))
          else
            u_part <- u_part +
              (bu %*% t(private$design_matrix) + sweep(
                sweep(
                  c_data, 2,
                  apply(self$in_data, 1,
                        self$u_sigma), "*"), 1,
                apply(x, 1, self$u_sigma), "*")) %*% private$u_exp_modifier
        }
      }
      if (length(self$beta_mu) == 1) return(c(beta_part + u_part))
      return(beta_part + u_part)
    },
    get_cov_d = function(x, p1, xp = NULL, p2 = NULL, full = FALSE) {
      if (nrow(x) > 2000 && !full) {
        k <- ceiling(nrow(x)/2000)
        m <- ceiling(nrow(x)/k)
        s_df <- split(x, rep(1:k, each = m, length.out = nrow(x)))
        return(unlist(purrr::map(
          s_df, ~self$get_cov_d(., p1, xp, p2, FALSE)), use.names = FALSE))
      }
      deriv_func <- tryCatch(
        get(paste0(self$corr$corr_name, "_d")),
        error = function(e) return(NULL)
      )
      if(is.null(deriv_func))
        stop("No derivative expression available for correlation function.")
      range_scale1 <- diff(self$ranges[[p1]])/2
      x <- x[, names(self$ranges)[names(self$ranges) %in% names(x)]]
      x <- eval_funcs(scale_input, x, self$ranges)
      if (is.null(p2)) p2 <- p1
      range_scale2 <- diff(self$ranges[[p2]])/2
      if (!p1 %in% names(self$ranges) || !p2 %in% names(self$ranges)) {
        if (!full || is.null(xp)) return(rep(0, nrow(x)))
        return(matrix(0, nrow = nrow(x), ncol = nrow(xp)))
      }
      p1_ind <- which(names(self$ranges) == p1)
      p2_ind <- which(names(self$ranges) == p2)
      if (is.null(self$model_terms)) {
        model_terms <- purrr::map_chr(self$basis_f,
                                      function_to_names,
                                      names(self$ranges), FALSE)
        g_d1 <- purrr::map(model_terms, ~D(parse(text = .), p1))
        g_d2 <- purrr::map(model_terms, ~D(parse(text = .), p2))
      }
      else {
        model_terms <- self$model_terms
        g_d1 <- purrr::map(self$model_terms,
                           ~D(parse(text = sub("I\\((\\w*\\^\\d*)\\)",
                                               "\\1", gsub(":", "*", .))), p1))
        g_d2 <- purrr::map(self$model_terms,
                           ~D(parse(text = sub("I\\((\\w*\\^\\d*)\\)",
                                               "\\1", gsub(":", "*", .))), p2))
      }
      g_x <- t(
        matrix(
          unlist(
            do.call(
              'rbind', purrr::map(
                seq_along(x[,1]),
                function(y) purrr::map(g_d1,
                                       ~eval(., envir = x[y,]))))),
          nrow = nrow(x)))
      gp_x <- t(
        matrix(
          unlist(
            do.call(
              'rbind', purrr::map(
                seq_along(x[,1]),
                function(y) purrr::map(g_d2,
                                       ~eval(., envir = x[y,]))))),
          nrow = nrow(x)))
      x <- data.matrix(x)
      bupart_x <- apply(x, 1, self$beta_u_cov)
      null_flag <- FALSE
      if (is.null(xp)) {
        null_flag <- TRUE
        xp <- x
        bupart_xp <- bupart_x
      }
      else {
        xp <- eval_funcs(
          scale_input,
          xp[, names(self$ranges)[names(self$ranges) %in% names(xp)]],
          self$ranges)
        gp_x <- t(
          matrix(
            unlist(
              do.call(
                'rbind',
                purrr::map(
                  seq_along(xp[,1]),
                  function(y) purrr::map(g_d2,
                                         ~eval(., envir = xp[y,]))))),
            nrow = nrow(xp)))
        xp <- data.matrix(xp)
        bupart_xp <- apply(xp, 1, self$beta_u_cov)
      }
      g_x <- g_x/range_scale1
      gp_x <- gp_x/range_scale2
      if (full || nrow(x) != nrow(xp)) {
          x_xp_c <- t(self$corr$get_corr_d(x, xp, p1_ind, p2_ind,
                                           self$active_vars))/
            (range_scale1*range_scale2)
        if (is.null(nrow(g_x))) beta_part <- g_x %*% self$beta_sigma %*% gp_x
        else beta_part <- t(g_x) %*% self$beta_sigma %*% gp_x
        if (is.numeric(self$u_sigma))
          u_part <- self$u_sigma^2 * x_xp_c
        else
          u_part <- sweep(
            sweep(
              x_xp_c, 2,
              apply(xp, 1,
                    self$u_sigma), "*"), 1,
            apply(x, 1, self$u_sigma), "*")
        if (!is.null(self$in_data)) {
          c_x <- t(self$corr$get_corr_d(x, self$in_data,
                                        p1_ind, p2 = NULL, self$active_vars))/
            range_scale1
          c_xp <- -1*t(self$corr$get_corr_d(xp, self$in_data,
                                      p2_ind, p2 = NULL, self$active_vars))/
            range_scale2
          if (is.numeric(self$u_sigma)) {
            c_x <- self$u_sigma^2 * c_x
            c_xp <- self$u_sigma^2 * c_xp
          }
          else {
            c_x <- sweep(
              sweep(
                matrix(
                  c_x, nrow = nrow(x)), 2,
                apply(self$in_data, 1,
                      self$u_sigma), "*"), 1,
              apply(x, 1, self$u_sigma), "*")
            c_xp <- sweep(
              sweep(
                matrix(
                  c_xp, nrow = nrow(xp)), 2,
                apply(self$in_data, 1, self$u_sigma), "*"), 1,
              apply(xp, 1, self$u_sigma), "*")
          }
          u_part <- u_part - c_x %*%
            (private$data_corrs - private$u_var_modifier) %*% t(c_xp)
          bupart_x <- bupart_x -
            private$beta_u_cov_modifier %*%
            (private$design_matrix %*% bupart_x + t(c_x))
          bupart_xp <- if(null_flag) bupart_x else bupart_xp -
            private$beta_u_cov_modifier %*%
            (private$design_matrix %*% bupart_xp + t(c_xp))
          if (is.null(nrow(g_x))) {
            bupart_x <- t(bupart_x)
            bupart_xp <- t(bupart_xp)
          }
        }
        if (is.null(nrow(g_x))) bupart <- outer(g_x, bupart_xp, "*") +
            outer(bupart_x, gp_x, "*")
        else bupart <- t(g_x) %*% bupart_xp + t(bupart_x) %*% gp_x
      }
      else {
        point_seq <- seq_len(nrow(x))
        if (is.null(nrow(g_x))) beta_part <- diag(diag(self$beta_sigma) *
                                                    outer(g_x, gp_x))
        else beta_part <- purrr::map_dbl(
          point_seq, ~g_x[,.] %*% self$beta_sigma %*% gp_x[,.])
        if (is.numeric(self$u_sigma))
          u_part <- self$u_sigma^2 * diag(
            self$corr$get_corr_d(x, xp, p1_ind, p2_ind, self$active_vars))/
            (range_scale1*range_scale2)
        else u_part <- purrr::map_dbl(
          point_seq,
          ~self$u_sigma(x[.,])*self$u_sigma(x[.,])) *
            diag(self$corr$get_corr_d(x, xp, p1_ind,
                                      p2_ind, self$active_vars))/
            (range_scale1*range_scale2)
        if (!is.null(self$in_data)) {
          c_x <- t(self$corr$get_corr_d(x, self$in_data,
                                        p1_ind, p2 = NULL, self$active_vars))/
            range_scale1
          c_xp <- -1*t(self$corr$get_corr_d(xp, self$in_data, p2_ind,
                                      p2 = NULL, self$active_vars))/
            range_scale2
          if (is.numeric(self$u_sigma)) {
            c_x <- self$u_sigma^2 * c_x
            c_xp <- self$u_sigma^2 * c_xp
          }
          else {
            c_x <- sweep(
              sweep(
                matrix(c_x, nrow = nrow(x)), 2,
                apply(self$in_data, 1, self$u_sigma), "*"), 1,
              apply(x, 1, self$u_sigma), "*")
            c_xp <- sweep(
              sweep(
                matrix(c_xp, nrow = nrow(xp)), 2,
                apply(self$in_data, 1, self$u_sigma), "*"), 1,
              apply(xp, 1, self$u_sigma), "*")
          }
          u_part <- u_part -
            rowSums((c_x %*%
                       (private$data_corrs - private$u_var_modifier)) * c_xp)
          bupart_x <- bupart_x -
            private$beta_u_cov_modifier %*%
            (private$design_matrix %*% bupart_x + t(c_x))
          bupart_xp <- if(null_flag) bupart_x else bupart_xp -
            private$beta_u_cov_modifier %*%
            (private$design_matrix %*% bupart_xp + t(c_xp))
        }
        if(is.null(nrow(g_x)))
          bupart <- purrr::map_dbl(
            point_seq,
            ~t(purrr::map_dbl(
              self$basis_f, purrr::exec, x[.,])) *
              bupart_xp[.] + t(bupart_x[.] * purrr::map_dbl(
                self$basis_f, purrr::exec, xp[.,])))
        else bupart <- purrr::map_dbl(
          point_seq,
          ~t(purrr::map_dbl(self$basis_f, purrr::exec, x[.,])) %*%
            bupart_xp[,.] + t(bupart_x[,.] %*% purrr::map_dbl(
              self$basis_f, purrr::exec, xp[.,])))
      }
      return(round(beta_part + u_part + bupart, 10))
    },
    implausibility = function(x, z, cutoff = NULL) {
      if (is.null(nrow(x))) x <- setNames(
        data.frame(matrix(x, ncol = 1)), names(self$ranges))
      if (nrow(x) > 2000) {
        k <- ceiling(nrow(x)/2000)
        m <- ceiling(nrow(x)/k)
        s_df <- split(x, rep(1:k, each = m, length.out = nrow(x)))
        return(unlist(
          purrr::map(
            s_df, ~self$implausibility(., z = z, cutoff = cutoff)),
          use.names = FALSE))
      }
      temp_scale_x <- x[, names(self$ranges)[names(self$ranges) %in% names(x)]]
      temp_scale_x <- eval_funcs(scale_input, temp_scale_x, self$ranges)
      temp_scale_x <- data.matrix(temp_scale_x)
      corr_x <- self$corr$get_corr(self$in_data, temp_scale_x, self$active_vars)
      disc_quad <- sum(purrr::map_dbl(self$disc, ~.^2))
      if (!is.numeric(z) && !is.null(z$val)) {
        imp_var <- self$get_cov(x, c_x = corr_x, c_xp = corr_x) + z$sigma^2 + disc_quad
        imp <- sqrt((z$val - self$get_exp(x, c_data = corr_x))^2/imp_var)
      }
      else {
        pred <- self$get_exp(x, c_data = corr_x)
        bound_check <- tryCatch(
            {purrr::map_dbl(pred, function(y) {
            if (y <= z[2] && y >= z[1]) return(0)
            if (y < z[1]) return(-1)
            if (y > z[2]) return(1)
          })},
          error = function(e) {
            stop(
              paste("Problem with form of targets: perhaps target missing?", e)
              )
          }
        )
        which_compare <- purrr::map_dbl(bound_check, function(y) {
          if (y < 1) return(z[1])
          return(z[2])
        })
        uncerts <- self$get_cov(x, c_x = corr_x, c_xp = corr_x) + disc_quad
        uncerts[uncerts <= 0] <- 0.0001
        imp <- bound_check * (pred - which_compare)/sqrt(uncerts)
      }
      if (is.null(cutoff)) return(imp)
      return(imp <= cutoff)
    },
    adjust = function(data, out_name) {
      this_data_in <- data.matrix(eval_funcs(
        scale_input, data[,names(self$ranges)], self$ranges))
      this_data_out <- data[,out_name]
      if (all(eigen(self$beta_sigma)$values == 0)) {
        new_beta_var <- self$beta_sigma
        new_beta_exp <- self$beta_mu
      }
      else {
        G <- apply(this_data_in, 1,
                   function(x) purrr::map_dbl(self$basis_f, purrr::exec, x))
        Ot <- self$corr$get_corr(this_data_in, actives = self$active_vars)
        O <- tryCatch(
          chol2inv(chol(Ot)),
          error = function(x) {
            MASS::ginv(Ot)
          }
        )
        siginv <- tryCatch(
          chol2inv(chol(self$beta_sigma)),
          error = function(e) {
            MASS::ginv(self$beta_sigma)
          }
        )
        new_beta_var <- tryCatch(
          chol2inv(chol(G %*% O %*% (if(is.null(nrow(G))) G else t(G)) + siginv)),
          error = function(e) {
            MASS::ginv(G %*% O %*% (if(is.null(nrow(G))) G else t(G)) + siginv)
          }
        )
        new_beta_exp <- new_beta_var %*%
          (siginv %*% self$beta_mu + G %*% O %*% this_data_out)
      }
      new_em <- Emulator$new(self$basis_f,
                             beta = list(mu = new_beta_exp, sigma = new_beta_var),
                             u = list(sigma = self$u_sigma, corr = self$corr),
                             ranges = self$ranges,
                             data = data[, c(names(self$ranges), out_name)],
                             original_em = self, out_name = out_name,
                             model = self$model, a_vars = self$active_vars,
                             discs = self$disc, multiplier = self$multiplier)
      return(new_em)
    },
    set_sigma = function(sigma) {
      if (is.null(self$o_em)) {
        new_em <- self$clone()
        new_em$u_sigma <- sigma
        return(new_em)
      }
      new_o_em <- self$o_em$clone()
      new_o_em$u_sigma <- sigma
      dat <- setNames(
        data.frame(
          cbind(
            eval_funcs(scale_input,
                       data.frame(self$in_data), self$ranges, FALSE),
            self$out_data)), c(names(self$ranges), self$output_name))
      return(new_o_em$adjust(dat, self$output_name))
    },
    mult_sigma = function(m) {
      if (is.null(self$o_em)) {
        new_em <- self$clone()
        new_em$multiplier <- new_em$multiplier * m
        return(new_em)
      }
      new_o_em <- self$o_em$clone()
      new_o_em$multiplier <- new_o_em$multiplier * m
      dat <- setNames(
        data.frame(
          cbind(
            eval_funcs(scale_input,
                       data.frame(self$in_data), self$ranges, FALSE),
            self$out_data)), c(names(self$ranges), self$output_name))
      return(new_o_em$adjust(dat, self$output_name))
    },
    set_hyperparams = function(hp, nugget = self$corr$nugget) {
      current_u <- self$corr
      if (all(names(current_u$hyper_p) != names(hp)))
        stop("Hyperparameter specification does not match current correlation function.")
      if (is.null(self$o_em)) {
        new_em <- self$clone()
        new_em$corr <- new_em$corr$set_hyper_p(hp, nugget)
        return(new_em)
      }
      new_o_em <- self$o_em$clone()
      new_o_em$corr <- new_o_em$corr$set_hyper_p(hp, nugget)
      dat <- setNames(
        data.frame(
          cbind(
            eval_funcs(
              scale_input,
              data.frame(self$in_data), self$ranges, FALSE), self$out_data)),
        c(names(self$ranges), self$output_name))
      return(new_o_em$adjust(dat, self$output_name))
    },
    print = function(...) {
      cat("Parameters and ranges: ",
          paste(names(self$ranges),
                paste0(purrr::map(self$ranges, round, 4)),
                sep = ": ", collapse = ": "), "\n")
      cat("Specifications: \n")
        if (!is.null(self$model))
          cat("\t Basis functions: ",
              paste0(names(self$model$coefficients), collapse="; "), "\n")
        else if (!is.null(self$o_em$model))
          cat("\t Basis Functions: ",
              paste0(names(self$o_em$model$coefficients), collapse="; "), "\n")
        else
          cat("\t Basis functions: ",
              paste0(purrr::map_chr(
                self$basis_f, function_to_names, names(self$ranges), FALSE),
                collapse = "; "), "\n")
      cat("\t Active variables",
          paste0(names(self$ranges)[self$active_vars], collapse = "; "), "\n")
      cat("\t Regression Surface Expectation: ",
          paste(round(self$beta_mu, 4), collapse = "; "), "\n")
      cat("\t Regression surface Variance (eigenvalues): ",
          paste(round(eigen(self$beta_sigma)$values, 4), collapse = "; "), "\n")
      cat("Correlation Structure: \n")
      if (!is.null(private$data_corrs))
        cat("Bayes-adjusted emulator - prior specifications listed. \n")
      cat("\t Variance (Representative): ",
          if (is.numeric(self$u_sigma))
            self$u_sigma^2
          else
            self$u_sigma(purrr::map_dbl(self$ranges, mean))^2, "\n")
      cat("\t Expectation: ", self$u_mu(rep(0, length(ranges))), "\n")
      self$corr$print(prepend = "\t")
      cat("Mixed covariance: ", self$beta_u_cov(rep(0, length(ranges))), "\n")
    },
    plot = function(...) {
      emulator_plot(self, ...)
    }
  )
)
