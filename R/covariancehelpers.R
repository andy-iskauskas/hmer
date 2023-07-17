EmulatedMatrix <- R6::R6Class(
  "EmulatorMatrix",
  public = list(
    mat = NULL,
    theta_t = NULL,
    rho = NULL,
    logmatrix = FALSE,
    initialize = function(em_mat, thet, rho_mat, logged) {
      self$mat <- em_mat
      self$theta_t <- thet
      self$rho <- rho_mat
      self$logmatrix <- logged
    },
    get_matrix = function() return(self$mat),
    get_exp = function(x, check_positive = TRUE) {
      t_out <- array(t(sapply(seq_along(self$mat), function(i) {
        if (!"Emulator" %in% class(self$mat[[i]])) output <- rep(NA, nrow(x))
        else if (self$mat[[i]]$output_name == "zero_emulator") {
          if (self$logmatrix) output <- rep(-Inf, nrow(x))
          output <- rep(0, nrow(x))
        }
        else if (check_positive && !self$logmatrix) output <- self$mat[[i]]$get_exp(x, check_neg = FALSE)
        else output <- self$mat[[i]]$get_exp(x)
      })), c(dim(self$mat), nrow(x)))
      t_out[is.na(t_out)] <- 0
      for (i in 1:nrow(x)) {
        t_out[,,i] <- t_out[,,i] + t(t_out[,,i]) - diag(diag(t_out[,,i]))
      }
      if (check_positive && !self$logmatrix) {
        t_out <- array(apply(t_out, 3, function(mat) {
          es <- eigen(mat)
          eve <- es$vectors
          eva <- es$values
          eva[eva < 0] <- 0
          eve %*% diag(eva) %*% t(eve)
        }), dim = c(dim(self$mat), nrow(x)))
      }
      if (self$logmatrix) {
        t_out <- array(apply(t_out, 3, function(mat) {
          es <- eigen(mat)
          eve <- es$vectors
          eva <- es$values
          eva <- exp(eva)
          eve %*% diag(eva) %*% t(eve)
        }), dim = c(dim(self$mat), nrow(x)))
      }
      return(t_out)
    },
    get_cov = function(x, check_positive = TRUE) {
      t_out <- array(t(sapply(seq_along(self$mat), function(i) {
        if (!"Emulator" %in% class(self$mat[[i]])) output <- rep(NA, nrow(x))
        else if (self$mat[[i]]$output_name == "zero_emulator") {
          output <- rep(0, nrow(x))
        }
        else if (check_positive) output <- self$mat[[i]]$get_cov(x, check_neg = FALSE)
        else output <- self$mat[[i]]$get_cov(x)
      })), c(dim(self$mat), nrow(x)))
      t_out[is.na(t_out)] <- 0
      for (i in 1:nrow(x)) {
        t_out[,,i] <- t_out[,,i] + t(t_out[,,i]) - diag(diag(t_out[,,i]))
      }
      if (check_positive) {
        t_out <- array(apply(t_out, 3, function(mat) {
          es <- eigen(mat)
          eve <- es$vectors
          eva <- es$values
          eva[eva < 0] <- 0
          eve %*% diag(eva) %*% t(eve)
        }), dim = c(dim(self$mat), nrow(x)))
      }
      return(t_out)
    },
    get_uncertainty = function(x, mean_ems, check_positive = TRUE) {
      prior_uncert <- self$get_exp(x, check_positive)
      cor_mat <- array(apply(prior_uncert, 3, cov2cor), dim = c(dim(self$mat), nrow(x)))
      e_mat <- matrix(nrow = length(mean_ems), ncol = length(mean_ems))
      e_mat[row(e_mat) == col(e_mat)] <- mean_ems
      e_mat <- matrix(e_mat, nrow = length(mean_ems), ncol = length(mean_ems))
      cov_out <- array(t(sapply(seq_along(e_mat), function(i) {
        if ("logical" %in% class(e_mat[[i]])) return(rep(0, nrow(x)))
        return(as.numeric(e_mat[[i]]$get_cov(x)))
      })), dim = c(dim(e_mat), nrow(x)))
      return(array(sapply(seq_len(nrow(x)), function(i) {
        sqrt(cov_out[,,i]) %*% cor_mat[,,i] %*% sqrt(cov_out[,,i])
      }), dim = c(dim(e_mat), nrow(x))))
    },
    print = function(...) {
      prior_var_matrix <- matrix(sapply(self$mat, function(x) {
        if (!"Emulator" %in% class(x)) return(0)
        if (x$output_name == "zero_emulator") return(0)
        return(x$u_sigma^2)
      }), nrow = nrow(self$mat))
      rownames(prior_var_matrix) <- colnames(prior_var_matrix) <- purrr::map_chr(diag(self$mat), ~.$output_name)
      prior_var_matrix <- prior_var_matrix + t(prior_var_matrix) - diag(diag(prior_var_matrix))
      cat("Emulated covariance matrix\n")
      cat(paste0("Outputs: ", paste0(rownames(prior_var_matrix), collapse = ", "), "\n"))
      cat("Prior uncertainty:\n")
      print(prior_var_matrix)
      cat("Between-output prior structure:\n")
      cat(paste("theta_t:", round(self$theta_t, 3), "\n"))
      cat("Correlation matrix rho:\n")
      print(round(self$rho, 3))
      if (self$logmatrix) cat("(Log-covariance emulators)")
      invisible(self)
    }
  )
)

partition_by_output <- function(data, out_names, by_time = TRUE, return_label = FALSE) {
  input_names <- names(data)[!names(data) %in% out_names]
  out_times <- as.numeric(sub(".*[^\\d](\\d+)$", "\\1", out_names, perl = TRUE))
  out_labels <- sub("(.*[^\\d])\\d+$", "\\1", out_names, perl = TRUE)
  if (by_time) {
    unique_times <- unique(out_times)
    if (return_label)
      partitioned <- purrr::map(unique_times, function(tm) {
        out_names[which(out_times == tm)]
      })
    else
      partitioned <- purrr::map(unique_times, function(tm) {
        data[, c(input_names, out_names[which(out_times == tm)])]
      })
  }
  else {
    unique_labels <- unique(out_labels)
    if (return_label)
      partitioned <- purrr::map(unique_labels, function(lb) {
        out_names[which(out_labels == lb)]
      })
    else
      partitioned <- purrr::map(unique_labels, function(lb) {
        data[, c(input_names, out_names[which(out_labels == lb)])]
      })
  }
  return(partitioned)
}
#' Estimate rho matrix
#'
#' @importFrom stats var
#' @importFrom dplyr group_by across
#'
#' @keywords internal
#' @noRd
get_mpc_rho_est <- function(data, out_names, ...) {
  part_data <- partition_by_output(data, out_names, ...)
  part_covs <- purrr::map(part_data, function(dat) {
    g_by_point <- dat |> group_by(across(all_of(names(data)[!names(data) %in% out_names])))
    vars <- do.call('rbind.data.frame', purrr::map(dplyr::group_rows(g_by_point), function(gp) {
      apply(g_by_point[gp,names(dat) %in% out_names], 2, var)
    }))
    vars <- setNames(vars, names(dat)[names(dat) %in% out_names])
    return(cor(vars))
  })
  rho_mat_lengths <- unique(purrr::map_dbl(part_covs, nrow))
  rho_mat_subset <- purrr::map(seq_len(max(rho_mat_lengths)), function(l) {
    part_covs[purrr::map_lgl(part_covs, ~nrow(.) == l)]
  })
  rho_mat_subset <- rho_mat_subset[purrr::map_lgl(rho_mat_subset, ~length(.) > 0)]
  rho_mat <- purrr::map(rho_mat_subset, ~Reduce("+", .)/length(.))
  out_mat <- rho_mat[[1]]
  if (length(rho_mat) > 1) {
    for (i in 2:length(rho_mat)) {
      out_mat <- rbind(cbind(out_mat, matrix(0, nrow = nrow(out_mat), ncol = ncol(rho_mat[[i]]))),
                       cbind(matrix(0, nrow = nrow(rho_mat[[i]]), ncol = ncol(out_mat)), rho_mat[[i]]))
    }
  }
  rownames(out_mat) <- colnames(out_mat) <- unique(sub("(.*[^\\d])\\d+$", "\\1", out_names, perl = TRUE))
  return(out_mat)
}
# Estimate theta-time value
get_mpc_theta_est <- function(data, out_names, ems, rho = NULL) {
  if (is.null(rho)) rho <- get_mpc_rho_est(data, out_names)
  part_indices <- partition_by_output(data, out_names, FALSE, TRUE)
  g_by_point <- data |> group_by(across(all_of(names(data)[!names(data) %in% out_names])))
  vars <- do.call('rbind.data.frame', purrr::map(dplyr::group_rows(g_by_point), function(gp) {
    apply(g_by_point[gp,names(data) %in% out_names], 2, var)
  })) |> setNames(names(data)[names(data) %in% out_names])
  all_covs <- cov(vars)
  all_cors <- cor(vars)
  unique_inputs <- unique(data[,names(data)[!names(data) %in% out_names]])
  em_preds <- sqrt(do.call('cbind.data.frame', purrr::map(ems, ~.$get_cov(unique_inputs))) |>
                     setNames(purrr::map_chr(ems, "output_name")))
  all_vals <- do.call('c', purrr::map(seq_len(nrow(em_preds)), function(k) {
    theta_ests <- c()
    for (i in 1:(ncol(all_cors)-1)) {
      for (j in (i+1):ncol(all_cors)) {
        nm1 <- out_names[i]
        nm2 <- out_names[j]
        group1 <- which(purrr::map_lgl(part_indices, ~nm1 %in% .))
        group2 <- which(purrr::map_lgl(part_indices, ~nm2 %in% .))
        times <- as.numeric(sub(".*[^\\d](\\d+)$", "\\1", c(nm1, nm2), perl = TRUE))
        rho_val <- rho[group1, group2]
        theta_ests <- suppressWarnings(c(theta_ests, sqrt(-diff(times)^2/log(all_cors[i,j]/(rho_val*em_preds[k,i]*em_preds[k,j])))))
      }
    }
    return(theta_ests)
  }))
  all_vals <- all_vals[!is.na(c(all_vals)) & !is.nan(c(all_vals))]
  return(mean(all_vals))
}
