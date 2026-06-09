pv_symmetrize <- function(x) {
  0.5 * (x + t(x))
}

pv_check_numeric_matrix <- function(x, name) {
  if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x))) {
    pv_abort(sprintf("`%s` must be a finite numeric matrix.", name))
  }
  x
}

pv_validate_rubin_u_array <- function(U, M, p) {
  if (is.list(U)) {
    if (length(U) != M) {
      pv_abort("`U` list length must equal the number of imputations.")
    }
    out <- array(NA_real_, dim = c(p, p, M))
    for (m in seq_len(M)) {
      label <- sprintf("U[[%d]]", m)
      Um <- pv_check_numeric_matrix(U[[m]], label)
      if (!identical(dim(Um), c(p, p))) {
        pv_abort("Every `U` matrix must be p x p.")
      }
      if (!is.null(rownames(Um)) && !is.null(colnames(Um)) &&
          !identical(rownames(Um), colnames(Um))) {
        pv_abort("When present, row and column names of each `U` matrix must match.")
      }
      if (max(abs(Um - t(Um))) > 1e-10) {
        pv_abort("Every `U` matrix must be symmetric.")
      }
      Um <- pv_symmetrize(Um)
      eigenvalues <- eigen(Um, symmetric = TRUE, only.values = TRUE)$values
      if (min(eigenvalues) < -1e-10) {
        pv_abort("Every `U` matrix must be positive semidefinite.")
      }
      out[, , m] <- Um
    }
    return(out)
  }

  if (is.array(U) && length(dim(U)) == 3L) {
    if (!identical(dim(U), c(p, p, M))) {
      pv_abort("`U` array must have dimension p x p x M.")
    }
    if (!is.numeric(U) || any(!is.finite(U))) {
      pv_abort("`U` array must be finite numeric.")
    }
    for (m in seq_len(M)) {
      Um <- U[, , m]
      if (max(abs(Um - t(Um))) > 1e-10) {
        pv_abort("Every `U` matrix must be symmetric.")
      }
      Um <- pv_symmetrize(Um)
      eigenvalues <- eigen(Um, symmetric = TRUE, only.values = TRUE)$values
      if (min(eigenvalues) < -1e-10) {
        pv_abort("Every `U` matrix must be positive semidefinite.")
      }
      U[, , m] <- Um
    }
    return(U)
  }

  pv_abort("`U` must be a list of p x p matrices or a p x p x M array.")
}

pv_validate_rubin_rho <- function(rho) {
  if (!is.numeric(rho) || length(rho) < 1L || any(!is.finite(rho)) ||
      any(rho < 0) || any(rho > 1 + sqrt(.Machine$double.eps))) {
    pv_abort("`rho` must be finite and in [0, 1].")
  }
  rho
}

rubin_df_classic <- function(M, rho) {
  M <- pv_assert_scalar_number(M, "M", integer = TRUE, lower = 1)
  rho <- pv_validate_rubin_rho(rho)

  out <- if (M > 1L) {
    (M - 1L) / pmax(rho^2, .Machine$double.eps)
  } else {
    rep(Inf, length(rho))
  }
  names(out) <- names(rho)
  out
}

rubin_df_barnard_rubin <- function(M, rho, df_complete) {
  nu_old <- rubin_df_classic(M, rho)
  rho <- pv_validate_rubin_rho(rho)

  if (missing(df_complete) || is.null(df_complete)) {
    pv_abort("`df_complete` is required when `df_method = \"barnard_rubin\"`.")
  }
  if (!is.numeric(df_complete) || length(df_complete) < 1L ||
      any(is.na(df_complete)) || any(df_complete <= 0)) {
    pv_abort("`df_complete` must be positive.")
  }
  if (length(df_complete) == 1L) {
    df_complete <- rep(df_complete, length(rho))
  }
  if (length(df_complete) != length(rho)) {
    pv_abort("`df_complete` must be scalar or match the length of `rho`.")
  }
  if (any(is.infinite(df_complete))) {
    out <- nu_old
    finite <- is.finite(df_complete)
    if (!any(finite)) {
      return(out)
    }
    out[finite] <- rubin_df_barnard_rubin(M, rho[finite], df_complete[finite])
    return(out)
  }

  nu_obs <- ((df_complete + 1) / (df_complete + 3)) *
    df_complete * pmax(1 - rho, 0)
  1 / (1 / nu_old + 1 / pmax(nu_obs, .Machine$double.eps))
}

rubin_pool_matrix <- function(
  beta,
  U,
  orientation = c("rows_pv", "cols_pv"),
  conf_level = 0.95,
  allow_m1 = FALSE,
  df_method = c("classic", "barnard_rubin"),
  df_complete = NULL
) {
  orientation <- match.arg(orientation)
  df_method <- match.arg(df_method)
  conf_level <- pv_assert_probability(conf_level, "conf_level")
  allow_m1 <- pv_assert_scalar_logical(allow_m1, "allow_m1")
  beta <- pv_check_numeric_matrix(beta, "beta")

  if (orientation == "rows_pv") {
    M <- nrow(beta)
    p <- ncol(beta)
    beta_pxm <- t(beta)
  } else {
    p <- nrow(beta)
    M <- ncol(beta)
    beta_pxm <- beta
  }

  if (M < 2L && !allow_m1) {
    pv_abort("Need M >= 2 imputations unless `allow_m1 = TRUE`.")
  }
  if (p < 1L) {
    pv_abort("`beta` must contain at least one parameter.")
  }
  param_names <- if (orientation == "rows_pv") colnames(beta) else rownames(beta)
  if (is.null(param_names)) {
    param_names <- paste0("theta", seq_len(p))
  }

  U_array <- pv_validate_rubin_u_array(U, M = M, p = p)
  beta_bar <- rowMeans(beta_pxm)
  names(beta_bar) <- param_names
  U_bar <- pv_symmetrize(apply(U_array, c(1, 2), mean))

  if (M > 1L) {
    beta_centered <- sweep(beta_pxm, 1L, beta_bar, FUN = "-")
    B <- (beta_centered %*% t(beta_centered)) / (M - 1L)
  } else {
    B <- matrix(0, p, p)
  }
  B <- pv_symmetrize(B)
  T_MI <- pv_symmetrize(U_bar + (1 + 1 / M) * B)

  lambda <- diag((1 + 1 / M) * B) / pmax(diag(T_MI), .Machine$double.eps)
  riv <- diag((1 + 1 / M) * B) / pmax(diag(U_bar), .Machine$double.eps)
  names(lambda) <- names(riv) <- param_names
  df_complete_for_df <- df_complete
  if (!is.null(df_complete_for_df) && length(df_complete_for_df) > 1L &&
      !is.null(names(df_complete_for_df))) {
    if (anyDuplicated(names(df_complete_for_df)) ||
        !setequal(names(df_complete_for_df), param_names)) {
      pv_abort("`df_complete` names must match parameter names.")
    }
    df_complete_for_df <- df_complete_for_df[param_names]
  }
  df <- if (identical(df_method, "barnard_rubin")) {
    rubin_df_barnard_rubin(M, lambda, df_complete_for_df)
  } else {
    rubin_df_classic(M, lambda)
  }

  se <- sqrt(diag(T_MI))
  alpha <- 1 - conf_level
  t_quantile <- stats::qt(1 - alpha / 2, df = df)

  dimnames(U_bar) <- list(names(beta_bar), names(beta_bar))
  dimnames(B) <- list(names(beta_bar), names(beta_bar))
  dimnames(T_MI) <- list(names(beta_bar), names(beta_bar))
  names(se) <- names(df) <- names(lambda) <- names(riv) <- names(beta_bar)

  if (is.null(df_complete_for_df)) {
    df_complete_out <- rep(NA_real_, p)
  } else if (length(df_complete_for_df) == 1L) {
    df_complete_out <- rep(df_complete_for_df, p)
  } else {
    df_complete_out <- df_complete_for_df
  }
  names(df_complete_out) <- names(beta_bar)

  list(
    beta = beta_bar,
    beta_bar = beta_bar,
    U_bar = U_bar,
    B = B,
    T_MI = T_MI,
    total_var = T_MI,
    se = se,
    df = df,
    df_classic = rubin_df_classic(M, lambda),
    df_method = df_method,
    df_complete = df_complete_out,
    lambda = lambda,
    fmi = lambda,
    riv = riv,
    rho = lambda,
    ci_low = beta_bar - t_quantile * se,
    ci_high = beta_bar + t_quantile * se,
    M = M,
    p = p,
    orientation = orientation,
    conf_level = conf_level
  )
}

rubin_pool <- rubin_pool_matrix

rubin_pool_scalar <- function(
  q,
  u,
  conf_level = 0.95,
  allow_m1 = FALSE,
  df_method = c("classic", "barnard_rubin"),
  df_complete = NULL
) {
  df_method <- match.arg(df_method)
  if (!is.numeric(q) || !is.numeric(u) || length(q) != length(u) ||
      length(q) < 1L || any(!is.finite(q)) || any(!is.finite(u))) {
    pv_abort("`q` and `u` must be finite numeric vectors of equal positive length.")
  }
  if (any(u < 0)) {
    pv_abort("`u` must contain non-negative within-imputation variances.")
  }
  U <- lapply(u, function(x) matrix(x, 1L, 1L))
  out <- rubin_pool_matrix(
    beta = matrix(q, ncol = 1L, dimnames = list(NULL, "theta")),
    U = U,
    orientation = "rows_pv",
    conf_level = conf_level,
    allow_m1 = allow_m1,
    df_method = df_method,
    df_complete = df_complete
  )
  list(
    q_bar = unname(out$beta[[1L]]),
    u_bar = unname(out$U_bar[1L, 1L]),
    b = unname(out$B[1L, 1L]),
    total_var = unname(out$T_MI[1L, 1L]),
    se = unname(out$se[[1L]]),
    df = unname(out$df[[1L]]),
    df_classic = unname(out$df_classic[[1L]]),
    df_method = out$df_method,
    df_complete = unname(out$df_complete[[1L]]),
    ci_low = unname(out$ci_low[[1L]]),
    ci_high = unname(out$ci_high[[1L]]),
    M = out$M
  )
}
