ccc_as_draw_matrix <- function(draws) {
  if (is.matrix(draws)) {
    out <- draws
  } else if (is.data.frame(draws)) {
    out <- as.matrix(draws)
  } else {
    pv_abort("`draws` must be a numeric matrix or data frame.")
  }
  if (!is.numeric(out) || any(!is.finite(out))) {
    pv_abort("`draws` must contain only finite numeric values.")
  }
  if (is.null(colnames(out)) || any(!nzchar(colnames(out)))) {
    pv_abort("`draws` must have non-empty column names.")
  }
  if (anyDuplicated(colnames(out))) {
    pv_abort("`draws` must have unique column names.")
  }
  if (nrow(out) < 2L) {
    pv_abort("`draws` must contain at least two posterior draws.")
  }
  out
}

ccc_forbidden_target_text <- function(target) {
  values <- c(
    target$target_source,
    target$source,
    target$policy,
    target$target_mode,
    target$target_id,
    attr(target, "target_source", exact = TRUE),
    attr(target, "source", exact = TRUE)
  )
  values <- values[!vapply(values, is.null, logical(1))]
  tolower(paste(unlist(values, use.names = FALSE), collapse = " "))
}

ccc_assert_external_target <- function(target) {
  text <- ccc_forbidden_target_text(target)
  blocked <- c(
    "mock",
    "sham",
    "fallback",
    "self",
    "diagnostic_self_target",
    "raw posterior",
    "raw_posterior",
    "stacked posterior",
    "stacked_posterior",
    "proportional"
  )
  bad <- blocked[vapply(blocked, function(term) grepl(term, text, fixed = TRUE), logical(1))]
  if (length(bad) > 0L) {
    pv_abort(sprintf("CCC target provenance is forbidden for confirmatory calibration: %s.", paste(unique(bad), collapse = ", ")))
  }
  if (!grepl("external", text, fixed = TRUE)) {
    pv_abort("CCC calibration requires an explicit external target provenance.")
  }
  if (!identical(target$target_source, "external_brr_fay_rubin")) {
    pv_abort("CCC calibration requires `target_source = \"external_brr_fay_rubin\"`.")
  }
  invisible(TRUE)
}

ccc_resolve_target <- function(target, fe_names, raw_fe_cov = NULL) {
  if (!is.list(target)) {
    pv_abort("`target` must be a `pvstackr_brr_target` object or a compatible named list.")
  }
  pv_assert_named_list(target, "target")
  required <- c("beta", "T_MI", "fe_names", "target_source")
  missing <- setdiff(required, names(target))
  if (length(missing) > 0L) {
    pv_abort(sprintf("CCC target is missing required field(s): %s.", paste(missing, collapse = ", ")))
  }
  ccc_assert_external_target(target)
  if (!identical(target$fe_names, fe_names)) {
    pv_abort("CCC target `fe_names` must exactly match fixed-effect draw names in order.")
  }

  beta <- target$beta
  if (!is.numeric(beta) || any(!is.finite(beta)) || is.null(names(beta))) {
    pv_abort("CCC target `beta` must be a finite named numeric vector.")
  }
  missing_beta <- setdiff(fe_names, names(beta))
  if (length(missing_beta) > 0L) {
    pv_abort(sprintf("CCC target `beta` is missing fixed-effect name(s): %s.", paste(missing_beta, collapse = ", ")))
  }
  beta <- beta[fe_names]

  T_MI <- target$T_MI
  pv_validate_target_matrix(T_MI, fe_names, raw_fe_cov = raw_fe_cov)
  T_MI <- T_MI[fe_names, fe_names, drop = FALSE]

  list(
    beta = beta,
    T_MI = T_MI,
    target_source = target$target_source,
    target_hash = target$target_hash %||% NA_character_
  )
}

ccc_resolve_param_map <- function(draws, target, fe_names = NULL, param_map = NULL) {
  all_names <- colnames(draws)
  if (!is.null(param_map)) {
    if (!is.list(param_map)) {
      pv_abort("`param_map` must be a list when supplied.")
    }
    fe_ref <- param_map$fe_idx %||% param_map$fe_names
    if (is.null(fe_ref)) {
      pv_abort("`param_map` must contain `fe_idx` or `fe_names`.")
    }
  } else {
    fe_ref <- fe_names %||% target$fe_names %||% names(target$beta)
  }
  if (is.null(fe_ref)) {
    pv_abort("Fixed-effect names must be supplied or available from `target`.")
  }

  if (is.character(fe_ref)) {
    fe_names_out <- pv_validate_unique_character(fe_ref, "fe_names")
    fe_idx <- match(fe_names_out, all_names)
    if (anyNA(fe_idx)) {
      pv_abort(sprintf("Fixed-effect draw column(s) not found: %s.", paste(fe_names_out[is.na(fe_idx)], collapse = ", ")))
    }
  } else {
    fe_idx <- as.integer(fe_ref)
    if (length(fe_idx) < 1L || anyNA(fe_idx) || any(fe_idx < 1L) || any(fe_idx > length(all_names)) || anyDuplicated(fe_idx)) {
      pv_abort("`fe_idx` must identify unique draw columns.")
    }
    fe_names_out <- all_names[fe_idx]
  }

  vc_ref <- if (!is.null(param_map)) param_map$vc_idx %||% param_map$vc_names else NULL
  if (is.null(vc_ref)) {
    vc_idx <- setdiff(seq_along(all_names), fe_idx)
  } else if (length(vc_ref) == 0L) {
    vc_idx <- integer()
  } else if (is.character(vc_ref)) {
    vc_idx <- match(vc_ref, all_names)
    if (anyNA(vc_idx)) {
      pv_abort(sprintf("Variance-component draw column(s) not found: %s.", paste(vc_ref[is.na(vc_idx)], collapse = ", ")))
    }
  } else {
    vc_idx <- as.integer(vc_ref)
    if (anyNA(vc_idx) || any(vc_idx < 1L) || any(vc_idx > length(all_names)) || anyDuplicated(vc_idx)) {
      pv_abort("`vc_idx` must identify unique draw columns.")
    }
  }
  if (length(intersect(fe_idx, vc_idx)) > 0L) {
    pv_abort("Fixed-effect and variance-component draw columns must not overlap.")
  }

  list(
    fe_idx = fe_idx,
    vc_idx = vc_idx,
    fe_names = fe_names_out,
    vc_names = all_names[vc_idx]
  )
}

ccc_chol <- function(matrix, label) {
  tryCatch(
    chol(pv_symmetrize(matrix), pivot = FALSE),
    error = function(e) pv_abort(sprintf("`%s` must be positive definite for CCC calibration.", label))
  )
}

ccc_calibration_matrix <- function(Sigma_raw, Sigma_target) {
  U_raw <- ccc_chol(Sigma_raw, "Sigma_raw")
  U_target <- ccc_chol(Sigma_target, "Sigma_target")
  L_raw <- t(U_raw)
  L_target <- t(U_target)
  L_raw_inv <- forwardsolve(L_raw, diag(nrow(L_raw)))
  L_target %*% L_raw_inv
}

ccc_center_thresholds <- function() {
  list(
    warn = 1e-2,
    block = 5e-2
  )
}

ccc_conditioning_thresholds <- function() {
  list(
    kappa_A = list(warn = 1e6, block = 1e8)
  )
}

ccc_conditioning_band <- function(value, thresholds) {
  if (value >= thresholds$block) {
    "red"
  } else if (value >= thresholds$warn) {
    "yellow"
  } else {
    "green"
  }
}

ccc_conditioning_diagnostics <- function(a_matrix_fro_rel, kappa_A) {
  thresholds <- ccc_conditioning_thresholds()
  kappa_band <- ccc_conditioning_band(kappa_A, thresholds$kappa_A)

  if (identical(kappa_band, "red")) {
    conditioning_status <- "blocked"
    reason_code <- "ccc_conditioning_red"
  } else if (identical(kappa_band, "yellow")) {
    conditioning_status <- "warning"
    reason_code <- "ccc_conditioning_yellow"
  } else {
    conditioning_status <- "ok"
    reason_code <- NA_character_
  }

  list(
    conditioning_status = conditioning_status,
    conditioning_reason_code = reason_code,
    conditioning_band = kappa_band,
    kappa_A_threshold_warn = thresholds$kappa_A$warn,
    kappa_A_threshold_block = thresholds$kappa_A$block,
    conditioning = list(
      gate_metric = "kappa_A",
      gate_value = unname(kappa_A),
      kappa_A = unname(kappa_A),
      a_matrix_fro_rel = unname(a_matrix_fro_rel),
      warn_threshold = thresholds$kappa_A$warn,
      block_threshold = thresholds$kappa_A$block,
      a_matrix_fro_rel_role = "descriptive_not_gated",
      status = conditioning_status,
      band = kappa_band,
      reason_code = reason_code
    )
  )
}

ccc_center_diagnostics <- function(raw_mean, target_mean, Sigma_target) {
  if (!identical(names(raw_mean), names(target_mean))) {
    pv_abort("CCC center diagnostics require aligned named centers.")
  }
  se_target <- sqrt(pmax(diag(Sigma_target), .Machine$double.eps))
  center_shift <- target_mean - raw_mean
  center_shift_se <- center_shift / se_target
  delta_c_by_term <- abs(center_shift_se)
  delta_c_rel <- sqrt(mean(center_shift_se^2))
  delta_c_max <- max(delta_c_by_term)
  thresholds <- ccc_center_thresholds()
  if (delta_c_max >= thresholds$block) {
    center_status <- "blocked"
    band <- "red"
    reason_code <- "center_separation_red"
  } else if (delta_c_max >= thresholds$warn) {
    center_status <- "warning"
    band <- "yellow"
    reason_code <- "center_separation_yellow"
  } else {
    center_status <- "ok"
    band <- "green"
    reason_code <- NA_character_
  }
  center_separation <- list(
    raw_center = raw_mean,
    target_center = target_mean,
    shift = center_shift,
    target_se = se_target,
    shift_rel = center_shift_se,
    max_abs_rel = unname(delta_c_max),
    gate_metric = "delta_c_max",
    gate_value = unname(delta_c_max),
    warn_threshold = thresholds$warn,
    block_threshold = thresholds$block,
    band = band,
    reason_code = reason_code
  )

  list(
    delta_c_rel = unname(delta_c_rel),
    delta_c_max = unname(delta_c_max),
    delta_c_by_term = delta_c_by_term,
    center_shift = center_shift,
    center_shift_se = center_shift_se,
    center_threshold_warn = thresholds$warn,
    center_threshold_block = thresholds$block,
    center_status = center_status,
    center_reason_code = reason_code,
    center_separation = center_separation
  )
}

ccc_center_warnings <- function(center_diagnostics) {
  status <- center_diagnostics$center_status
  if (identical(status, "blocked")) {
    return(sprintf(
      "Center separation diagnostic delta_c_max = %.6g (RMS delta_c_rel = %.6g) exceeds the block threshold %.6g.",
      center_diagnostics$delta_c_max,
      center_diagnostics$delta_c_rel,
      center_diagnostics$center_threshold_block
    ))
  }
  if (identical(status, "warning")) {
    return(sprintf(
      "Center separation diagnostic delta_c_max = %.6g (RMS delta_c_rel = %.6g) exceeds the warning threshold %.6g.",
      center_diagnostics$delta_c_max,
      center_diagnostics$delta_c_rel,
      center_diagnostics$center_threshold_warn
    ))
  }
  character()
}

ccc_conditioning_warnings <- function(diagnostics) {
  status <- diagnostics$conditioning_status
  if (identical(status, "blocked")) {
    return(sprintf(
      "CCC conditioning diagnostic kappa_A = %.6g exceeds the block threshold %.6g; reportable stack_direct estimates are blocked.",
      diagnostics$kappa_A,
      diagnostics$kappa_A_threshold_block
    ))
  }
  if (identical(status, "warning")) {
    return(sprintf(
      "CCC conditioning diagnostic kappa_A = %.6g exceeds the warning threshold %.6g.",
      diagnostics$kappa_A,
      diagnostics$kappa_A_threshold_warn
    ))
  }
  character()
}

ccc_diagnostics <- function(A, Sigma_raw, Sigma_target, Sigma_cal_emp_raw, raw_mean, target_mean) {
  matrix_residual <- pv_symmetrize(A %*% Sigma_raw %*% t(A) - Sigma_target)
  empirical_residual <- pv_symmetrize(Sigma_cal_emp_raw - Sigma_target)
  fro_target <- sqrt(sum(Sigma_target^2))
  sv <- svd(A)$d
  kappa_A <- if (min(sv) > 0) max(sv) / min(sv) else Inf
  a_matrix_fro_rel <- sqrt(sum((A - diag(nrow(A)))^2)) / sqrt(nrow(A))
  center <- ccc_center_diagnostics(raw_mean, target_mean, Sigma_target)
  conditioning <- ccc_conditioning_diagnostics(
    a_matrix_fro_rel = a_matrix_fro_rel,
    kappa_A = kappa_A
  )
  c(
    list(
      a_matrix_fro_rel = a_matrix_fro_rel
    ),
    center,
    conditioning,
    list(
    rho1 = if (fro_target > 0) sqrt(sum(matrix_residual^2)) / fro_target else sqrt(sum(matrix_residual^2)),
    rho2 = max(abs(diag(matrix_residual)) / pmax(abs(diag(Sigma_target)), .Machine$double.eps)),
    matrix_residual = matrix_residual,
    empirical_residual = empirical_residual,
    empirical_fro_rel = if (fro_target > 0) sqrt(sum(empirical_residual^2)) / fro_target else sqrt(sum(empirical_residual^2)),
    kappa_A = kappa_A
    )
  )
}

ccc_block_matrix <- function(A, all_names, fe_idx) {
  A_full <- diag(length(all_names))
  dimnames(A_full) <- list(all_names, all_names)
  A_full[fe_idx, fe_idx] <- A
  A_full
}

ccc_calibrate <- function(
  draws,
  target,
  fe_names = NULL,
  param_map = NULL,
  center = c("target", "posterior")
) {
  draws <- ccc_as_draw_matrix(draws)
  center <- match.arg(center)
  map <- ccc_resolve_param_map(draws, target, fe_names = fe_names, param_map = param_map)

  draws_fe <- draws[, map$fe_idx, drop = FALSE]
  Sigma_raw <- pv_symmetrize(stats::cov(draws_fe))
  resolved_target <- ccc_resolve_target(target, map$fe_names, raw_fe_cov = Sigma_raw)
  Sigma_target <- resolved_target$T_MI
  A <- ccc_calibration_matrix(Sigma_raw, Sigma_target)
  A_full <- ccc_block_matrix(A, colnames(draws), map$fe_idx)

  raw_mean <- colMeans(draws_fe)
  output_mean <- if (identical(center, "target")) resolved_target$beta else raw_mean
  centered <- sweep(draws_fe, 2L, raw_mean, FUN = "-")
  draws_fe_cal <- sweep(centered %*% t(A), 2L, output_mean, FUN = "+")
  colnames(draws_fe_cal) <- map$fe_names

  draws_calibrated <- draws
  draws_calibrated[, map$fe_idx] <- draws_fe_cal
  Sigma_cal_emp_raw <- pv_symmetrize(stats::cov(draws_fe_cal))
  diagnostics <- ccc_diagnostics(
    A,
    Sigma_raw,
    Sigma_target,
    Sigma_cal_emp_raw,
    raw_mean = raw_mean,
    target_mean = resolved_target$beta
  )
  warnings <- c(ccc_center_warnings(diagnostics), ccc_conditioning_warnings(diagnostics))

  out <- list(
    draws_calibrated = draws_calibrated,
    draws_fe_cal = draws_fe_cal,
    A = A,
    A_full = A_full,
    psi_hat = output_mean,
    psi_raw = raw_mean,
    psi_target = resolved_target$beta,
    Sigma_raw = Sigma_raw,
    Sigma_target = Sigma_target,
    Sigma_cal_emp = Sigma_cal_emp_raw,
    Sigma_cal_emp_raw = Sigma_cal_emp_raw,
    diagnostics = diagnostics,
    flags = list(
      nearpd_raw = FALSE,
      nearpd_target = FALSE,
      nearpd_cal = FALSE,
      target_repaired = FALSE,
      center_separation_warning = identical(diagnostics$center_status, "warning"),
      center_separation_blocked = identical(diagnostics$center_status, "blocked"),
      conditioning_warning = identical(diagnostics$conditioning_status, "warning"),
      conditioning_blocked = identical(diagnostics$conditioning_status, "blocked"),
      kappa_a_warning = identical(diagnostics$conditioning_status, "warning"),
      kappa_a_blocked = identical(diagnostics$conditioning_status, "blocked")
    ),
    param_map = map,
    control = list(
      center = center,
      allow_target_nearpd = FALSE
    ),
    target_source = resolved_target$target_source,
    target_hash = resolved_target$target_hash,
    center = center,
    policy = list(
      fixed_effects_only = TRUE,
      target_repair = "forbidden",
      vc_passthrough = length(map$vc_idx) > 0L
    ),
    ccc_status = "ok",
    schema_version = "0.1.0",
    provenance = list(
      function_name = "ccc_calibrate",
      package = "pvstackr"
    ),
    warnings = warnings
  )
  class(out) <- c("pvstackr_ccc", "list")
  out
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
