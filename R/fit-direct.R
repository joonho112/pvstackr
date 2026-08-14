pv_fit_direct_estimates <- function(ccc, target, conf_level) {
  fe_names <- target$fe_names
  estimate <- ccc$psi_hat[fe_names]
  se <- sqrt(diag(ccc$Sigma_target))[fe_names]
  df <- target$df[fe_names]
  if (!is.numeric(estimate) || any(!is.finite(estimate))) {
    pv_abort("Direct-fit calibrated estimates must be finite.")
  }
  if (!is.numeric(se) || any(!is.finite(se)) || any(se <= 0)) {
    pv_abort("Direct-fit target standard errors must be finite and positive.")
  }
  if (!is.numeric(df) || any(is.na(df)) || any(df <= 0)) {
    pv_abort("Direct-fit target degrees of freedom must be positive.")
  }
  df_complete <- target$df_complete
  if (is.null(df_complete)) {
    df_complete <- stats::setNames(rep(NA_real_, length(fe_names)), fe_names)
  } else {
    df_complete <- df_complete[fe_names]
  }
  alpha <- 1 - conf_level
  crit <- stats::qt(1 - alpha / 2, df = df)
  data.frame(
    term = fe_names,
    estimate = unname(estimate),
    se = unname(se),
    std.error = unname(se),
    df = unname(df),
    df_method = target$df_method,
    df_complete = unname(df_complete),
    conf_level = conf_level,
    conf_low = unname(estimate - crit * se),
    conf_high = unname(estimate + crit * se),
    conf.low = unname(estimate - crit * se),
    conf.high = unname(estimate + crit * se),
    interval_role = target$interval_role,
    coverage_claim_allowed = target$coverage_claim_allowed,
    parameter_scope = "fixed_effect",
    target_source = target$target_source,
    target_hash = target$target_hash,
    stringsAsFactors = FALSE
  )
}

pv_fit_direct_status <- function(stack_fit, ccc, sampler_gate = NULL) {
  sampler_gate <- sampler_gate %||% list(
    status = "ok",
    reason_codes = character(),
    warnings = character()
  )
  status <- sampler_gate$status
  reason_codes <- sampler_gate$reason_codes
  warnings <- sampler_gate$warnings

  promote <- function(current, candidate) {
    ranks <- c(ok = 0L, warning = 1L, blocked = 2L)
    if (ranks[[candidate]] > ranks[[current]]) candidate else current
  }

  prior_policy <- stack_fit$meta$prior_policy %||% list()
  if (isTRUE(prior_policy$explicit_prior_warning)) {
    prior_reason_code <- prior_policy$reason_code
    if (!is.character(prior_reason_code) ||
        length(prior_reason_code) != 1L ||
        is.na(prior_reason_code) ||
        !nzchar(prior_reason_code)) {
      prior_reason_code <- "explicit_prior_warning"
    }
    reason_codes <- c(reason_codes, prior_reason_code)
    warnings <- c(warnings, prior_policy$warning)
    status <- promote(status, "warning")
  }

  center_status <- ccc$diagnostics$center_status %||% "ok"
  if (identical(center_status, "warning")) {
    reason_codes <- c(reason_codes, ccc$diagnostics$center_reason_code %||% "center_separation_warning")
    warnings <- c(warnings, ccc$warnings)
    status <- promote(status, "warning")
  }
  if (identical(center_status, "blocked")) {
    reason_codes <- c(reason_codes, ccc$diagnostics$center_reason_code %||% "center_separation_blocked")
    warnings <- c(warnings, ccc$warnings)
    status <- promote(status, "blocked")
  }

  conditioning_status <- ccc$diagnostics$conditioning_status %||% "ok"
  if (identical(conditioning_status, "warning")) {
    reason_codes <- c(reason_codes, ccc$diagnostics$conditioning_reason_code %||% "ccc_conditioning_warning")
    warnings <- c(warnings, ccc$warnings)
    status <- promote(status, "warning")
  }
  if (identical(conditioning_status, "blocked")) {
    reason_codes <- c(reason_codes, ccc$diagnostics$conditioning_reason_code %||% "ccc_conditioning_blocked")
    warnings <- c(warnings, ccc$warnings)
    status <- promote(status, "blocked")
  }

  reason_codes <- unique(reason_codes[nzchar(reason_codes)])
  warnings <- unique(warnings[nzchar(warnings)])
  if (identical(status, "warning") && length(warnings) == 0L) {
    warnings <- "Direct fit produced warning-level diagnostics."
  }
  list(
    status = status,
    reason_codes = reason_codes,
    warnings = warnings
  )
}

pv_fit_direct_blocked_preflight <- function(preflight) {
  keep <- c(
    "formula_string", "rhs_string", "fe_names", "target_hash",
    "target_source", "policy", "binding_proof"
  )
  missing <- setdiff(keep, names(preflight))
  if (length(missing) > 0L) {
    pv_abort(sprintf(
      "Stack-direct preflight is missing blocked-snapshot field(s): %s.",
      paste(missing, collapse = ", ")
    ))
  }
  pv_binding_proof_validate(preflight$binding_proof)
  out <- list(
    formula_string = unname(as.character(preflight$formula_string)),
    rhs_string = unname(as.character(preflight$rhs_string)),
    fe_names = unname(as.character(preflight$fe_names)),
    target_hash = unname(as.character(preflight$target_hash)),
    target_source = unname(as.character(preflight$target_source)),
    policy = list(
      fixed_effects_only = isTRUE(preflight$policy$fixed_effects_only),
      target_repair = unname(as.character(preflight$policy$target_repair))
    ),
    binding_proof = preflight$binding_proof
  )
  class(out) <- c("pvstackr_stack_direct_preflight_snapshot", "list")
  out
}

pv_fit_direct_retained_preflight <- function(preflight, keep_data) {
  keep_data <- pv_assert_scalar_logical(keep_data, "keep_data")
  out <- preflight
  out$formula <- tryCatch(
    stats::as.formula(out$formula_string, env = baseenv()),
    error = function(error) NULL
  )
  if (is.null(out$formula)) {
    pv_abort("A redacted stack_direct preflight requires a canonical formula string.")
  }
  class(out) <- c("pvstackr_stack_direct_preflight", "list")
  out
}

pv_fit_direct_plain_named_numeric <- function(x) {
  out <- as.numeric(x)
  names(out) <- unname(as.character(names(x)))
  out
}

pv_fit_direct_plain_matrix <- function(x) {
  rn <- rownames(x)
  cn <- colnames(x)
  dimnames <- if (is.null(rn) && is.null(cn)) {
    NULL
  } else {
    list(
      if (is.null(rn)) NULL else unname(as.character(rn)),
      if (is.null(cn)) NULL else unname(as.character(cn))
    )
  }
  matrix(
    as.numeric(x),
    nrow = nrow(x),
    ncol = ncol(x),
    dimnames = dimnames
  )
}

pv_fit_direct_numeric_equal <- function(
  expected,
  observed,
  allow_na = FALSE,
  allow_positive_infinity = FALSE
) {
  if (!is.numeric(expected) || !is.numeric(observed) ||
      !identical(typeof(expected), typeof(observed)) ||
      !identical(attributes(expected), attributes(observed)) ||
      length(expected) != length(observed)) {
    return(FALSE)
  }
  expected_nan <- is.nan(expected)
  observed_nan <- is.nan(observed)
  if (any(expected_nan) || any(observed_nan)) {
    return(FALSE)
  }
  expected_na <- is.na(expected)
  observed_na <- is.na(observed)
  if (!identical(expected_na, observed_na) ||
      (any(expected_na) && !allow_na)) {
    return(FALSE)
  }
  expected_infinite <- is.infinite(expected)
  observed_infinite <- is.infinite(observed)
  if (!identical(expected_infinite, observed_infinite)) {
    return(FALSE)
  }
  if (any(expected_infinite) &&
      (!allow_positive_infinity ||
       any(expected[expected_infinite] < 0) ||
       !identical(
         expected[expected_infinite],
         observed[observed_infinite]
       ))) {
    return(FALSE)
  }
  finite <- !(expected_na | expected_infinite)
  difference <- abs(expected[finite] - observed[finite])
  scale <- pmax(abs(expected[finite]), abs(observed[finite]))
  all(difference <= 1e-12 + 1e-10 * scale)
}

pv_fit_direct_require_numeric_equal <- function(
  expected,
  observed,
  label,
  allow_na = FALSE,
  allow_positive_infinity = FALSE
) {
  if (!pv_fit_direct_numeric_equal(
    expected,
    observed,
    allow_na = allow_na,
    allow_positive_infinity = allow_positive_infinity
  )) {
    pv_abort(sprintf(
      "Blocked target `%s` must match its canonical per-PV reconstruction.",
      label
    ))
  }
  invisible(TRUE)
}

pv_fit_direct_independent_target_item <- function(
  item,
  declared_pv_col,
  declared_rep_weight_cols,
  declared_R,
  declared_fay_k,
  declared_multiplier
) {
  declared_pv_col <- pv_assert_scalar_string(
    declared_pv_col,
    "blocked target declared pv_col"
  )
  if (!identical(item$pv_col, declared_pv_col)) {
    pv_abort("Blocked target per-PV `pv_col` must exactly align with declared `pv_cols`.")
  }
  beta <- pv_fit_direct_plain_named_numeric(item$beta)
  if (!identical(unname(as.character(item$fe_names)), names(beta))) {
    pv_abort("Blocked target per-PV fixed-effect names must align with beta.")
  }
  item_R <- pv_assert_scalar_number(
    item$R,
    "blocked target per-PV R",
    integer = TRUE,
    lower = 2
  )
  if (!identical(as.integer(item_R), declared_R)) {
    pv_abort("Blocked target per-PV replicate count must exactly match top-level R.")
  }
  pv_fit_direct_require_numeric_equal(
    declared_fay_k,
    item$fay_k,
    "per_pv$fay_k"
  )
  pv_fit_direct_require_numeric_equal(
    declared_multiplier,
    item$fay_variance_multiplier,
    "per_pv$fay_variance_multiplier"
  )
  if (!is.matrix(item$replicate_diff) ||
      !is.numeric(item$replicate_diff) ||
      any(!is.finite(item$replicate_diff))) {
    pv_abort("Blocked target per-PV replicate differences must be a finite numeric matrix.")
  }
  replicate_diff <- pv_fit_direct_plain_matrix(item$replicate_diff)
  declared_rep_weight_cols <- unname(as.character(declared_rep_weight_cols))
  if (!identical(dim(replicate_diff), c(length(beta), declared_R)) ||
      !identical(rownames(replicate_diff), names(beta)) ||
      !identical(colnames(replicate_diff), declared_rep_weight_cols)) {
    pv_abort("Blocked target per-PV replicate differences must exactly align with fixed effects and declared replicate weights.")
  }
  U <- declared_multiplier * tcrossprod(replicate_diff)
  dimnames(U) <- list(names(beta), names(beta))
  pv_fit_direct_require_numeric_equal(U, item$U, "per_pv$U")
  replicate_beta <- sweep(replicate_diff, 1L, beta, FUN = "+")
  pv_fit_direct_require_numeric_equal(
    replicate_beta,
    item$replicate_beta,
    "per_pv$replicate_beta"
  )
  list(
    beta = beta,
    U = U,
    fe_names = unname(as.character(item$fe_names)),
    pv_col = unname(as.character(declared_pv_col)),
    R = declared_R,
    fay_k = declared_fay_k,
    fay_variance_multiplier = declared_multiplier,
    replicate_beta = replicate_beta,
    replicate_diff = replicate_diff
  )
}

pv_fit_direct_independent_target_provenance <- function(provenance) {
  function_name <- pv_assert_scalar_string(
    provenance$function_name,
    "blocked target provenance function_name"
  )
  assembled_at <- pv_assert_scalar_string(
    provenance$assembled_at,
    "blocked target provenance assembled_at"
  )
  package <- pv_assert_scalar_string(
    provenance$package,
    "blocked target provenance package"
  )
  timestamp_pattern <-
    "^[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}(\\.[0-9]{1,6})?$"
  timestamp_base <- substr(assembled_at, 1L, 19L)
  parsed_time <- suppressWarnings(as.POSIXct(
    timestamp_base,
    format = "%Y-%m-%d %H:%M:%S",
    tz = "UTC"
  ))
  timestamp_roundtrip <- if (is.na(parsed_time)) {
    NA_character_
  } else {
    format(parsed_time, "%Y-%m-%d %H:%M:%S", tz = "UTC")
  }
  if (!identical(function_name, "pv_brr_target") ||
      !identical(package, "pvstackr") ||
      !grepl(timestamp_pattern, assembled_at) ||
      !identical(timestamp_roundtrip, timestamp_base)) {
    pv_abort("Blocked target provenance must use the canonical pv_brr_target scalar identity and timestamp.")
  }
  list(
    function_name = function_name,
    assembled_at = assembled_at,
    package = package
  )
}

pv_fit_direct_independent_target_v01 <- function(target) {
  formula_string <- unname(as.character(target$formula_string))
  pv_cols <- unname(as.character(target$pv_cols))
  if (length(pv_cols) != target$M ||
      anyNA(pv_cols) || any(!nzchar(pv_cols)) || anyDuplicated(pv_cols)) {
    pv_abort("Blocked target `pv_cols` must be unique, non-missing, and aligned with `M`.")
  }
  if (!identical(target$engine, "lm") ||
      !identical(
        target$policy$replicate_weight_role,
        "external_design_variance_only"
      ) ||
      !identical(target$schema_version, "0.1.0")) {
    pv_abort("Blocked target must use the current canonical BRR-Fay engine, policy, and schema identity.")
  }
  M <- as.integer(target$M)
  R <- as.integer(target$R)
  fay_k <- as.numeric(target$fay_k)
  rep_weight_cols <- unname(as.character(target$rep_weight_cols))
  if (length(rep_weight_cols) != R || anyNA(rep_weight_cols) ||
      any(!nzchar(rep_weight_cols)) || anyDuplicated(rep_weight_cols)) {
    pv_abort("Blocked target top-level R must equal the unique declared replicate-weight count.")
  }
  fay_variance_multiplier <- 1 / (R * (1 - fay_k)^2)
  pv_fit_direct_require_numeric_equal(
    fay_variance_multiplier,
    target$fay_variance_multiplier,
    "fay_variance_multiplier"
  )
  per_pv <- lapply(
    seq_along(target$per_pv),
    function(index) pv_fit_direct_independent_target_item(
      target$per_pv[[index]],
      pv_cols[[index]],
      rep_weight_cols,
      R,
      fay_k,
      fay_variance_multiplier
    )
  )
  beta_rows <- do.call(rbind, lapply(per_pv, `[[`, "beta"))
  colnames(beta_rows) <- target$fe_names
  pooled <- tryCatch(
    rubin_pool_matrix(
      beta = beta_rows,
      U = lapply(per_pv, `[[`, "U"),
      orientation = "rows_pv",
      conf_level = 0.95,
      allow_m1 = M == 1L,
      df_method = target$df_method,
      df_complete = target$df_complete
    ),
    error = function(error) pv_abort(
      "Blocked target canonical Rubin recomputation failed."
    )
  )
  expected_target_fields <- list(
    beta = pooled$beta,
    beta_bar = pooled$beta,
    U_bar = pooled$U_bar,
    B = pooled$B,
    T_MI = pooled$T_MI,
    total_var = pooled$total_var,
    se = pooled$se,
    df = pooled$df,
    df_classic = pooled$df_classic,
    df_complete = pooled$df_complete,
    lambda = pooled$lambda,
    fmi = pooled$fmi,
    riv = pooled$riv
  )
  for (field in names(expected_target_fields)) {
    pv_fit_direct_require_numeric_equal(
      expected_target_fields[[field]],
      target[[field]],
      field,
      allow_na = identical(field, "df_complete"),
      allow_positive_infinity = field %in% c(
        "df", "df_classic", "df_complete"
      )
    )
  }
  out <- list(
    beta = pooled$beta,
    beta_bar = pooled$beta,
    U_bar = pooled$U_bar,
    B = pooled$B,
    T_MI = pooled$T_MI,
    total_var = pooled$total_var,
    se = pooled$se,
    df = pooled$df,
    df_classic = pooled$df_classic,
    df_method = unname(as.character(target$df_method)),
    df_complete = pooled$df_complete,
    interval_role = unname(as.character(target$interval_role)),
    coverage_claim_allowed = isTRUE(target$coverage_claim_allowed),
    lambda = pooled$lambda,
    fmi = pooled$fmi,
    riv = pooled$riv,
    fe_names = unname(as.character(target$fe_names)),
    per_pv = per_pv,
    formula = stats::as.formula(formula_string, env = baseenv()),
    formula_string = formula_string,
    rhs_string = unname(as.character(target$rhs_string)),
    M = M,
    R = R,
    fay_k = fay_k,
    fay_variance_multiplier = fay_variance_multiplier,
    pv_cols = pv_cols,
    weight_col = unname(as.character(target$weight_col)),
    rep_weight_cols = rep_weight_cols,
    id_cols = unname(as.character(target$id_cols)),
    design_hash = unname(as.character(target$design_hash)),
    target_source = unname(as.character(target$target_source)),
    target_hash = unname(as.character(target$target_hash)),
    engine = unname(as.character(target$engine)),
    policy = list(
      replicate_weight_role = "external_design_variance_only",
      target_repair = unname(as.character(target$policy$target_repair)),
      fixed_effects_only = isTRUE(target$policy$fixed_effects_only),
      df_method = unname(as.character(target$policy$df_method)),
      interval_role = unname(as.character(target$policy$interval_role)),
      coverage_claim_allowed = isTRUE(target$policy$coverage_claim_allowed)
    ),
    schema_version = "0.1.0",
    provenance = pv_fit_direct_independent_target_provenance(
      target$provenance
    ),
    warnings = character()
  )
  class(out) <- c("pvstackr_brr_target", "list")
  validate_pvstackr_brr_target(out)
  out
}

pv_fit_direct_independent_target_v02 <- function(target) {
  validate_pvstackr_brr_target(target)
  out <- target
  out$formula <- stats::as.formula(
    unname(as.character(target$formula_string)),
    env = baseenv()
  )
  class(out) <- c("pvstackr_brr_target", "list")
  validate_pvstackr_brr_target(out)
  out
}

pv_fit_direct_independent_target <- function(target) {
  if (identical(target$schema_version, "0.2.0")) {
    return(pv_fit_direct_independent_target_v02(target))
  }
  pv_fit_direct_independent_target_v01(target)
}

pv_fit_direct_blocked_control <- function(control) {
  pv_fit_blocked_control(control)
}

pv_fit_direct_blocked_ccc_diagnostics <- function(stack_fit, ccc) {
  validate_pvstackr_stack_fit(stack_fit)
  validate_pvstackr_ccc(ccc)
  diagnostics <- ccc$diagnostics
  prior_policy <- stack_fit$meta$prior_policy %||% pv_stack_prior_policy(NULL)
  list(
    source = "ccc_reportability_gate",
    target_hash = unname(as.character(ccc$target_hash)),
    prior = list(
      explicit_warning = isTRUE(prior_policy$explicit_prior_warning),
      reason_code = if (isTRUE(prior_policy$explicit_prior_warning)) {
        "explicit_prior_warning"
      } else {
        NA_character_
      }
    ),
    center = list(
      status = unname(as.character(diagnostics$center_status)),
      reason_code = unname(as.character(diagnostics$center_reason_code)),
      delta_c_rel = unname(as.numeric(diagnostics$delta_c_rel)),
      delta_c_max = unname(as.numeric(diagnostics$delta_c_max)),
      warn_threshold = unname(as.numeric(diagnostics$center_threshold_warn)),
      block_threshold = unname(as.numeric(diagnostics$center_threshold_block))
    ),
    conditioning = list(
      status = unname(as.character(diagnostics$conditioning_status)),
      reason_code = unname(as.character(diagnostics$conditioning_reason_code)),
      band = unname(as.character(diagnostics$conditioning_band)),
      kappa_A = unname(as.numeric(diagnostics$kappa_A)),
      a_matrix_fro_rel = unname(as.numeric(diagnostics$a_matrix_fro_rel)),
      warn_threshold = unname(as.numeric(diagnostics$kappa_A_threshold_warn)),
      block_threshold = unname(as.numeric(diagnostics$kappa_A_threshold_block))
    ),
    residual = list(
      rho1 = unname(as.numeric(diagnostics$rho1)),
      rho2 = unname(as.numeric(diagnostics$rho2)),
      empirical_fro_rel = unname(as.numeric(diagnostics$empirical_fro_rel))
    )
  )
}

pv_fit_direct_slim_blocked_status <- function(sampler_gate, ccc_diagnostics) {
  status <- sampler_gate$status
  reason_codes <- sampler_gate$reason_codes
  warnings <- sampler_gate$warnings
  promote <- function(current, candidate) {
    ranks <- c(ok = 0L, warning = 1L, blocked = 2L)
    if (ranks[[candidate]] > ranks[[current]]) candidate else current
  }

  if (isTRUE(ccc_diagnostics$prior$explicit_warning)) {
    reason_codes <- c(reason_codes, ccc_diagnostics$prior$reason_code)
    warnings <- c(warnings, pv_stack_prior_warnings(TRUE))
    status <- promote(status, "warning")
  }
  diagnostic_bridge <- list(
    center_status = ccc_diagnostics$center$status,
    center_reason_code = ccc_diagnostics$center$reason_code,
    delta_c_rel = ccc_diagnostics$center$delta_c_rel,
    delta_c_max = ccc_diagnostics$center$delta_c_max,
    center_threshold_warn = ccc_diagnostics$center$warn_threshold,
    center_threshold_block = ccc_diagnostics$center$block_threshold,
    conditioning_status = ccc_diagnostics$conditioning$status,
    conditioning_reason_code = ccc_diagnostics$conditioning$reason_code,
    kappa_A = ccc_diagnostics$conditioning$kappa_A,
    kappa_A_threshold_warn = ccc_diagnostics$conditioning$warn_threshold,
    kappa_A_threshold_block = ccc_diagnostics$conditioning$block_threshold
  )
  if (ccc_diagnostics$center$status %in% c("warning", "blocked")) {
    reason_codes <- c(reason_codes, ccc_diagnostics$center$reason_code)
    warnings <- c(warnings, ccc_center_warnings(diagnostic_bridge))
    status <- promote(status, ccc_diagnostics$center$status)
  }
  if (ccc_diagnostics$conditioning$status %in% c("warning", "blocked")) {
    reason_codes <- c(reason_codes, ccc_diagnostics$conditioning$reason_code)
    warnings <- c(warnings, ccc_conditioning_warnings(diagnostic_bridge))
    status <- promote(status, ccc_diagnostics$conditioning$status)
  }
  list(
    status = status,
    reason_codes = unique(reason_codes[nzchar(reason_codes)]),
    warnings = unique(warnings[nzchar(warnings)])
  )
}

pv_fit_direct_blocked_redaction <- function(
  source = c("ccc_reportability_gate", "sampler_gate")
) {
  source <- match.arg(source)
  list(
    status = "withheld",
    policy = "generic_blocked_fail_closed",
    source = source,
    withheld = c(
      "design", "stack_fit", "ccc", "backend_fit", "prepared_data",
      "raw_draws", "calibrated_draws", "nuisance_draws", "log_lik",
      "A", "A_full", "beta", "psi", "Sigma", "target_moment_copies"
    )
  )
}

pv_fit_not_implemented <- function(method) {
  pv_abort(sprintf(
    "`method = \"%s\"` is recognized but not implemented in this package stage. Use `method = \"stack_direct\"`, `method = \"per_pv\"`, `method = \"stack_psis\"`, `pv_fit_direct()`, `pv_fit_reference()`, or `pv_fit_stack_psis()`.",
    method
  ))
}

#' Fit a pvstackr Method
#'
#' `pv_fit()` is the generic public fitting entry point. In this package stage,
#' `method = "stack_direct"`, `method = "per_pv"`, and
#' `method = "stack_psis"` are implemented.
#'
#' @details
#' ## Method dispatch
#' `pv_fit()` validates `method` and `control`, then forwards to a
#' method-specific fitter:
#' - `"stack_direct"` (default) dispatches to [pv_fit_direct()] and **requires**
#'   a `pvstackr_brr_target` `target` from [pv_brr_target()]; `pv_fit()` errors
#'   if `target` is `NULL`.
#' - `"per_pv"` dispatches to [pv_fit_reference()] (the orthodox per-PV
#'   reference); `target` is ignored.
#' - `"stack_psis"` dispatches to [pv_fit_stack_psis()] (the PSIS-reweighted
#'   stacked path); `target` is ignored.
#'
#' If `control` is `NULL` it is constructed with `pv_control(method = method)`;
#' otherwise `control$method` **must** equal `method` or the call errors. Any
#' arguments in `...` are forwarded verbatim to the dispatched fitter, so consult
#' that fitter's signature for the available extras (for example, the injected
#' backend adapters of [pv_fit_direct()], or caller-declared external PSIS weights of
#' [pv_fit_stack_psis()]).
#'
#' In this package stage, only `stack_direct` output is coverage-claimable: its
#' intervals are backed by the external Rubin/BRR-Fay target. `per_pv` and
#' `stack_psis` intervals are descriptive/reference, even with Barnard-Rubin
#' degrees of freedom.
#'
#' @param data Analysis data frame.
#' @param formula Two-sided formula with `OUTCOME` on the left-hand side.
#' @param target A method-specific target object. For `method = "stack_direct"`,
#'   this must be a `pvstackr_brr_target` from [pv_brr_target()] and is required.
#'   Ignored for `"per_pv"` and `"stack_psis"`. Default `NULL`.
#' @param method Public method identifier. Character scalar; one of
#'   `"stack_direct"`, `"stack_psis"`, or `"per_pv"`. Default `"stack_direct"`.
#' @param control Optional [pv_control()] object. If `NULL`, one is built with
#'   `pv_control(method = method)`. If supplied, `control$method` must equal
#'   `method`. Default `NULL`.
#' @param ... Additional arguments forwarded to the dispatched fitter:
#'   [pv_fit_direct()] for `"stack_direct"`, [pv_fit_reference()] for `"per_pv"`,
#'   or [pv_fit_stack_psis()] for `"stack_psis"`.
#'
#' @returns A `pvstackr_fit` object. Read it with the accessors
#'   ([get_estimates()], [get_target()], [get_draws()], [get_diagnostics()])
#'   rather than by `$`-indexing; user-facing components include `estimates`,
#'   `draws`, `target`, `diagnostics`, `status`, `reason_codes`, `warnings`, and
#'   `method`.
#' @examples
#' pisa_tiny <- read.csv(
#'   system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
#' )
#'
#' # Build the design and the external BRR-Fay target on the fixture.
#' design <- pv_design(
#'   pisa_tiny, formula = OUTCOME ~ x + female,
#'   pv_suffix = "READ", expected_M = 2L, expected_R = 4L, id_cols = "CNTSTUID"
#' )
#' target <- pv_brr_target(
#'   pisa_tiny, formula = OUTCOME ~ x + female,
#'   pv_cols = design$pv_cols, weight_col = design$weight_col,
#'   rep_weight_cols = design$rep_weight_cols, fay_k = design$fay_k,
#'   id_cols = design$id_cols
#' )
#'
#' # A live stack_direct fit needs an injected/precomputed backend:
#' \dontrun{
#' fit <- pv_fit(
#'   data = pisa_tiny, formula = OUTCOME ~ x + female,
#'   target = target, method = "stack_direct",
#'   fit_function = my_backend_adapter, draws_function = my_draws_adapter
#' )
#' }
#'
#' # Inspect the reportable object surface via the bundled cached fit instead.
#' path <- system.file(
#'   "extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr"
#' )
#' if (nzchar(path)) {
#'   fit <- readRDS(path)$fit     # a stack_direct pvstackr_fit
#'   fit                          # compact console print
#'   head(get_estimates(fit))     # reportable fixed-effect table
#' }
#' @family pvstackr-fitting
#' @seealso [pv_fit_direct()], [pv_fit_reference()], [pv_fit_stack_psis()],
#'   [pv_control()]; [pv_brr_target()], [pv_compare_methods()], [get_estimates()]
#' @export
pv_fit <- function(
  data,
  formula,
  target = NULL,
  method = "stack_direct",
  control = NULL,
  ...
) {
  method <- pv_validate_method(method)
  control <- if (is.null(control)) pv_control(method = method) else pv_validate_control(control)
  if (!identical(control$method, method)) {
    pv_abort("`control$method` must match `method`.")
  }
  if (identical(method, "per_pv")) {
    return(pv_fit_reference(
      data = data,
      formula = formula,
      control = control,
      ...
    ))
  }
  if (identical(method, "stack_psis")) {
    return(pv_fit_stack_psis(
      data = data,
      formula = formula,
      control = control,
      ...
    ))
  }
  if (!identical(method, "stack_direct")) {
    pv_fit_not_implemented(method)
  }
  if (is.null(target)) {
    pv_abort("`target` is required for `method = \"stack_direct\"`.")
  }
  pv_fit_direct(
    data = data,
    formula = formula,
    target = target,
    control = control,
    ...
  )
}

#' Fit the Direct Stacked Plausible-Value Model
#'
#' `pv_fit_direct()` runs the current `stack_direct` path: fixed-effect-only
#' compatibility preflight, one stacked plausible-value fit, CCC calibration to
#' an external BRR-Fay/Rubin target, and assembly of a reportable
#' `pvstackr_fit` object.
#'
#' @details
#' The current `stack_direct` contract is fixed-effect-only and requires an external
#' `pvstackr_brr_target` from [pv_brr_target()]. The formula RHS and derived
#' fixed-effect names must match the target exactly. Group terms such as
#' `(1 | school)` are rejected until a two-level target engine is implemented.
#'
#' The returned `pvstackr_fit` has `status = "ok"`, `"warning"`, or
#' `"blocked"`. Yellow CCC center separation and explicit priors produce
#' warning-status fits; red CCC center separation blocks reportable estimates.
#' Frozen sampler gates also promote R-hat or ESS findings to warning status and
#' block incomplete diagnostics, hard R-hat/ESS failures, configuration
#' mismatches, or any divergence before CCC calibration. A sampler-blocked fit
#' retains the independently valid external target but removes design data,
#' stack/backend fit, CCC payload, estimates, and draws. Its retained target is
#' an exact canonical allowlisted snapshot whose formula uses the safe base
#' environment and whose warnings are empty. Its preflight record is a
#' formula-free snapshot, so extensions and caller/private environments are not
#' retained.
#' CCC calibration-matrix conditioning is also gated by `kappa_A`: values at or
#' above `1e6` produce warning-status fits, and values at or above `1e8` block
#' reportable estimates. Reportable estimates and retained draws are
#' fixed-effect-only.
#'
#' With the default `control$center = "target"`, the reportable fixed-effect
#' estimate table is target-based: `estimate` is the external Rubin/BRR-Fay
#' coefficient after CCC centering, `se` is `sqrt(diag(target$T_MI))`, and the
#' degrees of freedom and interval metadata are inherited from `target`.
#' Therefore `df_method`, `df_complete`, `interval_role`, and
#' `coverage_claim_allowed` describe the external target policy, not a residual
#' df or interval policy estimated by the stacked backend fit. When the target
#' uses classic Rubin df, the reported df are the target's classic Rubin
#' imputation df; when it uses Barnard-Rubin df, `df_complete` is the explicit
#' complete-data df supplied to [pv_brr_target()]. The stacked fit supplies the
#' draw cloud used for calibration, retained calibrated draws when requested,
#' and center-separation agreement diagnostics; it is not a replacement source
#' for the headline Rubin/BRR-Fay fixed-effect numbers.
#'
#' Reportable `pv_fit_direct()` output requires `control$center = "target"`.
#' The lower-level CCC convention `center = "posterior"` is diagnostic and
#' exploratory: it leaves the output fixed-effect center at the raw stacked
#' posterior mean while pairing that center with the external target covariance.
#' `pv_fit_direct()` does not expose that hybrid as a reportable
#' `stack_direct` estimate table.
#'
#' See also [pvstackr_object_contracts].
#'
#' @param data Analysis data frame.
#' @param formula Two-sided formula with `OUTCOME` on the left-hand side.
#' @param target A `pvstackr_brr_target` object from [pv_brr_target()].
#' @param control A [pv_control()] object with `method = "stack_direct"`.
#' @param family Optional family representation. Reportable `stack_direct`
#'   accepts only Gaussian/identity semantics and passes a package-owned
#'   canonical `stats::gaussian("identity")` object to the backend; executable
#'   closures supplied by the caller are never forwarded.
#' @param prior Optional backend prior object. Under the materialized
#'   `stack_direct` design, only recognized invariant prior tables are passed
#'   through; coefficient-, intercept-, scoped-, and opaque priors fail before
#'   backend execution because their identities cannot be preserved exactly.
#'   Accepted explicit priors are still reported as warning-level diagnostics
#'   because the current identity result is scoped to MLE/flat-prior
#'   fixed-effect regimes.
#' @param fit_function Injected backend fitting function for this package stage.
#' @param draws_function Function extracting posterior draws from the injected
#'   backend fit.
#' @param param_map Optional explicit draw-column map passed to the stack-fit
#'   layer. Supply `fe_names` or `fe_idx` to identify fixed-effect columns, and
#'   optional `vc_names` or `vc_idx` for nuisance variance-component columns.
#'   Use `vc_names = character()` to drop all nuisance columns. Explicit maps
#'   are recommended when backend draw names do not follow the automatic `b_*`
#'   fixed-effect convention, including distributional names such as
#'   `b_sigma_*`.
#' @param diagnose_function Optional injected backend diagnostic extractor.
#'   Its existing named-list fields are preserved, while sampler fields are
#'   normalized under `stack_fit$diagnostics$sampler`. A missing, malformed, or
#'   failed extractor produces explicit incomplete diagnostics and reason codes.
#'   When the bundled brms adapter is used, pvstackr automatically obtains
#'   R-hat and bulk/tail ESS from posterior summaries and divergences from brms
#'   NUTS parameters. The normalized record includes total and per-chain ESS
#'   minima, chain count, post-warmup draws per chain, source, and completeness.
#'   Extraction records diagnostics but does not itself change fit status.
#' @param log_lik_function Optional log-likelihood extractor.
#' @param extract_log_lik Whether to extract log-likelihood draws.
#' @param cache_dir,cache_stem Cache directory and safe file stem. For the
#'   bundled brms backend, pvstackr creates the directory, verifies that it is
#'   writable, and uses brms `file_refit = "on_change"`. Set `cache_dir = NULL`
#'   to disable file caching. Injected adapters continue to own their cache
#'   implementation and receive these values only as fit arguments.
#' @param additional_args Additional named arguments passed to the injected fit
#'   function.
#'
#' @returns A `pvstackr_fit` object.
#' @section Reportable scope and coverage:
#' In this package stage, reportable output is **fixed-effect-only**; variance
#' components are fit but not calibrated to the target. Coverage claims are
#' enabled **only** for `stack_direct` rows backed by the external Rubin/BRR-Fay
#' target with Barnard-Rubin df and an explicit `df_complete`
#' (`interval_role = "coverage_barnard_rubin"`, `coverage_claim_allowed = TRUE`);
#' with classic Rubin df the rows are descriptive
#' (`interval_role = "descriptive_classic_rubin"`, `coverage_claim_allowed = FALSE`),
#' as in the bundled cached example. `per_pv` and `stack_psis` intervals are
#' **descriptive/reference** even with Barnard-Rubin degrees of freedom. "One
#' stacked fit" describes the computational topology, not a benchmarked speed
#' claim.
#'
#' @examples
#' pisa_tiny <- read.csv(
#'   system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
#' )
#'
#' # Build the design and the external BRR-Fay target on the fixture.
#' design <- pv_design(
#'   pisa_tiny, formula = OUTCOME ~ x + female,
#'   pv_suffix = "READ", expected_M = 2L, expected_R = 4L, id_cols = "CNTSTUID"
#' )
#' target <- pv_brr_target(
#'   pisa_tiny, formula = OUTCOME ~ x + female,
#'   pv_cols = design$pv_cols, weight_col = design$weight_col,
#'   rep_weight_cols = design$rep_weight_cols, fay_k = design$fay_k,
#'   id_cols = design$id_cols
#' )
#'
#' # A live stack_direct fit needs injected/precomputed backend adapters
#' # (e.g. brms/cmdstanr); not run in checks.
#' \dontrun{
#' fit <- pv_fit_direct(
#'   data = pisa_tiny, formula = OUTCOME ~ x + female,
#'   target = target,
#'   fit_function = my_backend_adapter, draws_function = my_draws_adapter
#' )
#' }
#'
#' # Inspect the reportable fixed-effect surface via the bundled cached fit.
#' path <- system.file(
#'   "extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr"
#' )
#' if (nzchar(path)) {
#'   fit <- readRDS(path)$fit     # a stack_direct pvstackr_fit
#'   head(get_estimates(fit))     # reportable fixed-effect table
#' }
#' @references
#' Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys.* Wiley.
#'
#' Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with
#' multiple imputation. *Biometrika*, 86(4), 948-955.
#' @family pvstackr-fitting
#' @seealso [pv_fit()], [pv_fit_reference()], [pv_fit_stack_psis()],
#'   [pv_control()]; [pv_brr_target()], [pv_compare_methods()], [get_estimates()],
#'   [pvstackr_object_contracts].
#' @export
pv_fit_direct <- function(
  data,
  formula,
  target,
  control = pv_control(method = "stack_direct"),
  family = NULL,
  prior = NULL,
  fit_function = NULL,
  draws_function = NULL,
  param_map = NULL,
  diagnose_function = NULL,
  log_lik_function = NULL,
  extract_log_lik = FALSE,
  cache_dir = "cache",
  cache_stem = "pvstackr-stack-direct",
  additional_args = list()
) {
  control <- if (is.null(control)) pv_control(method = "stack_direct") else pv_validate_control(control)
  if (!identical(control$method, "stack_direct")) {
    pv_abort("`control$method` must be `stack_direct` for `pv_fit_direct()`.")
  }
  if (!identical(control$center, "target")) {
    pv_abort("Reportable `pv_fit_direct()` output requires `control$center = \"target\"`; `center = \"posterior\"` is reserved for CCC diagnostic/exploratory checks.")
  }
  pv_validate_fit_data_retention_control(control)
  canonical_family <- pv_compatibility_stack_direct_family(family)
  resolved_preflight <- pv_stack_direct_preflight(
    data = data,
    formula = formula,
    target = target,
    family = canonical_family,
    return_model_bundle = TRUE
  )
  preflight <- resolved_preflight$preflight
  design <- new_pvstackr_design(
    data = data,
    formula = formula,
    pv_cols = target$pv_cols,
    weight_col = target$weight_col,
    rep_weight_cols = target$rep_weight_cols,
    fay_k = target$fay_k,
    id_cols = target$id_cols,
    roles = list(
      outcome_placeholder = "OUTCOME",
      method = "stack_direct"
    ),
    provenance = list(
      source = "pv_fit_direct",
      target_hash = target$target_hash,
      target_manifest_hash = preflight$binding_proof$target_manifest_hash,
      target_content_hash = target$target_content$target_content_hash,
      binding_verification_policy = preflight$binding_proof$verification_policy
    )
  )
  design <- pv_design_canonicalize_formula(design)
  stack_control <- control
  stack_control$keep_data <- FALSE
  stack_control$return_draws <- TRUE
  stack_fit <- pv_stack_fit(
    data = data,
    formula = formula,
    pv_cols = target$pv_cols,
    weight_col = target$weight_col,
    control = stack_control,
    family = canonical_family,
    prior = prior,
    fit_function = fit_function,
    draws_function = draws_function,
    param_map = param_map,
    diagnose_function = diagnose_function,
    log_lik_function = log_lik_function,
    extract_log_lik = extract_log_lik,
    cache_dir = cache_dir,
    cache_stem = cache_stem,
    additional_args = additional_args,
    resolved_model_bundle = resolved_preflight$model_bundle,
    resolved_binding_manifest = target$binding_manifest
  )
  runtime_control <- control
  runtime_control$return_draws <- TRUE
  stack_fit <- pv_stack_fit_composite_projection(
    stack_fit = stack_fit,
    control = runtime_control,
    canonicalize_formula = FALSE
  )
  sampler <- stack_fit$diagnostics$sampler
  sampler_gate <- pv_sampler_gate(
    sampler,
    expected_chains = control$chains,
    expected_post_warmup_draws_per_chain = control$iter - control$warmup
  )
  if (identical(sampler_gate$status, "blocked")) {
    blocked_target <- pv_fit_direct_independent_target(target)
    blocked_control <- pv_fit_direct_blocked_control(control)
    return(new_pvstackr_fit(
      method = "stack_direct",
      target = blocked_target,
      estimates = data.frame(),
      diagnostics = list(
        preflight = pv_fit_direct_blocked_preflight(preflight),
        sampler = sampler,
        sampler_gate = sampler_gate,
        redaction = pv_fit_direct_blocked_redaction("sampler_gate")
      ),
      status = "blocked",
      control = blocked_control,
      reason_codes = sampler_gate$reason_codes,
      provenance = list(
        wrapper_function = "pv_fit_direct",
        target_hash = target$target_hash,
        target_manifest_hash = preflight$binding_proof$target_manifest_hash,
        target_content_hash = target$target_content$target_content_hash,
        binding_verification_policy = preflight$binding_proof$verification_policy,
        independent_target_retained = TRUE,
        sampler_diagnostic_source = sampler$diagnostic_source
      ),
      warnings = sampler_gate$warnings
    ))
  }
  ccc <- ccc_calibrate(
    draws = stack_fit$stacked_draws,
    target = target,
    param_map = stack_fit$param_map,
    center = control$center,
    binding_proof = preflight$binding_proof
  )
  status_info <- pv_fit_direct_status(stack_fit, ccc, sampler_gate)
  if (identical(status_info$status, "blocked")) {
    blocked_target <- pv_fit_direct_independent_target(target)
    blocked_ccc <- pv_fit_direct_blocked_ccc_diagnostics(stack_fit, ccc)
    return(new_pvstackr_fit(
      method = "stack_direct",
      target = blocked_target,
      estimates = data.frame(),
      diagnostics = list(
        preflight = pv_fit_direct_blocked_preflight(preflight),
        sampler = sampler,
        sampler_gate = sampler_gate,
        ccc = blocked_ccc,
        redaction = pv_fit_direct_blocked_redaction("ccc_reportability_gate")
      ),
      status = "blocked",
      control = pv_fit_direct_blocked_control(control),
      reason_codes = status_info$reason_codes,
      provenance = list(
        wrapper_function = "pv_fit_direct",
        target_hash = blocked_target$target_hash,
        target_manifest_hash = blocked_target$binding_manifest$manifest_hash,
        target_content_hash = blocked_target$target_content$target_content_hash,
        binding_verification_policy = preflight$binding_proof$verification_policy,
        independent_target_retained = TRUE,
        sampler_diagnostic_source = sampler$diagnostic_source,
        ccc_target_hash = blocked_ccc$target_hash,
        reportability_policy = "generic_blocked_fail_closed"
      ),
      warnings = status_info$warnings
    ))
  }
  estimates <- pv_fit_direct_estimates(ccc, target, control$conf_level)
  draws <- if (isTRUE(control$return_draws)) {
    ccc$draws_fe_cal
  } else {
    NULL
  }

  retained_design <- if (isTRUE(control$keep_data)) {
    design
  } else {
    pv_design_target_data_free_snapshot(
      target,
      preflight$binding_proof
    )
  }
  retained_preflight <- pv_fit_direct_retained_preflight(
    preflight,
    keep_data = control$keep_data
  )
  fit_control <- control
  nested_control <- fit_control
  nested_control$return_draws <- FALSE
  stack_fit <- pv_stack_fit_composite_projection(
    stack_fit = stack_fit,
    control = nested_control,
    canonicalize_formula = FALSE
  )
  ccc <- pv_ccc_draw_projection(ccc, FALSE)

  new_pvstackr_fit(
    method = "stack_direct",
    design = retained_design,
    target = target,
    stack_fit = stack_fit,
    ccc = ccc,
    estimates = estimates,
    draws = draws,
    diagnostics = list(
      preflight = retained_preflight,
      sampler = sampler,
      sampler_gate = sampler_gate,
      stack_fit = stack_fit$diagnostics,
      stack_fit_warnings = stack_fit$warnings,
      ccc = ccc$diagnostics
    ),
    status = status_info$status,
    control = fit_control,
    reason_codes = status_info$reason_codes,
    provenance = list(
      wrapper_function = "pv_fit_direct",
      target_hash = target$target_hash,
      target_manifest_hash = preflight$binding_proof$target_manifest_hash,
      target_content_hash = target$target_content$target_content_hash,
      binding_verification_policy = preflight$binding_proof$verification_policy,
      ccc_target_hash = ccc$target_hash,
      stack_fit_long_data_hash = stack_fit$weight_summary$long_data_hash
    ),
    warnings = status_info$warnings
  )
}
