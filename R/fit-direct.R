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

#' Fit a model to plausible-value data
#'
#' `pv_fit()` fits a regression model to assessment data whose outcome is
#' given as plausible values, as in PISA, and returns estimates, standard
#' errors and intervals for the fixed effects. Each method has its own
#' function, whose help page lists the arguments that you pass through `...`:
#' [pv_fit_direct()] for `"stack_direct"`, [pv_fit_reference()] for
#' `"per_pv"` and [pv_fit_stack_psis()] for `"stack_psis"`.
#'
#' @details
#' ## Methods
#'
#' - `"stack_direct"` (default) fits one model to the stacked data (one copy
#'   of the data per plausible value), then applies the Cholesky calibration
#'   correction (CCC): the fixed-effect draws are transformed so that their
#'   mean and covariance equal those of `target`. The target comes from
#'   [pv_brr_target()], which fits a survey-weighted regression to each
#'   plausible value, estimates the sampling covariance of the coefficients
#'   with the BRR-Fay replicate weights, and combines the results with
#'   Rubin's rules.
#' - `"per_pv"` fits one model per plausible value and combines the results
#'   with Rubin's rules.
#' - `"stack_psis"` reweights one stacked fit toward each plausible value with
#'   importance weights and Pareto k-hat values that you supply (pvstackr
#'   checks them but does not compute them), then combines the reweighted
#'   results with Rubin's rules.
#'
#' Arguments in `...` are passed unchanged to the method's function. Besides
#' `data` and `formula`, which a method needs when it fits a model, each
#' method needs:
#'
#' - `"stack_direct"`: `target` and a model-fitting engine. Either set
#'   `pv_control(backend = "brms")`, which fits the stacked model with brms
#'   functions bundled with pvstackr (this needs the brms and posterior
#'   packages; if cmdstanr is installed, CmdStan must be configured or the
#'   fit stops with an error; without cmdstanr, brms uses rstan), or pass
#'   your own `fit_function`,
#'   `draws_function` and `diagnose_function` (without `diagnose_function`
#'   the fit is blocked). The plausible-value and weight columns are taken
#'   from `target`. See [pv_fit_direct()], also for `cache_dir`.
#' - `"per_pv"`: `pv_cols` with your own `fit_function` and `draws_function`,
#'   which are called once per plausible value and do not receive the survey
#'   weights automatically, or `per_pv_draws` computed elsewhere. There is no
#'   bundled engine for this method. See [pv_fit_reference()].
#' - `"stack_psis"`: one source of stacked draws (`stacked_draws`,
#'   `stack_fit`, or `fit_function` with `draws_function` and `pv_cols`);
#'   importance weights with Pareto k-hat values (`psis_weights` and
#'   `pareto_k`, or your `psis_function` applied to `log_ratios`); and
#'   `psis_producer` with `psis_producer_version`, which name the program
#'   that produced the weights. Without the last two, or when a k-hat value
#'   is not below `control$psis_k_threshold` (default 0.7), the fit is
#'   blocked. See [pv_fit_stack_psis()].
#'
#' When a method starts from draws computed elsewhere, pass `data = NULL` and
#' `formula = NULL` explicitly; `pv_fit()` has no defaults for them.
#'
#' ## What is reported
#'
#' Only the fixed effects are reported: the intercept and the slope
#' coefficients. Other model parameters, such as the residual standard
#' deviation, are not reported, and `"stack_direct"` does not calibrate them.
#' `"stack_direct"` requires `control$center = "target"` (the default): its
#' reported estimates, standard errors and degrees of freedom are those of
#' the target, and the stacked fit supplies the calibrated draws and the
#' diagnostics that compare the fit with the target.
#'
#' Each row of the estimate table has an interval (`conf_low`, `conf_high`)
#' at the level `conf_level` of [pv_control()] (default 0.95): the estimate
#' plus or minus a t quantile with `df` degrees of freedom times `se`.
#' `coverage_claim_allowed` says whether pvstackr's reporting rule lets you
#' read the interval as a confidence interval with nominal coverage (`TRUE`)
#' or labels it descriptive (`FALSE`), and `interval_role` names the case. The
#' label records how the interval was built; it does not certify its
#' coverage. pvstackr sets both columns by this rule:
#'
#' - `"stack_direct"` with a target built with Barnard-Rubin degrees of
#'   freedom (`df_method = "barnard_rubin"` and a `df_complete` value in
#'   [pv_brr_target()]): `interval_role = "coverage_barnard_rubin"` and
#'   `coverage_claim_allowed = TRUE`. These are the only intervals that the
#'   rule lets you read as confidence intervals with nominal coverage.
#' - `"stack_direct"` with the classic Rubin degrees of freedom
#'   (`df_method = "classic"`, the default): `"descriptive_classic_rubin"`
#'   and `FALSE`. The bundled example fit, computed from synthetic data, is
#'   of this kind.
#' - `"per_pv"` (`"reference_classic_rubin"`, `"reference_barnard_rubin"`)
#'   and `"stack_psis"` (`"psis_classic_rubin"`, `"psis_barnard_rubin"`):
#'   always `FALSE`, because their variance within each plausible value comes
#'   from the model's posterior draws, not from the replicate weights.
#'
#' These labels come from the target and the method, not from the status: a
#' fit with status `"warning"` keeps them, and a blocked fit has no rows.
#'
#' @param data A data frame with the plausible-value columns, the weight
#'   columns and the variables in `formula`, or `NULL` when the method starts
#'   from precomputed draws. For `"stack_direct"` it must be the data used to
#'   build `target` (the same rows in the same order, with the same plausible
#'   values, weights and covariates); otherwise `pv_fit()` stops with an
#'   error.
#' @param formula A two-sided formula with the placeholder `OUTCOME` on the
#'   left-hand side, for example `OUTCOME ~ x + female`. `OUTCOME` stands for
#'   the plausible values, so the model is written once. Use `NULL` when the
#'   method starts from precomputed draws. `"stack_direct"` and `"stack_psis"`
#'   stop with an error on random-effect terms such as `(1 | school)`. For
#'   `"stack_direct"` it must be the formula used to build `target`.
#' @param target For `"stack_direct"`, the `pvstackr_brr_target` object
#'   returned by [pv_brr_target()]; `pv_fit()` stops with an error if it is
#'   `NULL`. Ignored by `"per_pv"` and `"stack_psis"`. Default `NULL`.
#' @param method The fitting method: `"stack_direct"` (default), `"per_pv"` or
#'   `"stack_psis"`. See Details.
#' @param control A [pv_control()] object, or `NULL` (default) to use
#'   `pv_control(method = method)`. That default selects no fitting engine,
#'   so `"stack_direct"` stops with an error unless you set
#'   `pv_control(backend = "brms")` or pass your own `fit_function` (see
#'   Details). `pv_fit()` stops with an error if `control$method` differs
#'   from `method`.
#' @param ... Further arguments passed unchanged to the method's function (see
#'   Details).
#'
#' @returns A `pvstackr_fit` object. Read it with [get_estimates()],
#'   [get_target()], [get_draws()] and [get_diagnostics()] rather than with
#'   `$`; they stop with an error if the object was changed after it was
#'   created. Its `status` is `"ok"`, `"warning"` (estimates are returned and
#'   `warnings` says why) or `"blocked"` (the estimate table is empty);
#'   `reason_codes` gives the reasons for a warning or a block. A blocked fit
#'   does not stop with an error; in code, read the status as
#'   `summary(fit)$status` and the reasons as `summary(fit)$reason_codes`
#'   (`summary()` checks the fit as the reading functions do).
#'   [get_estimates()] returns one row per fixed effect with columns such as
#'   `term`, `estimate`, `se`, `df`, `conf_low`, `conf_high`, `interval_role`
#'   and `coverage_claim_allowed`.
#' @examples
#' pisa_tiny <- read.csv(
#'   system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
#' )
#'
#' # Declare the columns and build the BRR-Fay target from the bundled
#' # synthetic data.
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
#' # Fit with the bundled brms engine. This samples with Stan, so it is not
#' # run here.
#' \dontrun{
#' fit <- pv_fit(
#'   data = pisa_tiny, formula = OUTCOME ~ x + female,
#'   target = target, method = "stack_direct",
#'   control = pv_control(method = "stack_direct", backend = "brms")
#' )
#'
#' # With your own engine, pass all three functions (see ?pv_fit_direct).
#' # Without diagnose_function the fit is blocked.
#' fit <- pv_fit(
#'   data = pisa_tiny, formula = OUTCOME ~ x + female,
#'   target = target, method = "stack_direct",
#'   fit_function = my_fit_function, draws_function = my_draws_function,
#'   diagnose_function = my_diagnose_function, cache_dir = NULL
#' )
#' }
#'
#' # Read the example fit that ships with pvstackr. It was made from the same
#' # data with fitting functions that return draws around the target, not
#' # with a sampler.
#' path <- system.file(
#'   "extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr"
#' )
#' if (nzchar(path)) {
#'   fit <- readRDS(path)$fit     # a stack_direct fit (class pvstackr_fit)
#'   print(fit)                   # method, status and interval note
#'   head(get_estimates(fit))     # the fixed-effect estimate table
#' }
#' @family pvstackr-fitting
#' @seealso [pv_brr_target()] to build the target, [pv_compare_methods()] to
#'   compare fits, and [pvstackr_object_contracts] for the parts of a fit.
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

#' Fit a stacked model and calibrate it to a BRR-Fay target
#'
#' `pv_fit_direct()` runs the `"stack_direct"` method of [pv_fit()]. It fits
#' one model to the stacked data (one copy of `data` per plausible value),
#' checks the sampler diagnostics of that fit, and then applies the Cholesky
#' calibration correction (CCC): the fixed-effect draws are transformed so
#' that their mean and covariance equal the estimates and the covariance
#' (`T_MI`) of `target`. The reported estimates, standard errors and degrees
#' of freedom are those of `target`.
#'
#' @details
#' ## The stacked fit
#'
#' The stacked data hold one copy of `data` for each plausible value, with
#' that plausible value as the outcome. Each row is weighted by the final
#' survey weight named in `target`, divided by its mean and by the number of
#' plausible values, and the model is a Gaussian linear model with the
#' identity link. pvstackr adds the columns of the design matrix of
#' `formula`, the intercept included, to the stacked data as ordinary
#' covariates, so the fitting function receives the formula
#' `.pvstackr_y | weights(.pvstackr_weight) ~ 0 + pvstackrMM001 + ...`, in
#' which the weights enter through `weights()`. With the bundled brms engine
#' and `prior = NULL`, brms therefore gives every fixed effect, the intercept
#' included, its default flat prior. The residual standard deviation `sigma`
#' is estimated but neither calibrated nor reported.
#'
#' ## What is reported
#'
#' Only the fixed effects are reported. By pvstackr's reporting rule, the
#' intervals of a `stack_direct` fit can be read as confidence intervals with
#' nominal coverage only when `target` uses Barnard-Rubin degrees of freedom; [pv_fit()] gives the full
#' rule.
#'
#' `pv_fit_direct()` requires `control$center = "target"` (the default of
#' [pv_control()]), so the estimate table holds the numbers of the target:
#' `estimate` is the target estimate, `se` is `sqrt(diag(target$T_MI))`, and
#' `df`, `df_method`, `df_complete`, `interval_role` and
#' `coverage_claim_allowed` are copied from `target`: the degrees of freedom
#' follow the rule chosen in [pv_brr_target()] and are
#' not residual degrees of freedom estimated from the stacked fit. The
#' stacked fit supplies the draws that are calibrated and the diagnostics
#' described below. The other value, `center = "posterior"`, would keep the
#' mean of the stacked draws and calibrate only their covariance; it is
#' accepted by [pv_control()], but `pv_fit_direct()` stops with an error when
#' it is set.
#'
#' ## Fitting engine
#'
#' With `pv_control(backend = "brms")` and `fit_function = NULL`, the brms
#' engine bundled with pvstackr fits the stacked model. It needs the brms and
#' posterior packages, and [pv_backend_brms_fit_function()] describes how it
#' chooses the Stan backend. A `draws_function` that you pass replaces its
#' draws step, but a `diagnose_function` is ignored: the bundled engine
#' computes the sampler diagnostics itself. Otherwise, as with the default
#' `backend = "none"`, pass all three functions. A `fit_function` that you
#' pass is used whatever the `backend` value, which is only passed on to it.
#' There is no argument for draws computed elsewhere.
#'
#' ## Checks and status
#'
#' Before the model is fitted, `pv_fit_direct()` stops with an error when
#' `data` or `formula` differ from those used to build `target`, when
#' `formula` has a random-effect term such as `(1 | school)`, when `family`
#' is not Gaussian with the identity link, or when `prior` is not accepted.
#' After the stacked fit, pvstackr's checks give the fit the status `"ok"`,
#' `"warning"` or `"blocked"`:
#'
#' - Sampler diagnostics, checked first. The six values listed under
#'   `diagnose_function`, from the bundled engine or from your function, must
#'   all be available, and `chains` and `post_warmup_draws_per_chain` must
#'   equal `chains` and `iter - warmup` in `control`. The largest R-hat and
#'   the smallest bulk and tail effective sample sizes can give a warning or
#'   a block, and any divergent transition blocks the fit. A fit blocked here
#'   is not calibrated.
#' - Calibration. `delta_c_max`, the largest difference between a target
#'   estimate and the mean of the stacked draws in target standard errors,
#'   and `kappa_A`, the condition number of the calibration matrix, can each
#'   give a warning or a block.
#' - Priors. Any explicit `prior` gives a warning.
#'
#' The thresholds are listed under "Status and checks" in
#' [pvstackr_object_contracts]. The status is `"ok"` when no check gave a
#' warning or a block; it does not show that the sampler converged. A fit
#' with status `"warning"` has estimates, and `reason_codes` and `warnings`
#' name the checks. A blocked fit has an empty estimate table and no draws;
#' it keeps the target and the diagnostics that explain the block.
#'
#' @param data The data frame used to build `target`: the same rows in the
#'   same order, with the same plausible values, weights and covariates.
#'   Otherwise `pv_fit_direct()` stops with an error.
#' @param formula The formula used to build `target`, with the placeholder
#'   `OUTCOME` on the left-hand side, such as `OUTCOME ~ x + female` (see
#'   [pv_fit()]). Random-effect terms such as `(1 | school)` stop with an
#'   error.
#' @param target The `pvstackr_brr_target` object returned by
#'   [pv_brr_target()] for `data` and `formula`.
#' @param control A [pv_control()] object with `method = "stack_direct"`.
#'   The default is `pv_control(method = "stack_direct")`, and `NULL` gives
#'   the same. Its `center` must be `"target"`.
#' @param family `NULL` (default), `"gaussian"`, or a family object or list
#'   with elements `family` and `link`, such as `gaussian()`. Only the
#'   Gaussian family with the identity link is accepted; any other family or
#'   link stops with an error. The fitting function always receives
#'   `stats::gaussian()` made by pvstackr, not the object you pass.
#' @param prior `NULL` (default) or a prior table: a data frame with
#'   character columns `class` and `coef`, such as the result of the brms
#'   function `set_prior()`. Only priors of class `"b"` or `"sigma"` with an
#'   empty `coef` are accepted; a prior on the intercept, on one coefficient
#'   or of another class stops with an error. A `"b"` prior is given to the
#'   slopes only, as in the original formula, so the intercept keeps the
#'   default prior. Any explicit prior, even a flat one, gives the fit the
#'   status `"warning"` (`explicit_prior_warning`).
#' @param fit_function `NULL` (default) or your function that fits the
#'   stacked model. `NULL` selects the bundled brms engine when
#'   `control$backend` is `"brms"` and stops with an error otherwise. The
#'   function is called with the arguments `formula`, `data`, `family`,
#'   `prior`, `chains`, `iter`, `warmup`, `cores`, `seed`, `backend`, `file`
#'   and `file_refit` (the stacked formula and data, the settings from
#'   `control` and the cache file; see Details and `cache_dir`), plus the
#'   elements of `additional_args`. There is no `weights` argument: the
#'   weights are in `formula` and in the column `.pvstackr_weight` of `data`,
#'   and a function that ignores them fits an unweighted model. It may
#'   return any object that `draws_function` and `diagnose_function` accept.
#' @param draws_function Your function that takes the object returned by the
#'   fitting function and returns its draws: a numeric matrix or data frame
#'   with one row per draw (at least two), finite values and unique column
#'   names. It may be `NULL` (default) only with the bundled engine. The
#'   fixed-effect columns must be named as in `target$fe_names` (such as
#'   `b_Intercept`) or after the columns of the stacked formula
#'   (`b_pvstackrMM001`, `b_pvstackrMM002`, ...); otherwise give `param_map`
#'   and list all fixed effects in the order of `target$fe_names`.
#' @param param_map `NULL` (default) or a named list that says which columns
#'   of the draws are the fixed effects. With `NULL`, columns whose names
#'   start with `b_` are the fixed effects; `sigma` and `sd_*` columns are
#'   other model parameters, which are not reported, and all other columns,
#'   such as `lp__`, are dropped. Give a list with `fe_names` (names) or
#'   `fe_idx` (positions) when the fixed-effect columns have other names or
#'   when other columns start with `b_`, such as distributional parameters
#'   named `b_sigma_*`. Its optional elements `vc_names` and `vc_idx` select
#'   the other model parameters and change nothing that is reported.
#' @param diagnose_function Your function that takes the object returned by
#'   the fitting function and returns its sampler diagnostics as a named list
#'   with the elements `rhat_max` (the largest R-hat), `ess_bulk_min` and
#'   `ess_tail_min` (the smallest bulk and tail effective sample sizes, totals
#'   over all chains), `divergences` (the number of divergent transitions),
#'   `chains` and `post_warmup_draws_per_chain`; the list may also sit in an
#'   element `sampler`. For the per-chain check, pvstackr divides each
#'   effective sample size by `chains`; a returned `ess_bulk_per_chain_min`
#'   or `ess_tail_per_chain_min` must equal that ratio. If the function is
#'   missing or fails, or a value is missing or invalid, the fit is blocked
#'   (`sampler_diagnostics_incomplete`), and `diagnostic_reason_codes` in
#'   `get_diagnostics(fit)$sampler` says why. The bundled engine ignores this
#'   argument: it takes R-hat and the effective sample sizes over all
#'   variables of the fit from the posterior package and the number of
#'   divergent transitions from the brms sampler record.
#' @param log_lik_function `NULL` (default) or your function that takes the
#'   object returned by the fitting function and returns the pointwise
#'   log-likelihood: a finite numeric matrix with one row per draw and one
#'   column per row of the stacked data. It is called only when
#'   `extract_log_lik = TRUE`; the bundled engine does not supply one.
#' @param extract_log_lik Whether to call `log_lik_function`, which is then
#'   required: `TRUE` or `FALSE` (default). The matrix is kept in the fit's
#'   `stack_fit` only with `pv_control(keep_log_lik = TRUE)`, for your own
#'   use: no pvstackr function uses it. `keep_log_lik = TRUE` stops with an
#'   error unless `extract_log_lik = TRUE`.
#' @param cache_dir,cache_stem The folder and the file name (without
#'   extension) for a cached stacked fit; defaults `"cache"` and
#'   `"pvstackr-stack-direct"`. `cache_stem` must be a plain file name,
#'   without a path separator. For the bundled brms engine, pvstackr creates
#'   the folder, verifies that it is writable, and uses brms
#'   `file_refit = "on_change"`. Set `cache_dir = NULL` to disable file
#'   caching. Your own `fit_function` receives
#'   `file = file.path(path.expand(cache_dir), cache_stem)` and
#'   `file_refit = "on_change"`, or `file = NULL` and `file_refit = "never"`
#'   when `cache_dir = NULL`; pvstackr does not create the folder for it.
#'
#'   brms refits a cached model only when its Stan code or data change, not
#'   when `chains`, `iter`, `warmup` or `seed` change. The sampler check then
#'   blocks a cached fit whose number of chains or of draws per chain
#'   differs from `control`, and a changed seed returns the cached draws
#'   without a warning. Use `cache_dir = NULL` or a new `cache_stem` when you
#'   vary these settings, as in a seed-stability check.
#' @param additional_args A named list of further arguments for the fitting
#'   function, `list()` by default. They are passed to your `fit_function`
#'   or, with the bundled engine, to the brms function `brm()`, for example
#'   `list(control = list(adapt_delta = 0.95))`. They may not repeat the
#'   arguments that pvstackr sets (see `fit_function`).
#'
#' @returns A `pvstackr_fit` object with `method = "stack_direct"`. Read it
#'   with [get_estimates()], [get_target()], [get_draws()] and
#'   [get_diagnostics()]; [pv_fit()] explains the status and the estimate
#'   table, and [pvstackr_object_contracts] describes every part of the
#'   object. The estimate table has one row per fixed effect, with the
#'   estimates, standard errors and degrees of freedom of `target`.
#'   [get_draws()] returns the calibrated fixed-effect draws (one row per
#'   draw of the stacked fit, one column per fixed effect) when
#'   `control$return_draws` is `TRUE`, the default. [get_diagnostics()]
#'   returns `preflight` (the check that `data`, `formula` and `target`
#'   match), `sampler` and `sampler_gate` (the sampler diagnostics and their
#'   check), `stack_fit`, `stack_fit_warnings` (notes from the stacked fit,
#'   such as dropped draw columns) and `ccc` (the calibration diagnostics).
#'   For a blocked fit it returns only `preflight`, `sampler`,
#'   `sampler_gate`, `redaction` (what was removed) and, when the
#'   calibration check blocked the fit, `ccc`, with its values grouped as
#'   `center`, `conditioning`, `residual` and `prior` (for example
#'   `get_diagnostics(fit)$ccc$center$delta_c_max`).
#'
#' @examples
#' pisa_tiny <- read.csv(
#'   system.file("extdata", "pisa_tiny.csv", package = "pvstackr")
#' )
#'
#' # Declare the columns and build the BRR-Fay target from the bundled
#' # synthetic data.
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
#' # A live fit samples with Stan, so it is not run here. With
#' # backend = "brms" and no fit_function, the bundled brms engine is used.
#' \dontrun{
#' fit <- pv_fit_direct(
#'   data = pisa_tiny, formula = OUTCOME ~ x + female, target = target,
#'   control = pv_control(method = "stack_direct", backend = "brms")
#' )
#'
#' # With the default backend = "none", fit_function and draws_function are
#' # required, and without diagnose_function the fit is blocked.
#' fit <- pv_fit_direct(
#'   data = pisa_tiny, formula = OUTCOME ~ x + female, target = target,
#'   fit_function = my_fit_function, draws_function = my_draws_function,
#'   diagnose_function = my_diagnose_function, cache_dir = NULL
#' )
#' }
#'
#' # Read the example fit that ships with pvstackr. It was made from the same
#' # data with fitting functions that return draws around the target, not
#' # with a sampler.
#' path <- system.file(
#'   "extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr"
#' )
#' if (nzchar(path)) {
#'   fit <- readRDS(path)$fit     # a stack_direct fit (class pvstackr_fit)
#'   head(get_estimates(fit))     # the target's estimates, se and df
#' }
#' @references
#' Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys.* Wiley.
#'
#' Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with
#' multiple imputation. *Biometrika*, 86(4), 948-955.
#' @family pvstackr-fitting
#' @seealso [pv_brr_target()] to build the target,
#'   [pv_backend_brms_fit_function()] for the bundled brms engine, and
#'   [pvstackr_object_contracts] for the parts of a fit and the thresholds of
#'   the checks.
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
