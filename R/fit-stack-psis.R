pv_psis_normalize_weights <- function(weights, n_draws = NULL) {
  weights <- as.matrix(weights)
  if (!is.numeric(weights) || any(!is.finite(weights)) || any(weights < 0)) {
    pv_abort("`psis_weights` must be a finite non-negative numeric matrix.")
  }
  if (!is.null(n_draws) && nrow(weights) != n_draws) {
    pv_abort("`psis_weights` row count must match the stacked draw count.")
  }
  if (ncol(weights) < 2L) {
    pv_abort("`psis_weights` must contain at least two plausible-value columns.")
  }
  sums <- colSums(weights)
  if (any(!is.finite(sums)) || any(sums <= 0)) {
    pv_abort("Every `psis_weights` column must have positive total weight.")
  }
  normalized <- sweep(weights, 2L, sums, FUN = "/")
  matrix(
    as.numeric(normalized),
    nrow = nrow(normalized),
    ncol = ncol(normalized),
    dimnames = list(NULL, colnames(normalized))
  )
}

pv_psis_validate_log_ratios <- function(log_ratios) {
  if (!is.matrix(log_ratios) || !is.numeric(log_ratios) ||
      length(dim(log_ratios)) != 2L || nrow(log_ratios) < 2L ||
      ncol(log_ratios) < 2L || any(!is.finite(log_ratios))) {
    pv_abort("`log_ratios` must be a finite numeric matrix with at least two rows and two plausible-value columns.")
  }
  column_names <- colnames(log_ratios)
  if (!is.null(column_names) &&
      (anyNA(column_names) || any(!nzchar(column_names)) ||
        anyDuplicated(column_names))) {
    pv_abort("Named `log_ratios` columns must be unique and nonempty.")
  }
  matrix(
    as.numeric(log_ratios),
    nrow = nrow(log_ratios),
    ncol = ncol(log_ratios),
    dimnames = list(NULL, column_names)
  )
}

pv_psis_normalize_log_ratios <- function(log_ratios) {
  log_ratios <- pv_psis_validate_log_ratios(log_ratios)
  shifted <- sweep(log_ratios, 2L, apply(log_ratios, 2L, max), FUN = "-")
  pv_psis_normalize_weights(exp(shifted))
}

pv_psis_producer_record <- function(producer = NULL, producer_version = NULL) {
  if (is.null(producer) && is.null(producer_version)) {
    return(list(
      weight_method = "unspecified_external",
      producer = NA_character_,
      producer_version = NA_character_
    ))
  }
  if (is.null(producer) || is.null(producer_version)) {
    pv_abort("`psis_producer` and `psis_producer_version` must be supplied together.")
  }
  values <- list(
    psis_producer = producer,
    psis_producer_version = producer_version
  )
  limits <- c(psis_producer = 128L, psis_producer_version = 64L)
  for (field in names(values)) {
    value <- values[[field]]
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !nzchar(trimws(value)) || !identical(value, trimws(value)) ||
        !is.null(attributes(value)) ||
        nchar(value, type = "bytes") > limits[[field]] ||
        grepl("[[:cntrl:]]", value)) {
      pv_abort(sprintf("`%s` must be a bare, nonempty bounded scalar string.", field))
    }
  }
  list(
    weight_method = "caller_declared_external_psis",
    producer = producer,
    producer_version = producer_version
  )
}

pv_psis_result <- function(
  psis_weights = NULL,
  pareto_k = NULL,
  log_ratios = NULL,
  psis_function = NULL,
  psis_producer = NULL,
  psis_producer_version = NULL
) {
  if (!is.null(psis_weights)) {
    if (!is.null(log_ratios) || !is.null(psis_function)) {
      pv_abort("Supply exactly one weight route; direct `psis_weights` cannot be combined with `log_ratios` or `psis_function`.")
    }
    if (is.null(pareto_k)) {
      pv_abort("`pareto_k` is required when supplying `psis_weights` directly.")
    }
    producer <- pv_psis_producer_record(psis_producer, psis_producer_version)
    return(c(list(
      weights = psis_weights,
      pareto_k = pareto_k,
      source = "supplied_psis_weights",
      pareto_k_source = "supplied"
    ), producer))
  }
  if (!is.null(psis_function)) {
    if (!is.null(pareto_k)) {
      pv_abort("An injected `psis_function` must return `pareto_k`; do not also supply `pareto_k`.")
    }
    if (!is.function(psis_function)) {
      pv_abort("`psis_function` must be a function.")
    }
    if (is.null(log_ratios)) {
      pv_abort("`log_ratios` is required when using `psis_function`.")
    }
    log_ratios <- pv_psis_validate_log_ratios(log_ratios)
    out <- psis_function(log_ratios)
    pv_assert_named_list(out, "psis_function output")
    if (!all(c("weights", "pareto_k") %in% names(out))) {
      pv_abort("`psis_function` output must include `weights` and `pareto_k`.")
    }
    out_weights <- as.matrix(out$weights)
    if (!identical(dim(out_weights), dim(log_ratios))) {
      pv_abort("`psis_function` weights must preserve the complete log-ratio row and plausible-value dimensions.")
    }
    log_ratio_names <- colnames(log_ratios)
    output_names <- colnames(out_weights)
    if (!is.null(log_ratio_names) && !is.null(output_names) &&
        !identical(output_names, log_ratio_names)) {
      pv_abort("`psis_function` weight columns must preserve log-ratio column order and labels.")
    }
    if (is.null(output_names) && !is.null(log_ratio_names)) {
      colnames(out_weights) <- log_ratio_names
    }
    producer <- pv_psis_producer_record(psis_producer, psis_producer_version)
    return(c(list(
      weights = out_weights,
      pareto_k = out$pareto_k,
      source = "injected_psis_function",
      pareto_k_source = "injected_function_output"
    ), producer))
  }
  if (!is.null(log_ratios) && !is.null(pareto_k)) {
    if (!is.null(psis_producer) || !is.null(psis_producer_version)) {
      pv_abort("Self-normalized raw log ratios cannot be relabeled as externally produced PSIS weights.")
    }
    return(list(
      weights = pv_psis_normalize_log_ratios(log_ratios),
      pareto_k = pareto_k,
      source = "self_normalized_log_ratios",
      pareto_k_source = "supplied",
      weight_method = "self_normalized_raw_importance",
      producer = "pvstackr",
      producer_version = "0.2.0"
    ))
  }
  pv_abort("Supply PSIS diagnostics through `psis_weights` plus `pareto_k`, or through `psis_function` plus `log_ratios`.")
}

pv_psis_align_columns <- function(weights, pareto_k, pv_cols = NULL) {
  weights <- pv_psis_normalize_weights(weights)
  M <- ncol(weights)
  weight_names <- colnames(weights)
  if (is.null(pv_cols)) {
    if (is.null(weight_names)) {
      pv_cols <- paste0("PV", seq_len(M))
    } else {
      pv_cols <- pv_validate_unique_character(
        weight_names,
        "psis_weights column names",
        min_len = 2L
      )
    }
  } else {
    pv_cols <- pv_validate_unique_character(pv_cols, "pv_cols", min_len = 2L)
    if (length(pv_cols) != M) {
      pv_abort("`pv_cols` length must match the number of PSIS weight columns.")
    }
    if (!is.null(weight_names)) {
      if (anyNA(weight_names) || any(!nzchar(weight_names)) ||
          anyDuplicated(weight_names) || !setequal(weight_names, pv_cols)) {
        pv_abort("Named `psis_weights` columns must align exactly with `pv_cols`.")
      }
      weights <- weights[, pv_cols, drop = FALSE]
    }
  }
  colnames(weights) <- pv_cols

  if (!is.numeric(pareto_k) || length(pareto_k) != M) {
    pv_abort("`pareto_k` must be a numeric vector aligned to plausible values.")
  }
  if (is.null(names(pareto_k)) || all(!nzchar(names(pareto_k)))) {
    names(pareto_k) <- if (is.null(weight_names)) pv_cols else weight_names
    pareto_k <- pareto_k[pv_cols]
  } else if (any(!nzchar(names(pareto_k))) || anyDuplicated(names(pareto_k)) || !setequal(names(pareto_k), pv_cols)) {
    pv_abort("Named `pareto_k` values must align with `pv_cols`.")
  } else {
    pareto_k <- pareto_k[pv_cols]
  }
  list(weights = weights, pareto_k = pareto_k, pv_cols = pv_cols)
}

pv_weighted_mean_cov <- function(draws, weights) {
  if (!is.matrix(draws) || !is.numeric(draws)) {
    pv_abort("`draws` must be a numeric matrix.")
  }
  if (!is.numeric(weights) || length(weights) != nrow(draws) || any(!is.finite(weights)) || any(weights < 0)) {
    pv_abort("PSIS weights must align with draw rows.")
  }
  weights <- weights / sum(weights)
  center <- colSums(draws * weights)
  centered <- sweep(draws, 2L, center, FUN = "-")
  denom <- 1 - sum(weights^2)
  if (!is.finite(denom) || denom <= .Machine$double.eps) {
    pv_abort("PSIS weights are degenerate; weighted covariance cannot be computed.")
  }
  cov <- crossprod(centered * sqrt(weights)) / denom
  dimnames(cov) <- list(colnames(draws), colnames(draws))
  list(mean = center, cov = pv_symmetrize(cov))
}

pv_stack_psis_summarize <- function(stacked_draws, weights, param_map = NULL) {
  map <- pv_stack_param_map(stacked_draws, param_map = param_map)
  draws <- map$draws_selected[, map$fe_idx, drop = FALSE]
  weights <- as.matrix(weights)
  if (!is.numeric(weights) || nrow(weights) != nrow(draws) ||
      any(!is.finite(weights)) || any(weights < 0)) {
    pv_abort("PSIS weights must be a finite non-negative matrix aligned with stacked draws.")
  }
  if (any(abs(colSums(weights) - 1) > 1e-10)) {
    weights <- pv_psis_normalize_weights(weights, n_draws = nrow(draws))
  } else {
    weights <- matrix(
      as.numeric(weights),
      nrow = nrow(weights),
      ncol = ncol(weights),
      dimnames = list(NULL, colnames(weights))
    )
  }
  pv_cols <- colnames(weights)
  M <- ncol(weights)
  fe_names <- colnames(draws)
  beta <- matrix(NA_real_, nrow = M, ncol = length(fe_names), dimnames = list(pv_cols, fe_names))
  U <- vector("list", M)
  for (m in seq_len(M)) {
    out <- pv_weighted_mean_cov(draws, weights[, m])
    beta[m, ] <- out$mean
    U[[m]] <- out$cov
  }
  names(U) <- pv_cols
  list(
    pv_cols = pv_cols,
    fe_names = fe_names,
    beta = beta,
    U = U,
    proposal_draws = draws,
    weights = weights,
    param_map = list(
      fe_idx = map$fe_idx,
      vc_idx = integer(),
      fe_names = map$fe_names,
      vc_names = character(),
      dropped_names = character(),
      map_source = "fixed_effect_projection"
    )
  )
}

pv_stack_psis_estimates <- function(pool, psis, pooling_hash) {
  fe_names <- names(pool$beta)
  data.frame(
    term = fe_names,
    estimate = unname(pool$beta[fe_names]),
    se = unname(pool$se[fe_names]),
    std.error = unname(pool$se[fe_names]),
    df = unname(pool$df[fe_names]),
    df_method = pool$df_method,
    df_complete = unname(pool$df_complete[fe_names]),
    conf_level = pool$conf_level,
    conf_low = unname(pool$ci_low[fe_names]),
    conf_high = unname(pool$ci_high[fe_names]),
    conf.low = unname(pool$ci_low[fe_names]),
    conf.high = unname(pool$ci_high[fe_names]),
    interval_role = if (identical(pool$df_method, "barnard_rubin")) {
      "psis_barnard_rubin"
    } else {
      "psis_classic_rubin"
    },
    coverage_claim_allowed = FALSE,
    parameter_scope = "fixed_effect",
    target_source = "stack_psis_rubin_pooling",
    target_hash = pooling_hash,
    pooling_source = "stack_psis_rubin_pooling",
    pooling_hash = pooling_hash,
    psis_status = psis$status,
    pareto_k_max = psis$pareto_k_max,
    psis_k_threshold = psis$threshold,
    psis_source = psis$source,
    pareto_k_source = psis$pareto_k_source,
    weight_method = psis$weight_method,
    psis_producer = psis$producer,
    psis_producer_version = psis$producer_version,
    stringsAsFactors = FALSE
  )
}

pv_stack_psis_diagnostics <- function(pareto_k, threshold) {
  threshold <- pv_validate_psis_k_threshold(threshold)
  if (!is.numeric(pareto_k) || length(pareto_k) < 2L) {
    pv_abort("`pareto_k` must be numeric and aligned to plausible values.")
  }
  if (any(!is.finite(pareto_k))) {
    bad <- !is.finite(pareto_k)
    return(list(
      pareto_k = pareto_k,
      threshold = threshold,
      pareto_k_max = NA_real_,
      status = "not_evaluated",
      bad_pv_cols = names(pareto_k)[bad],
      reason_code = "psis_k_not_evaluated"
    ))
  }
  bad <- pareto_k >= threshold
  list(
    pareto_k = pareto_k,
    threshold = threshold,
    pareto_k_max = max(pareto_k),
    status = if (any(bad)) "failed" else "ok",
    bad_pv_cols = names(pareto_k)[bad],
    reason_code = if (any(bad)) "psis_k_too_high" else NA_character_
  )
}

pv_stack_psis_weight_diagnostics <- function(weights, pv_cols) {
  weights <- as.matrix(weights)
  if (!is.numeric(weights) || any(!is.finite(weights)) || any(weights < 0) ||
      any(abs(colSums(weights) - 1) > 1e-10)) {
    pv_abort("Weight diagnostics require finite non-negative column-normalized weights.")
  }
  if (!identical(colnames(weights), pv_cols)) {
    pv_abort("Normalized weight diagnostics must align with the declared plausible values.")
  }
  n_draws <- as.integer(nrow(weights))
  weight_ess_iid <- stats::setNames(
    as.numeric(1 / colSums(weights^2)),
    pv_cols
  )
  list(
    normalization = "column_sum_one",
    n_draws = n_draws,
    ess_definition = "kish_iid_normalized_weights_v1",
    weight_ess_iid = weight_ess_iid,
    weight_ess_fraction = stats::setNames(
      as.numeric(weight_ess_iid / n_draws),
      pv_cols
    ),
    max_normalized_weight = stats::setNames(
      as.numeric(apply(weights, 2L, max)),
      pv_cols
    )
  )
}

pv_stack_psis_apply_weight_gate <- function(psis, weight_method) {
  if (!identical(psis$status, "ok")) {
    return(psis)
  }
  if (identical(weight_method, "self_normalized_raw_importance")) {
    psis$status <- "unsmoothed"
    psis$reason_code <- "psis_smoothing_not_applied"
    return(psis)
  }
  if (identical(weight_method, "unspecified_external")) {
    psis$status <- "provenance_incomplete"
    psis$reason_code <- "psis_weight_provenance_incomplete"
  }
  psis
}

pv_stack_psis_blocked_warning <- function(psis) {
  if (identical(psis$status, "not_evaluated")) {
    return(sprintf(
      "PSIS Pareto-k was not fully evaluated for: %s.",
      paste(psis$bad_pv_cols, collapse = ", ")
    ))
  }
  if (identical(psis$status, "failed")) {
    return(sprintf(
      "PSIS Pareto-k met or exceeded threshold %.3f for: %s.",
      psis$threshold,
      paste(psis$bad_pv_cols, collapse = ", ")
    ))
  }
  if (identical(psis$status, "unsmoothed")) {
    return("PSIS smoothing was not applied; self-normalized raw importance weights are diagnostic-only.")
  }
  if (identical(psis$status, "provenance_incomplete")) {
    return("External PSIS weight provenance is incomplete; declare both producer and producer version.")
  }
  pv_abort("Blocked PSIS diagnostics have an unsupported status.")
}

pv_stack_psis_pool_hash <- function(summary, pool, psis) {
  pv_hash_payload(list(
    pv_cols = summary$pv_cols,
    fe_names = summary$fe_names,
    beta = summary$beta,
    U = summary$U,
    pareto_k = psis$pareto_k,
    threshold = psis$threshold,
    weight_source = psis$source,
    weight_method = psis$weight_method,
    producer = psis$producer,
    producer_version = psis$producer_version,
    weight_ess_iid = psis$weight_ess_iid,
    max_normalized_weight = psis$max_normalized_weight,
    pooled_beta = pool$beta,
    T_MI = pool$T_MI
  ))
}

pv_stack_psis_source <- function(stack_fit, stacked_draws, fit_function) {
  supplied <- c(!is.null(stack_fit), !is.null(stacked_draws), !is.null(fit_function))
  if (sum(supplied) != 1L) {
    pv_abort("Supply exactly one stacked source: `stack_fit`, `stacked_draws`, or `fit_function`.")
  }
  if (!is.null(stack_fit)) "stack_fit" else if (!is.null(stacked_draws)) "stacked_draws" else "injected_fit"
}

pv_stack_psis_blocked_redaction <- function() {
  list(
    status = "withheld",
    policy = "immutable_psis_fail_closed",
    withheld = c(
      "stack_fit", "pooling", "weighted", "beta", "U", "T_MI",
      "se", "df", "weights", "draws", "proposal_draws"
    )
  )
}

pv_stack_psis_reject_group_rhs <- function(rhs) {
  if (pv_formula_has_random_effect_bar(rhs)) {
    pv_abort("Random-effect/group terms are not supported for `stack_psis` in v0.1.")
  }
  invisible(FALSE)
}

pv_stack_psis_reject_group_formula <- function(formula) {
  if (is.null(formula)) {
    return(invisible(FALSE))
  }
  rhs <- pv_formula_rhs_checked(formula)
  pv_stack_psis_reject_group_rhs(rhs)
}

#' Fit the PSIS-Reweighted Stacked Method
#'
#' `pv_fit_stack_psis()` implements the `stack_psis` API for this package stage:
#' one stacked draw source, externally produced per-plausible-value weights and
#' Pareto-k diagnostics, model-based Rubin pooling of weighted fixed-effect
#' summaries, and explicit provenance and diagnostic gates.
#'
#' @details
#' This function does not depend on `loo`, `brms`, or `cmdstanr` directly and
#' does not run a live PSIS routine from the `loo` package by default. Users may
#' supply weights plus Pareto-k values, or inject a
#' `psis_function(log_ratios)` that returns `weights` and `pareto_k`. Numeric
#' weights or a function name cannot prove that Pareto smoothing occurred.
#' Reportable output therefore also requires both `psis_producer` and
#' `psis_producer_version`; these fields are caller-declared provenance, not
#' package verification. In v0.2, reportable output is fixed-effect-only; group terms such as
#' `(1 | school)` and `(1 || school)` are rejected.
#'
#' `stack_psis` interval metadata is diagnostic/reference vocabulary, not a
#' design-coverage claim. Classic Rubin pooling is labeled
#' `interval_role = "psis_classic_rubin"`; Barnard-Rubin pooling is labeled
#' `interval_role = "psis_barnard_rubin"` and uses the supplied `df_complete`.
#' Both roles set `coverage_claim_allowed = FALSE` because the pooled covariance
#' is model-based weighted covariance of the stacked draws. The pooling
#' provenance is labeled `stack_psis_rubin_pooling`. Diagnostics separate the
#' input route (`supplied_psis_weights`, `injected_psis_function`, or
#' `self_normalized_log_ratios`), the Pareto-k source, and the weight method.
#' The self-normalized path does not run Pareto smoothing and is always blocked
#' from estimates. Supplied or injected weights without producer/version
#' provenance are likewise diagnostic-only.
#'
#' Every path records per-PV weight-concentration diagnostics after column
#' normalization: `weight_ess_iid = 1 / sum(w^2)`, its draw-count fraction, and
#' the largest normalized weight. This is a Kish-style iid weight diagnostic,
#' not MCMC ESS, the `loo` relative-efficiency diagnostic, or an autocorrelation-adjusted PSIS
#' effective sample size. These quantities explain weight concentration but do
#' not relax or replace the Pareto-k gate.
#'
#' Reportability uses an immutable fail-closed gate: every plausible value must
#' have a finite Pareto-k strictly below `control$psis_k_threshold`, and that
#' threshold cannot exceed `0.7`. Failed, unevaluated, unsmoothed, or
#' provenance-incomplete input always returns a blocked `pvstackr_fit` with no
#' reportable estimates, raw weights, or draws; bounded ESS and maximum-weight
#' diagnostics remain available. The legacy
#' `fallback = "warn"` argument remains accepted for call compatibility and
#' emits a deprecation warning, but it cannot relax the gate and is recorded
#' only as the requested policy.
#'
#' When `control$return_draws = TRUE`, retained normalized weights are used to
#' recompute the ESS, maximum-weight, and weighted-summary diagnostics during
#' deep validation (`weight_diagnostic_authority =
#' "retained_weights_recomputed"`). When weights are intentionally redacted,
#' these bounded diagnostics are protected by the package-owned validation
#' stamp and cross-field feasibility checks
#' (`"owned_stamp_bounded_projection"`), but the original weight vector cannot
#' be independently reconstructed from the compact object. This is an explicit
#' portability and threat-model boundary, not a stronger statistical claim.
#'
#' @param data Optional analysis data frame for the injected stacked-fit route.
#' @param formula Optional formula with `OUTCOME` on the left-hand side for the
#'   injected stacked-fit route.
#' @param pv_cols Plausible-value columns.
#' @param control A [pv_control()] object with `method = "stack_psis"`.
#' @param family,prior Optional backend arguments passed to an injected stacked
#'   fit.
#' @param fit_function,draws_function Injected stacked backend fit and draw
#'   extractor.
#' @param stack_fit Optional existing `pvstackr_stack_fit`.
#' @param stacked_draws Optional stacked posterior draw matrix.
#' @param param_map Optional explicit draw-column map. Supply `fe_names` or
#'   `fe_idx` to identify fixed-effect columns, and optional `vc_names` or
#'   `vc_idx` for nuisance variance-component columns. Use
#'   `vc_names = character()` to drop all nuisance columns. Explicit maps are
#'   recommended when backend draw names do not follow the automatic `b_*`
#'   fixed-effect convention, including distributional names such as
#'   `b_sigma_*`.
#' @param psis_weights Optional normalized or unnormalized external weight
#'   matrix, one column per plausible value. The numeric matrix alone does not
#'   establish that PSIS smoothing occurred.
#' @param pareto_k Pareto-k diagnostics aligned to plausible values.
#' @param log_ratios Optional log-ratio matrix used with `psis_function`, or
#'   self-normalized when `pareto_k` is supplied.
#' @param psis_function Optional function returning `weights` and `pareto_k`.
#' @param psis_producer,psis_producer_version Optional bounded scalar strings
#'   identifying the external PSIS producer and its version. Both are required
#'   together for supplied or injected weights to be reportable. They record a
#'   caller declaration; package-verified `loo` execution is deferred.
#' @param fallback Requested behavior when Pareto-k fails: `"block"` or the
#'   legacy `"warn"`. Both choices now fail closed; `"warn"` cannot make failed
#'   PSIS output reportable.
#' @param weight_col,rep_weight_cols,fay_k,id_cols Optional design metadata.
#' @param df_method,df_complete,allow_m1 Rubin pooling options.
#' @param cache_dir,cache_stem Cache metadata for injected stacked fits.
#' @param additional_args Additional named arguments passed to injected fits.
#'
#' @returns A `pvstackr_fit` object with `method = "stack_psis"`.
#' @section Reportable scope and coverage:
#' In this package stage, reportable output is **fixed-effect-only**; variance
#' components are fit but not calibrated to the target. Coverage claims are
#' enabled **only** for `stack_direct` rows backed by the external Rubin/BRR-Fay
#' target (`interval_role = "coverage_barnard_rubin"`,
#' `coverage_claim_allowed = TRUE`). `stack_psis` intervals are
#' **descriptive/reference** even with Barnard-Rubin degrees of freedom. "One
#' stacked fit" describes the computational topology, not a benchmarked speed
#' claim. A small Pareto-\eqn{\hat k} does **not** certify correct variance:
#' reusing one draw cloud correlates the per-PV imputations, so intervals can run
#' narrow; \eqn{\hat k} is specification-dependent. Treat `stack_psis` as a
#' cross-check, never the deliverable.
#'
#' @examples
#' set.seed(1)
#' # Fail-closed stack_psis demo: equal placeholder weights and arbitrary
#' # Pareto-k values do not establish that PSIS ran. This example deliberately
#' # omits producer/version, so it remains diagnostic-only.
#' # Fixed-effect
#' # columns are named like `b_*`; weights have one column per plausible value.
#' M <- 2L
#' stacked_draws <- matrix(
#'   rnorm(400 * 2), ncol = 2,
#'   dimnames = list(NULL, c("b_Intercept", "b_x"))
#' )
#' psis_weights <- matrix(1 / 400, nrow = 400, ncol = M)  # equal weights
#' pareto_k     <- rep(0.2, M)                             # all "good" (< 0.5)
#' fit_psis <- pv_fit_stack_psis(
#'   stacked_draws = stacked_draws,
#'   pv_cols       = paste0("PV", seq_len(M)),
#'   psis_weights  = psis_weights,
#'   pareto_k      = pareto_k,
#'   control       = pv_control(method = "stack_psis")
#' )
#' fit_psis                    # blocked: provenance_incomplete
#' get_diagnostics(fit_psis)$psis[c("status", "weight_method")]
#' @references
#' Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC. *Statistics and
#' Computing*, 27(5), 1413-1432.
#'
#' Vehtari, A., Simpson, D., Gelman, A., Yao, Y., & Gabry, J. (2024). Pareto
#' smoothed importance sampling. *Journal of Machine Learning Research*, 25(72),
#' 1-58.
#'
#' Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys.* Wiley.
#' @family pvstackr-fitting
#' @seealso [pv_fit()], [pv_fit_direct()], [pv_fit_reference()],
#'   [pv_control()]; [pv_compare_methods()], [get_estimates()],
#'   [pvstackr_object_contracts].
#' @export
pv_fit_stack_psis <- function(
  data = NULL,
  formula = NULL,
  pv_cols = NULL,
  control = pv_control(method = "stack_psis"),
  family = NULL,
  prior = NULL,
  fit_function = NULL,
  draws_function = NULL,
  stack_fit = NULL,
  stacked_draws = NULL,
  param_map = NULL,
  psis_weights = NULL,
  pareto_k = NULL,
  log_ratios = NULL,
  psis_function = NULL,
  psis_producer = NULL,
  psis_producer_version = NULL,
  fallback = c("block", "warn"),
  weight_col = NULL,
  rep_weight_cols = NULL,
  fay_k = 0.5,
  id_cols = NULL,
  df_method = c("classic", "barnard_rubin"),
  df_complete = NULL,
  allow_m1 = FALSE,
  cache_dir = "cache",
  cache_stem = "pvstackr-stack-psis",
  additional_args = list()
) {
  fallback <- match.arg(fallback)
  if (identical(fallback, "warn")) {
    warning(
      paste(
        "`fallback = \"warn\"` is deprecated and behaves as `\"block\"`;",
        "it cannot make failed PSIS output reportable."
      ),
      call. = FALSE
    )
  }
  df_method <- match.arg(df_method)
  control <- if (is.null(control)) pv_control(method = "stack_psis") else pv_validate_control(control)
  if (!identical(control$method, "stack_psis")) {
    pv_abort("`control$method` must be `stack_psis` for `pv_fit_stack_psis()`.")
  }
  pv_validate_fit_data_retention_control(control)
  pv_stack_psis_reject_group_formula(formula)

  stacked_source <- pv_stack_psis_source(stack_fit, stacked_draws, fit_function)
  if (identical(stacked_source, "injected_fit")) {
    if (is.null(data) || is.null(formula) || is.null(pv_cols)) {
      pv_abort("`data`, `formula`, and `pv_cols` are required when using an injected `fit_function`.")
    }
    stack_control <- control
    stack_control$keep_data <- FALSE
    stack_control$return_draws <- TRUE
    stack_fit <- pv_stack_fit(
      data = data,
      formula = formula,
      pv_cols = pv_cols,
      weight_col = weight_col,
      control = stack_control,
      family = family,
      prior = prior,
      fit_function = fit_function,
      draws_function = draws_function,
      param_map = param_map,
      cache_dir = cache_dir,
      cache_stem = cache_stem,
      additional_args = additional_args
    )
    stacked_draws <- stack_fit$stacked_draws
  } else if (identical(stacked_source, "stack_fit")) {
    validate_pvstackr_stack_fit(stack_fit)
    if (!identical(stack_fit$schema_version, "0.2.0")) {
      pv_abort(
        "Legacy stack-fit objects are inspection-only and cannot enter a current stack_psis fit."
      )
    }
    pv_stack_psis_reject_group_rhs(stack_fit$formula[[3L]])
    if (is.null(stack_fit$stacked_draws)) {
      pv_abort(
        "A supplied stack-fit must retain stacked draws to enter `pv_fit_stack_psis()`; refit it with `return_draws = TRUE`."
      )
    }
    stacked_draws <- stack_fit$stacked_draws
    if (is.null(param_map)) {
      param_map <- list(
        fe_idx = stack_fit$param_map$fe_idx,
        vc_idx = stack_fit$param_map$vc_idx
      )
    }
  }

  if (!is.null(log_ratios)) {
    log_ratios <- pv_psis_validate_log_ratios(log_ratios)
    if (is.null(dim(stacked_draws)) || length(dim(stacked_draws)) != 2L ||
        nrow(log_ratios) != nrow(stacked_draws)) {
      pv_abort("`log_ratios` rows must align with the complete stacked draw rows before any injected function is called.")
    }
    if (!is.null(pv_cols)) {
      declared_pv_cols <- pv_validate_unique_character(
        pv_cols,
        "pv_cols",
        min_len = 2L
      )
      if (ncol(log_ratios) != length(declared_pv_cols)) {
        pv_abort("`log_ratios` columns must align with the complete declared plausible-value universe.")
      }
      log_ratio_names <- colnames(log_ratios)
      if (is.null(log_ratio_names)) {
        colnames(log_ratios) <- declared_pv_cols
      } else {
        if (!setequal(log_ratio_names, declared_pv_cols)) {
          pv_abort("Named `log_ratios` columns must align exactly with `pv_cols`.")
        }
        log_ratios <- log_ratios[, declared_pv_cols, drop = FALSE]
      }
    }
  }

  psis_raw <- pv_psis_result(
    psis_weights = psis_weights,
    pareto_k = pareto_k,
    log_ratios = log_ratios,
    psis_function = psis_function,
    psis_producer = psis_producer,
    psis_producer_version = psis_producer_version
  )
  aligned <- pv_psis_align_columns(psis_raw$weights, psis_raw$pareto_k, pv_cols = pv_cols)
  if (is.null(dim(stacked_draws)) || length(dim(stacked_draws)) != 2L ||
      nrow(aligned$weights) != nrow(stacked_draws)) {
    pv_abort("Every weight column must align with the complete stacked draw rows.")
  }
  psis <- pv_stack_psis_diagnostics(aligned$pareto_k, control$psis_k_threshold)
  psis <- pv_stack_psis_apply_weight_gate(psis, psis_raw$weight_method)
  weight_diagnostics <- pv_stack_psis_weight_diagnostics(
    aligned$weights,
    aligned$pv_cols
  )
  failed <- !identical(psis$status, "ok")
  weight_diagnostic_authority <- if (!failed && isTRUE(control$return_draws)) {
    "retained_weights_recomputed"
  } else {
    "owned_stamp_bounded_projection"
  }
  psis_diagnostics <- c(psis, list(
    pv_cols = aligned$pv_cols,
    source = psis_raw$source,
    pareto_k_source = psis_raw$pareto_k_source,
    weight_method = psis_raw$weight_method,
    producer = psis_raw$producer,
    producer_version = psis_raw$producer_version
  ), weight_diagnostics, list(
    weight_diagnostic_authority = weight_diagnostic_authority,
    fallback_requested = fallback,
    fallback_effective = if (failed) "block" else fallback
  ))
  if (failed) {
    blocked_control <- pv_fit_blocked_control(control)
    warnings <- pv_stack_psis_blocked_warning(psis)
    return(new_pvstackr_fit(
      method = "stack_psis",
      design = NULL,
      target = NULL,
      stack_fit = NULL,
      ccc = NULL,
      estimates = data.frame(),
      draws = NULL,
      diagnostics = list(
        psis = psis_diagnostics,
        redaction = pv_stack_psis_blocked_redaction()
      ),
      status = "blocked",
      control = blocked_control,
      reason_codes = psis$reason_code,
      provenance = list(
        wrapper_function = "pv_fit_stack_psis",
        stacked_source = stacked_source,
        psis_source = psis_raw$source,
        reportability_policy = "immutable_psis_fail_closed"
      ),
      warnings = warnings
    ))
  }
  summary <- pv_stack_psis_summarize(
    stacked_draws = stacked_draws,
    weights = aligned$weights,
    param_map = param_map
  )
  if (!identical(summary$pv_cols, aligned$pv_cols)) {
    pv_abort("PSIS weight columns must align with plausible-value labels.")
  }
  pool <- rubin_pool_matrix(
    beta = summary$beta,
    U = summary$U,
    orientation = "rows_pv",
    conf_level = control$conf_level,
    allow_m1 = allow_m1,
    df_method = df_method,
    df_complete = df_complete
  )
  pooling_hash <- pv_stack_psis_pool_hash(summary, pool, psis_diagnostics)
  estimates <- pv_stack_psis_estimates(pool, psis_diagnostics, pooling_hash)

  design <- NULL
  if (!is.null(data) && !is.null(formula) && !is.null(aligned$pv_cols)) {
    design <- new_pvstackr_design(
      data = data,
      formula = formula,
      pv_cols = aligned$pv_cols,
      weight_col = weight_col,
      rep_weight_cols = rep_weight_cols,
      fay_k = fay_k,
      id_cols = id_cols,
      roles = list(outcome_placeholder = "OUTCOME", method = "stack_psis"),
      provenance = list(source = "pv_fit_stack_psis", pooling_hash = pooling_hash)
    )
    design <- pv_design_canonicalize_formula(design)
    if (!isTRUE(control$keep_data)) {
      design <- pv_design_data_free_snapshot(design)
    }
  }

  nested_control <- control
  nested_control$return_draws <- FALSE
  stack_fit <- pv_stack_fit_composite_projection(
    stack_fit = stack_fit,
    control = nested_control,
    canonicalize_formula = TRUE
  )

  new_pvstackr_fit(
    method = "stack_psis",
    design = design,
    target = NULL,
    stack_fit = stack_fit,
    ccc = NULL,
    estimates = estimates,
    draws = NULL,
    diagnostics = list(
      psis = psis_diagnostics,
      pooling = list(
        beta = pool$beta,
        U_bar = pool$U_bar,
        B = pool$B,
        T_MI = pool$T_MI,
        lambda = pool$lambda,
        df = pool$df,
        df_classic = pool$df_classic,
        df_method = pool$df_method,
        df_complete = pool$df_complete,
        pooling_hash = pooling_hash,
        pooling_source = "stack_psis_rubin_pooling"
      ),
      weighted = list(
        beta = summary$beta,
        U = summary$U,
        proposal_draws = if (isTRUE(control$return_draws)) {
          summary$proposal_draws
        } else {
          NULL
        },
        weights = if (isTRUE(control$return_draws)) summary$weights else NULL,
        param_map = summary$param_map
      )
    ),
    status = "ok",
    control = control,
    reason_codes = character(),
    provenance = list(
      wrapper_function = "pv_fit_stack_psis",
      stacked_source = stacked_source,
      psis_source = psis_raw$source,
      pooling_hash = pooling_hash
    ),
    warnings = character()
  )
}
