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
  sweep(weights, 2L, sums, FUN = "/")
}

pv_psis_normalize_log_ratios <- function(log_ratios) {
  log_ratios <- as.matrix(log_ratios)
  if (!is.numeric(log_ratios) || any(!is.finite(log_ratios))) {
    pv_abort("`log_ratios` must be a finite numeric matrix.")
  }
  shifted <- sweep(log_ratios, 2L, apply(log_ratios, 2L, max), FUN = "-")
  pv_psis_normalize_weights(exp(shifted))
}

pv_psis_result <- function(psis_weights = NULL, pareto_k = NULL, log_ratios = NULL, psis_function = NULL) {
  if (!is.null(psis_weights)) {
    if (is.null(pareto_k)) {
      pv_abort("`pareto_k` is required when supplying `psis_weights` directly.")
    }
    return(list(weights = psis_weights, pareto_k = pareto_k, source = "supplied_weights"))
  }
  if (!is.null(psis_function)) {
    if (!is.function(psis_function)) {
      pv_abort("`psis_function` must be a function.")
    }
    if (is.null(log_ratios)) {
      pv_abort("`log_ratios` is required when using `psis_function`.")
    }
    out <- psis_function(log_ratios)
    pv_assert_named_list(out, "psis_function output")
    if (!all(c("weights", "pareto_k") %in% names(out))) {
      pv_abort("`psis_function` output must include `weights` and `pareto_k`.")
    }
    return(list(weights = out$weights, pareto_k = out$pareto_k, source = "psis_function"))
  }
  if (!is.null(log_ratios) && !is.null(pareto_k)) {
    return(list(weights = pv_psis_normalize_log_ratios(log_ratios), pareto_k = pareto_k, source = "log_ratios_self_normalized"))
  }
  pv_abort("Supply PSIS diagnostics through `psis_weights` plus `pareto_k`, or through `psis_function` plus `log_ratios`.")
}

pv_psis_align_columns <- function(weights, pareto_k, pv_cols = NULL) {
  weights <- pv_psis_normalize_weights(weights)
  M <- ncol(weights)
  if (is.null(pv_cols)) {
    pv_cols <- colnames(weights)
    if (is.null(pv_cols) || any(!nzchar(pv_cols))) {
      pv_cols <- paste0("PV", seq_len(M))
    }
  } else {
    pv_cols <- pv_validate_unique_character(pv_cols, "pv_cols", min_len = 2L)
    if (length(pv_cols) != M) {
      pv_abort("`pv_cols` length must match the number of PSIS weight columns.")
    }
  }
  colnames(weights) <- pv_cols

  if (!is.numeric(pareto_k) || length(pareto_k) != M) {
    pv_abort("`pareto_k` must be a numeric vector aligned to plausible values.")
  }
  if (is.null(names(pareto_k)) || all(!nzchar(names(pareto_k)))) {
    names(pareto_k) <- pv_cols
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
  weights <- pv_psis_normalize_weights(weights, n_draws = nrow(draws))
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
    weights = weights,
    param_map = list(
      fe_idx = map$fe_idx,
      vc_idx = map$vc_idx,
      fe_names = map$fe_names,
      vc_names = map$vc_names,
      dropped_names = map$dropped_names,
      map_source = map$map_source
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
    stringsAsFactors = FALSE
  )
}

pv_stack_psis_diagnostics <- function(pareto_k, threshold) {
  threshold <- pv_assert_scalar_number(
    threshold,
    "psis_k_threshold",
    lower = 0,
    upper = 1,
    inclusive_lower = FALSE,
    inclusive_upper = TRUE
  )
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

pv_stack_psis_pool_hash <- function(summary, pool, psis) {
  pv_hash_payload(list(
    pv_cols = summary$pv_cols,
    fe_names = summary$fe_names,
    beta = summary$beta,
    U = summary$U,
    pareto_k = psis$pareto_k,
    threshold = psis$threshold,
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
#' one stacked draw source, supplied/precomputed or injected per-plausible-value
#' PSIS weights, Pareto-k diagnostics, model-based Rubin pooling of
#' PSIS-weighted fixed-effect summaries, and explicit gating when Pareto-k
#' diagnostics fail.
#'
#' @details
#' This function does not depend on `loo`, `brms`, or `cmdstanr` directly and
#' does not run a live PSIS routine from the `loo` package by default. Users may
#' supply PSIS weights plus Pareto-k values, or inject a
#' `psis_function(log_ratios)` that returns `weights` and `pareto_k`.
#' In v0.1, reportable output is fixed-effect-only; group terms such as
#' `(1 | school)` and `(1 || school)` are rejected.
#'
#' `stack_psis` interval metadata is diagnostic/reference vocabulary, not a
#' design-coverage claim. Classic Rubin pooling is labeled
#' `interval_role = "psis_classic_rubin"`; Barnard-Rubin pooling is labeled
#' `interval_role = "psis_barnard_rubin"` and uses the supplied `df_complete`.
#' Both roles set `coverage_claim_allowed = FALSE` because the pooled covariance
#' is model-based weighted covariance of the stacked draws. The pooling
#' provenance is labeled `stack_psis_rubin_pooling`. PSIS diagnostics record
#' their weight source as `supplied_weights`, `psis_function`, or
#' `log_ratios_self_normalized`; the last path self-normalizes the log ratios
#' and does not by itself run Pareto smoothing.
#'
#' If any Pareto-k exceeds `control$psis_k_threshold`, the default
#' `fallback = "block"` returns a blocked `pvstackr_fit` with diagnostics but no
#' reportable estimates. `fallback = "warn"` retains estimates with warning
#' status for diagnostic comparison workflows.
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
#' @param psis_weights Optional normalized or unnormalized PSIS weight matrix,
#'   one column per plausible value.
#' @param pareto_k Pareto-k diagnostics aligned to plausible values.
#' @param log_ratios Optional log-ratio matrix used with `psis_function`, or
#'   self-normalized when `pareto_k` is supplied.
#' @param psis_function Optional function returning `weights` and `pareto_k`.
#' @param fallback Behavior when Pareto-k exceeds the threshold: `"block"` or
#'   `"warn"`.
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
#' # Injected stack_psis demo: a stacked draw matrix plus supplied PSIS weights
#' # and Pareto-k (the function does not run `loo` by default). Fixed-effect
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
#' fit_psis                    # a stack_psis pvstackr_fit
#' get_estimates(fit_psis)     # interval_role = "psis_*"; coverage_claim_allowed = FALSE
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
  df_method <- match.arg(df_method)
  control <- if (is.null(control)) pv_control(method = "stack_psis") else pv_validate_control(control)
  if (!identical(control$method, "stack_psis")) {
    pv_abort("`control$method` must be `stack_psis` for `pv_fit_stack_psis()`.")
  }
  pv_stack_psis_reject_group_formula(formula)

  stacked_source <- pv_stack_psis_source(stack_fit, stacked_draws, fit_function)
  if (identical(stacked_source, "injected_fit")) {
    if (is.null(data) || is.null(formula) || is.null(pv_cols)) {
      pv_abort("`data`, `formula`, and `pv_cols` are required when using an injected `fit_function`.")
    }
    stack_fit <- pv_stack_fit(
      data = data,
      formula = formula,
      pv_cols = pv_cols,
      weight_col = weight_col,
      control = control,
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
    pv_stack_psis_reject_group_rhs(stack_fit$formula[[3L]])
    stacked_draws <- stack_fit$stacked_draws
    if (is.null(param_map)) {
      param_map <- list(
        fe_idx = stack_fit$param_map$fe_idx,
        vc_idx = stack_fit$param_map$vc_idx
      )
    }
  }

  psis_raw <- pv_psis_result(
    psis_weights = psis_weights,
    pareto_k = pareto_k,
    log_ratios = log_ratios,
    psis_function = psis_function
  )
  aligned <- pv_psis_align_columns(psis_raw$weights, psis_raw$pareto_k, pv_cols = pv_cols)
  psis <- pv_stack_psis_diagnostics(aligned$pareto_k, control$psis_k_threshold)
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
  pooling_hash <- pv_stack_psis_pool_hash(summary, pool, psis)
  estimates <- pv_stack_psis_estimates(pool, psis, pooling_hash)

  failed <- !identical(psis$status, "ok")
  force_block <- identical(psis$status, "not_evaluated")
  status <- if (failed && (identical(fallback, "block") || force_block)) {
    "blocked"
  } else if (failed) {
    "warning"
  } else {
    "ok"
  }
  reason_codes <- if (failed) psis$reason_code else character()
  warnings <- if (failed && identical(psis$status, "not_evaluated")) {
    sprintf("PSIS Pareto-k was not fully evaluated for: %s.", paste(psis$bad_pv_cols, collapse = ", "))
  } else if (failed) {
    sprintf("PSIS Pareto-k exceeded threshold %.3f for: %s.", psis$threshold, paste(psis$bad_pv_cols, collapse = ", "))
  } else {
    character()
  }
  if (identical(status, "blocked")) {
    estimates <- data.frame()
  }

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
  }

  new_pvstackr_fit(
    method = "stack_psis",
    design = design,
    target = NULL,
    stack_fit = stack_fit,
    ccc = NULL,
    estimates = estimates,
    draws = NULL,
    diagnostics = list(
      psis = c(psis, list(source = psis_raw$source, fallback = fallback)),
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
        weights = if (isTRUE(control$return_draws)) summary$weights else NULL,
        param_map = summary$param_map
      )
    ),
    status = status,
    control = control,
    reason_codes = reason_codes,
    provenance = list(
      wrapper_function = "pv_fit_stack_psis",
      stacked_source = stacked_source,
      psis_source = psis_raw$source,
      pooling_hash = pooling_hash
    ),
    warnings = warnings
  )
}
