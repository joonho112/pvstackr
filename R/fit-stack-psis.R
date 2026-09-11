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
    pv_abort("Random-effect/group terms are not supported for `stack_psis`.")
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

#' Reweight stacked draws toward each plausible value
#'
#' `pv_fit_stack_psis()` runs the `"stack_psis"` method of [pv_fit()]. It
#' takes one set of posterior draws from a model fitted to the stacked data
#' (one copy of the data per plausible value) and, for each plausible value,
#' importance weights for these draws and a Pareto k-hat value, computed
#' outside pvstackr, for example by Pareto smoothed importance sampling
#' (PSIS) with the loo package. It reweights the draws toward the
#' posterior of each plausible value and combines the weighted results with
#' Rubin's rules. pvstackr checks the weights and the k-hat values but does
#' not compute them.
#'
#' @details
#' ## Inputs
#'
#' The stacked draws come from exactly one of `stacked_draws`,
#' `fit_function` with `draws_function`, or `stack_fit`. The weights come
#' from `psis_weights` with `pareto_k`, or from your `psis_function`, which
#' pvstackr calls on `log_ratios` and which returns both. The weights have
#' one row per stacked draw and one column per plausible value. Name the
#' program that produced them with `psis_producer` and
#' `psis_producer_version`: the numbers alone cannot show that the weights
#' were Pareto-smoothed.
#'
#' PSIS stabilizes importance weights with a generalized Pareto distribution
#' fitted to the upper tail of the importance ratios; the Pareto k-hat, the
#' estimated shape parameter of that distribution, shows how reliable the
#' reweighted estimates are (Vehtari et al. 2024). pvstackr itself
#' does not run Pareto smoothing or any other part of PSIS: it fits no
#' Pareto distribution and estimates no k-hat.
#'
#' ## Calculation
#'
#' For each plausible value \eqn{m}, pvstackr divides the weights by their
#' sum and computes the weighted mean \eqn{\hat\beta_m}{b_m} and the
#' weighted covariance matrix \eqn{U_m} of the fixed-effect draws (see
#' `param_map`), with the covariance divided by
#' \eqn{1 - \sum w^2}{1 - sum(w^2)}. Rubin's rules then combine the \eqn{M}
#' results as in [pv_brr_target()]: the total covariance is
#' \eqn{T_{\mathrm{MI}} = \bar U + (1 + 1/M) B}{T_MI = U_bar + (1 + 1/M) B},
#' the standard errors are the square roots of its diagonal, and the degrees
#' of freedom follow `df_method`.
#'
#' ## Checks and status
#'
#' By pvstackr's rule, the fit is blocked (status `"blocked"`, with no
#' estimates, weights or draws) unless every plausible value has a finite
#' k-hat strictly below `control$psis_k_threshold` (default 0.7, the largest
#' value that [pv_control()] accepts) and the weights come with
#' `psis_producer` and `psis_producer_version`. The reason code says which
#' part failed: `psis_k_not_evaluated` (a k-hat is missing or not finite),
#' `psis_k_too_high`, `psis_weight_provenance_incomplete` (no program named)
#' or `psis_smoothing_not_applied` (weights that pvstackr made from
#' `log_ratios` without smoothing). Otherwise the status is `"ok"`; a
#' `stack_psis` fit never has status `"warning"`, and pvstackr does not check
#' the sampler diagnostics of the stacked fit.
#'
#' Vehtari et al. (2024) recommend the cut-off
#' \eqn{\min(1 - 1/\log_{10} S, 0.7)}{min(1 - 1 / log10(S), 0.7)} for
#' \eqn{S} draws, for example 0.67 for 1000 draws; to use it, set
#' `psis_k_threshold` in [pv_control()].
#'
#' ## What is reported
#'
#' As for every method, only the fixed effects are reported ([pv_fit()]
#' gives the scope and the full rule for the intervals). The variance within
#' each plausible value comes from the weighted posterior draws, not from
#' replicate weights, so every row of the estimate table has
#' `interval_role = "psis_classic_rubin"` or `"psis_barnard_rubin"` and
#' `coverage_claim_allowed = FALSE`. A k-hat below the cut-off lets pvstackr
#' report the estimates, but the intervals remain descriptive. pvstackr
#' provides `stack_psis` as a check on a stacked fit, not as a replacement
#' for `stack_direct`.
#'
#' `get_diagnostics(fit)$psis` also records, for a blocked fit too, how
#' concentrated the normalized weights are, and `weight_diagnostic_authority`
#' says whether the weights were kept (`"retained_weights_recomputed"`, for a
#' fit with estimates and `control$return_draws = TRUE`) or removed
#' (`"owned_stamp_bounded_projection"`). These values do not change the
#' status; [pvstackr_object_contracts] explains them.
#'
#' @param data A data frame, or `NULL` (default). Required with
#'   `fit_function`, which receives the stacked data (see `fit_function`).
#'   Otherwise it is optional: given together with `formula`, it is used only
#'   to record the design (the fit's `design`), and the plausible-value names
#'   must then be columns of `data`.
#' @param formula A two-sided formula with the placeholder `OUTCOME` on the
#'   left-hand side, such as `OUTCOME ~ x + female` (see [pv_fit()]), or
#'   `NULL` (default). Required with `fit_function`. A random-effect term
#'   such as `(1 | school)` or a `weights()` term stops the function with an
#'   error, also when the formula only serves to record the design.
#' @param pv_cols A character vector with the names of the plausible values,
#'   at least two, or `NULL` (default). Required with `fit_function`: the
#'   plausible-value columns of `data`. Otherwise it is optional and names the
#'   columns of the weights; if the weights have column names, they must be
#'   the same. Without `pv_cols`, the column names of the weights are used,
#'   or `PV1`, `PV2`, ... when they have none.
#' @param control A [pv_control()] object with `method = "stack_psis"`. The
#'   default is `pv_control(method = "stack_psis")`, and `NULL` gives the
#'   same. Its `psis_k_threshold` sets the k-hat cut-off, and its sampler
#'   settings and `backend` are only passed on to `fit_function`.
#' @param family,prior Passed unchanged to `fit_function`; pvstackr does not
#'   check them. Default `NULL`.
#' @param fit_function Your function that fits the model to the stacked
#'   data, or `NULL` (default). pvstackr stacks `data` (one copy per
#'   plausible value, with the plausible value in the column `.pvstackr_y`)
#'   and calls the function once with the arguments `formula`, `data`,
#'   `family`, `prior`, `chains`, `iter`, `warmup`, `cores`, `seed`,
#'   `backend`, `file` and `file_refit` (the stacked formula and data, the
#'   settings from `control` and the cache file; see `cache_dir`), plus the
#'   elements of `additional_args`. The formula has the right-hand side of
#'   `formula`, for example
#'   `.pvstackr_y | weights(.pvstackr_weight) ~ x + female` for
#'   `OUTCOME ~ x + female`. The row weight
#'   `.pvstackr_weight` is `weight_col` divided by its mean and by the number
#'   of plausible values, or 1/M without `weight_col`; a function that
#'   ignores it fits an unweighted model. There is no bundled engine for this
#'   method, also with `pv_control(backend = "brms")`.
#' @param draws_function Your function that takes the object returned by
#'   `fit_function` and returns its draws: a numeric matrix or data frame with
#'   one row per draw (at least two), finite values and unique column names.
#'   Required with `fit_function`.
#' @param stack_fit A stacked-fit object of class `pvstackr_stack_fit` that
#'   holds its stacked draws, or `NULL` (default). The `stack_fit` component
#'   of a fit returned by pvstackr keeps no stacked draws and stops the
#'   function with an error, so in practice give `stacked_draws` or
#'   `fit_function` instead.
#' @param stacked_draws Posterior draws of a model fitted to the stacked
#'   data, computed elsewhere, or `NULL` (default): a numeric matrix or data
#'   frame with one row per draw (at least two), finite values and unique
#'   column names.
#' @param param_map `NULL` (default) or a named list that says which columns
#'   of the stacked draws are the fixed effects, by name (`fe_names`) or by
#'   position (`fe_idx`), as for [pv_fit_direct()]. With `NULL`, the columns
#'   whose names start with `b_` are the fixed effects; with `stack_fit`, the
#'   selection recorded in `stack_fit` is used. Only the fixed effects are
#'   used, and their column names become the terms of the estimate table, so
#'   they must start with `b_`. Give a list to leave out other columns that
#'   start with `b_`, such as distributional parameters named `b_sigma_*`.
#'   With `fit_function` or `stack_fit`, give the columns by name
#'   (`fe_names`), because positions can point to other columns there.
#' @param psis_weights A numeric matrix of importance weights, one row per
#'   stacked draw and one column per plausible value, or `NULL` (default).
#'   The values must be finite and not negative, and every column must have a
#'   positive sum; pvstackr divides each column by its sum. Give it with
#'   `pareto_k`.
#' @param pareto_k The Pareto k-hat values, one per plausible value: a
#'   numeric vector named by plausible value, or in the order of the weight
#'   columns. Required with `psis_weights`, and with `log_ratios` when there
#'   is no `psis_function`; leave it `NULL` (default) with `psis_function`,
#'   which returns it.
#' @param log_ratios A numeric matrix of log importance ratios for
#'   reweighting the stacked draws toward the posterior of each plausible
#'   value, with one row per stacked draw, one column per plausible value and
#'   finite values; or `NULL` (default). With `psis_function`, pvstackr passes
#'   it to that function. With `pareto_k` alone, pvstackr turns the ratios
#'   into normalized weights without smoothing and the fit is blocked
#'   (`psis_smoothing_not_applied`), so this gives only the diagnostics in
#'   `get_diagnostics(fit)$psis`.
#' @param psis_function Your function that takes `log_ratios` and returns a
#'   list with `weights` (a matrix of the same dimensions, with the columns in
#'   the same order) and `pareto_k`, or `NULL` (default). For example, a
#'   function that applies the function `psis()` of the loo package to each
#'   column.
#' @param psis_producer,psis_producer_version The name and the version of
#'   the program that produced the weights and the k-hat values, such as
#'   `"loo"` and `as.character(packageVersion("loo"))`: single strings of at
#'   most 128 and 64 bytes, without leading or trailing spaces or control
#'   characters. Give both or neither; without them the fit is blocked
#'   (`psis_weight_provenance_incomplete`). pvstackr records them in the
#'   estimate table and in `get_diagnostics(fit)$psis` as your statement of
#'   how the weights were made, but it does not check them. They cannot be
#'   given with weights that pvstackr makes from `log_ratios` without
#'   `psis_function`.
#' @param fallback `"block"` (default) or `"warn"`. Both block a fit that
#'   fails the check (see Details). `"warn"` is kept so that older code still
#'   runs: it gives a deprecation warning and is only recorded, as
#'   `fallback_requested` in `get_diagnostics(fit)$psis`.
#' @param weight_col,rep_weight_cols,fay_k,id_cols The final weight column,
#'   the replicate-weight columns, the Fay coefficient and the row-identifier
#'   columns, as in [pv_design()]. `weight_col` sets the row weights for
#'   `fit_function` (see there). With `data` and `formula`, all four are
#'   checked and recorded in the fit's `design`; the replicate weights and
#'   `fay_k` are not used in the calculation. Defaults `NULL`, `NULL`, `0.5`
#'   and `NULL`.
#' @param df_method The rule for the degrees of freedom: `"classic"`
#'   (default) or `"barnard_rubin"`, which needs `df_complete`. Both rules are
#'   defined in [pv_brr_target()]. With either rule the intervals are
#'   descriptive.
#' @param df_complete For `df_method = "barnard_rubin"`, the complete-data
#'   degrees of freedom, which you state: one positive number for all fixed
#'   effects, or one value per fixed effect, preferably named by the
#'   fixed-effect columns (such as `b_Intercept`); an unnamed vector is taken
#'   in the order of the fixed-effect columns. Leave it `NULL` (default) with
#'   `"classic"`: a value there changes nothing but the `df_complete` column
#'   of the estimate table.
#' @param allow_m1 Has no effect: a `stack_psis` fit needs weights and k-hat
#'   values for at least two plausible values. Default `FALSE`.
#' @param cache_dir,cache_stem The folder and the file name that are passed
#'   to `fit_function`; defaults `"cache"` and `"pvstackr-stack-psis"`. It
#'   receives `file = file.path(path.expand(cache_dir), cache_stem)` and
#'   `file_refit = "on_change"`, or `file = NULL` and `file_refit = "never"`
#'   when `cache_dir = NULL`. pvstackr does not create the folder.
#' @param additional_args A named list of further arguments passed to
#'   `fit_function`, `list()` by default. They may not repeat the arguments
#'   that pvstackr sets (see `fit_function`).
#'
#' @returns A `pvstackr_fit` object with `method = "stack_psis"`. Read it
#'   with [get_estimates()] and [get_diagnostics()]; [get_target()] and
#'   [get_draws()] return `NULL` for this method, and
#'   [pvstackr_object_contracts] describes every part of the object. The
#'   estimate table has one row per fixed effect with the combined estimate,
#'   its standard error and its degrees of freedom; besides the columns of
#'   every method, it has `pooling_source`, `pooling_hash` and the columns
#'   `psis_status`, `pareto_k_max`, `psis_k_threshold`, `psis_source`,
#'   `pareto_k_source`, `weight_method`, `psis_producer` and
#'   `psis_producer_version`. [get_diagnostics()] returns `psis` (the k-hat
#'   values and their check, the sources of the weights and the
#'   weight-concentration values), `pooling` (the Rubin's-rules combination,
#'   such as `U_bar`, `B` and `T_MI`) and `weighted` (the weighted estimate
#'   and covariance of each plausible value and, with
#'   `control$return_draws = TRUE`, the fixed-effect draws and the normalized
#'   weights). For a blocked fit the estimate table is empty, and
#'   [get_diagnostics()] returns only `psis` and `redaction`, which lists what
#'   was removed.
#'
#' @examples
#' set.seed(1)
#' # Simulated stacked draws with the fixed-effect columns b_Intercept and
#' # b_x, equal placeholder weights with one column per plausible value, and
#' # k-hat values below the 0.7 cut-off. These weights were not made by
#' # PSIS, so no program is named, and pvstackr blocks the fit.
#' M <- 2L
#' stacked_draws <- matrix(
#'   rnorm(400 * 2), ncol = 2,
#'   dimnames = list(NULL, c("b_Intercept", "b_x"))
#' )
#' psis_weights <- matrix(1 / 400, nrow = 400, ncol = M)  # equal weights
#' pareto_k     <- rep(0.2, M)                             # below 0.7
#' fit_psis <- pv_fit_stack_psis(
#'   stacked_draws = stacked_draws,
#'   pv_cols       = paste0("PV", seq_len(M)),
#'   psis_weights  = psis_weights,
#'   pareto_k      = pareto_k,
#'   control       = pv_control(method = "stack_psis")
#' )
#' # Status "blocked", reason code psis_weight_provenance_incomplete:
#' fit_psis
#' # The PSIS check has status "provenance_incomplete", and weight_method is
#' # "unspecified_external" because no program was named:
#' get_diagnostics(fit_psis)$psis[c("status", "weight_method")]
#' @references
#' Vehtari, A., Simpson, D., Gelman, A., Yao, Y., & Gabry, J. (2024). Pareto
#' smoothed importance sampling. *Journal of Machine Learning Research*, 25(72),
#' 1-58.
#'
#' Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys.* Wiley.
#' @family pvstackr-fitting
#' @seealso [pv_compare_methods()] to compare fits,
#'   [pv_migrate_legacy_psis_fit()] for `stack_psis` fits saved by earlier
#'   versions of pvstackr, and [pvstackr_object_contracts] for the parts of a
#'   fit.
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
