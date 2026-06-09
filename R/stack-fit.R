pv_stack_reserved_cols <- function() {
  c(".pvstackr_y", ".pvstackr_pv", ".pvstackr_row", ".pvstackr_weight")
}

pv_stack_formula_rhs <- function(formula) {
  rhs <- pv_formula_rhs_checked(formula)
  if (pv_formula_has_weights_call(rhs)) {
    pv_abort("Do not embed `weights()` in `formula`; stacked weights are constructed internally.")
  }
  paste(deparse(rhs, width.cutoff = 500L), collapse = "")
}

pv_stack_weight_vector <- function(data, weight_col, M) {
  n <- nrow(data)
  if (is.null(weight_col)) {
    return(list(weights = rep(1 / M, n), source = "constant_fractional", normalized_base = rep(1, n)))
  }
  weight_col <- pv_assert_scalar_string(weight_col, "weight_col")
  pv_validate_columns(data, weight_col, "weight_col")
  base <- data[[weight_col]]
  if (!is.numeric(base) || any(!is.finite(base))) {
    pv_abort("`weight_col` must contain finite numeric weights.")
  }
  if (any(base <= 0)) {
    pv_abort("`weight_col` must contain strictly positive weights.")
  }
  mean_base <- mean(base)
  if (!is.finite(mean_base) || mean_base <= 0) {
    pv_abort("`weight_col` must have a positive finite mean.")
  }
  normalized <- base / mean_base
  list(weights = normalized / M, source = weight_col, normalized_base = normalized)
}

pv_hash_stack_data <- function(data) {
  pv_hash_payload(data[pv_stack_reserved_cols()])
}

pv_prepare_stack_data <- function(
  data,
  formula,
  pv_cols,
  weight_col = NULL
) {
  if (!is.data.frame(data)) {
    pv_abort("`data` must be a data frame.")
  }
  reserved <- intersect(pv_stack_reserved_cols(), names(data))
  if (length(reserved) > 0L) {
    pv_abort(sprintf("`data` already contains reserved internal column(s): %s.", paste(reserved, collapse = ", ")))
  }
  pv_cols <- pv_validate_pv_columns(data, pv_cols)
  rhs <- pv_stack_formula_rhs(formula)
  M <- length(pv_cols)
  n <- nrow(data)
  weight_info <- pv_stack_weight_vector(data, weight_col, M)

  long_parts <- vector("list", M)
  for (m in seq_along(pv_cols)) {
    pv <- pv_cols[[m]]
    part <- data
    part$.pvstackr_y <- part[[pv]]
    part$.pvstackr_pv <- factor(pv, levels = pv_cols)
    part$.pvstackr_row <- seq_len(n)
    part$.pvstackr_weight <- weight_info$weights
    long_parts[[m]] <- part
  }
  long_data <- do.call(rbind, long_parts)
  rownames(long_data) <- NULL

  long_formula <- stats::as.formula(
    sprintf(".pvstackr_y | weights(.pvstackr_weight) ~ %s", rhs),
    env = environment(formula)
  )
  per_pv_weight_sum <- vapply(pv_cols, function(pv) {
    sum(long_data$.pvstackr_weight[long_data$.pvstackr_pv == pv])
  }, numeric(1))

  list(
    data = long_data,
    formula = long_formula,
    formula_string = paste(deparse(long_formula, width.cutoff = 500L), collapse = ""),
    pv_cols = pv_cols,
    weight_summary = list(
      topology = "single_long_fit",
      M = M,
      n_original = n,
      n_long = nrow(long_data),
      fractional_weight = 1 / M,
      weight_source = weight_info$source,
      weight_col = weight_col,
      mean_long_weight = mean(long_data$.pvstackr_weight),
      total_long_weight = sum(long_data$.pvstackr_weight),
      per_pv_weight_sum = per_pv_weight_sum,
      long_data_hash = pv_hash_stack_data(long_data),
      long_data_hash_columns = pv_stack_reserved_cols()
    )
  )
}

pv_stack_build_fit_args <- function(
  prepared,
  family = NULL,
  prior = NULL,
  chains = 4L,
  iter = 2000L,
  warmup = 1000L,
  cores = 1L,
  seed = NULL,
  backend = "none",
  cache_dir = "cache",
  cache_stem = "pvstackr-stack",
  additional_args = list()
) {
  if (!is.list(prepared) || !all(c("formula", "data") %in% names(prepared))) {
    pv_abort("`prepared` must come from `pv_prepare_stack_data()`.")
  }
  if (!is.list(additional_args) ||
      (length(additional_args) > 0L && is.null(names(additional_args))) ||
      any(!nzchar(names(additional_args)))) {
    pv_abort("`additional_args` must be a named list.")
  }
  if (length(additional_args) > 0L && anyDuplicated(names(additional_args))) {
    pv_abort("`additional_args` must have unique names.")
  }
  protected <- c(
    "formula", "data", "family", "prior", "chains", "iter", "warmup",
    "cores", "seed", "backend", "file", "file_refit"
  )
  overridden <- intersect(names(additional_args), protected)
  if (length(overridden) > 0L) {
    pv_abort(sprintf("`additional_args` may not override protected fit argument(s): %s.", paste(overridden, collapse = ", ")))
  }
  c(
    list(
      formula = prepared$formula,
      data = prepared$data,
      family = family,
      prior = prior,
      chains = chains,
      iter = iter,
      warmup = warmup,
      cores = cores,
      seed = seed,
      backend = backend,
      file = file.path(cache_dir, cache_stem),
      file_refit = "on_change"
    ),
    additional_args
  )
}

pv_stack_param_ref <- function(param_map, idx_field, names_field) {
  has_idx <- !is.null(param_map[[idx_field]])
  has_names <- !is.null(param_map[[names_field]])
  if (has_idx && has_names) {
    pv_abort(sprintf("`param_map` must not supply both `%s` and `%s`.", idx_field, names_field))
  }
  param_map[[idx_field]] %||% param_map[[names_field]]
}

pv_stack_param_ref_idx <- function(ref, all_names, field, min_len = 1L) {
  if (is.character(ref)) {
    values <- pv_validate_unique_character(ref, field, min_len = min_len)
    idx <- match(values, all_names)
    if (anyNA(idx)) {
      pv_abort(sprintf("Draw column(s) not found for `%s`: %s.", field, paste(values[is.na(idx)], collapse = ", ")))
    }
    return(idx)
  }
  if (!is.integer(ref) && !is.numeric(ref)) {
    pv_abort(sprintf("`%s` must be character names or integer indexes.", field))
  }
  if (length(ref) < min_len || anyNA(ref) || any(ref != floor(ref))) {
    pv_abort(sprintf("`%s` must identify unique draw columns.", field))
  }
  idx <- as.integer(ref)
  if (length(idx) > 0L && (any(idx < 1L) || any(idx > length(all_names)) || anyDuplicated(idx))) {
    pv_abort(sprintf("`%s` must identify unique draw columns.", field))
  }
  idx
}

pv_stack_param_drop_warning <- function(dropped_names, map_source) {
  if (length(dropped_names) == 0L) {
    return(character())
  }
  source_label <- if (identical(map_source, "explicit")) "explicit" else "automatic regex"
  sprintf(
    "Stacked draw column(s) were not selected by the %s parameter map and were dropped: %s.",
    source_label,
    paste(dropped_names, collapse = ", ")
  )
}

pv_stack_param_map <- function(draws, param_map = NULL) {
  draws <- ccc_as_draw_matrix(draws)
  all_names <- colnames(draws)
  if (is.null(param_map)) {
    fe_idx <- grep("^b_", all_names)
    vc_idx <- grep("^(sd_|sigma$)", all_names)
    map_source <- "auto_regex"
  } else {
    pv_assert_named_list(param_map, "param_map")
    fe_ref <- pv_stack_param_ref(param_map, "fe_idx", "fe_names")
    if (is.null(fe_ref)) {
      pv_abort("`param_map` must contain `fe_idx` or `fe_names`.")
    }
    fe_idx <- pv_stack_param_ref_idx(fe_ref, all_names, "fe_idx/fe_names", min_len = 1L)

    vc_ref <- pv_stack_param_ref(param_map, "vc_idx", "vc_names")
    if (is.null(vc_ref)) {
      vc_idx <- setdiff(grep("^(sd_|sigma$)", all_names), fe_idx)
    } else {
      vc_idx <- pv_stack_param_ref_idx(vc_ref, all_names, "vc_idx/vc_names", min_len = 0L)
    }
    map_source <- "explicit"
  }
  keep_idx <- c(fe_idx, vc_idx)
  if (length(fe_idx) == 0L) {
    pv_abort("Stacked draw matrix must contain fixed-effect columns, or `param_map` must identify them.")
  }
  if (length(intersect(fe_idx, vc_idx)) > 0L) {
    pv_abort("Stacked fixed-effect and variance-component draw columns must not overlap.")
  }
  dropped_idx <- setdiff(seq_along(all_names), keep_idx)
  dropped_names <- all_names[dropped_idx]
  list(
    draws_selected = draws[, keep_idx, drop = FALSE],
    fe_idx = seq_along(fe_idx),
    vc_idx = if (length(vc_idx) > 0L) seq.int(length(fe_idx) + 1L, length(keep_idx)) else integer(),
    fe_names = all_names[fe_idx],
    vc_names = all_names[vc_idx],
    original_fe_idx = fe_idx,
    original_vc_idx = vc_idx,
    original_dropped_idx = dropped_idx,
    dropped_names = dropped_names,
    map_source = map_source,
    warnings = pv_stack_param_drop_warning(dropped_names, map_source)
  )
}

pv_validate_log_lik <- function(log_lik, n_draws, n_obs) {
  if (!is.matrix(log_lik) || !is.numeric(log_lik) || any(!is.finite(log_lik))) {
    pv_abort("`log_lik` must be a finite numeric matrix.")
  }
  if (!identical(dim(log_lik), c(n_draws, n_obs))) {
    pv_abort(sprintf("`log_lik` must have dimensions %d x %d.", n_draws, n_obs))
  }
  log_lik
}

pv_stack_vc_policy <- function() {
  list(
    policy_id = "VC-STACK-01",
    status = "nuisance_pass_through_not_calibrated_not_reported",
    calibration_status = "not_calibrated",
    validation_status = "not_validated",
    reporting_status = "not_reported_pending_validation",
    confirmatory_reporting_allowed = FALSE,
    pipeline_consumption_allowed = FALSE,
    reportable_parameter_scope = "fixed_effects_only",
    reportable_parameter_regex = "^b_"
  )
}

pv_stack_prior_warnings <- function(prior) {
  if (is.null(prior)) {
    return(character())
  }
  "Explicit priors were supplied; pvstackr does not verify prior flatness, and stack_direct identity diagnostics assume an MLE/flat-prior fixed-effect regime."
}

pv_stack_prior_policy <- function(prior) {
  prior_warning <- pv_stack_prior_warnings(prior)
  list(
    identity_scope = "mle_or_flat_prior_fixed_effect_regime",
    explicit_prior_supplied = !is.null(prior),
    explicit_prior = !is.null(prior),
    flatness = if (is.null(prior)) "not_supplied" else "unknown_explicit",
    explicit_prior_warning = length(prior_warning) > 0L,
    warn_explicit_prior = length(prior_warning) > 0L,
    reason_code = if (length(prior_warning) > 0L) "explicit_prior_warning" else NA_character_,
    warning = prior_warning
  )
}

pv_stack_direct_fe_names <- function(data, formula) {
  rhs <- pv_formula_rhs_checked(formula)
  model_formula <- stats::as.formula(
    paste("~", pv_deparse_expr(rhs)),
    env = environment(formula)
  )
  frame <- tryCatch(
    stats::model.frame(model_formula, data = data, na.action = stats::na.fail),
    error = function(e) pv_abort(sprintf(
      "Stack-direct model frame contains missing or invalid values; check missing data, factor levels, and contrasts: %s",
      conditionMessage(e)
    ))
  )
  x <- tryCatch(
    stats::model.matrix(stats::terms(model_formula), data = frame),
    error = function(e) pv_abort(sprintf(
      "Stack-direct fixed-effect names could not be derived; check factor levels, contrasts, and formula/data consistency: %s",
      conditionMessage(e)
    ))
  )
  pv_normalize_fe_names(colnames(x))
}

pv_stack_direct_preflight <- function(
  data,
  formula,
  target
) {
  if (!is.data.frame(data)) {
    pv_abort("`data` must be a data frame.")
  }
  validate_pvstackr_brr_target(target)
  if (!identical(target$target_source, "external_brr_fay_rubin")) {
    pv_abort("`target` must be an external BRR-Fay Rubin target.")
  }
  if (!isTRUE(target$policy$fixed_effects_only)) {
    pv_abort("`target` must be fixed-effect-only for `stack_direct`.")
  }
  if (!identical(target$policy$target_repair, "forbidden")) {
    pv_abort("`target` repair policy must be forbidden for `stack_direct`.")
  }

  rhs <- pv_formula_rhs_checked(formula)
  if (pv_formula_has_weights_call(rhs)) {
    pv_abort("Do not embed `weights()` in `formula`; stacked weights are constructed internally.")
  }
  if (pv_formula_has_random_effect_bar(rhs)) {
    pv_abort("Random-effect/group terms are not supported for `stack_direct` BRR target compatibility in v0.1.")
  }
  rhs_string <- pv_deparse_expr(rhs)
  if (!identical(rhs_string, target$rhs_string)) {
    pv_abort(sprintf(
      "`stack_direct` formula RHS must match the BRR target formula RHS exactly and in the same order. Fit RHS: `%s`; target RHS: `%s`. Rebuild the target with the same formula used for fitting.",
      rhs_string,
      target$rhs_string
    ))
  }

  fe_names <- pv_stack_direct_fe_names(
    data = data,
    formula = formula
  )
  if (!identical(fe_names, target$fe_names)) {
    pv_abort(sprintf(
      "`stack_direct` fixed-effect names must match BRR target fixed-effect names. Derived names from current data/formula: %s. Target names: %s. Check factor levels, contrasts, and formula/data consistency, then rebuild the target if needed.",
      paste(fe_names, collapse = ", "),
      paste(target$fe_names, collapse = ", ")
    ))
  }

  out <- list(
    formula = formula,
    formula_string = pv_formula_string(formula),
    rhs_string = rhs_string,
    fe_names = fe_names,
    target_hash = target$target_hash,
    target_source = target$target_source,
    policy = list(
      fixed_effects_only = target$policy$fixed_effects_only,
      target_repair = target$policy$target_repair
    )
  )
  class(out) <- c("pvstackr_stack_direct_preflight", "list")
  out
}

pv_stack_fit <- function(
  data,
  formula,
  pv_cols,
  weight_col = NULL,
  control = pv_control(),
  family = NULL,
  prior = NULL,
  fit_function = NULL,
  draws_function = NULL,
  param_map = NULL,
  diagnose_function = NULL,
  log_lik_function = NULL,
  extract_log_lik = FALSE,
  cache_dir = "cache",
  cache_stem = "pvstackr-stack",
  additional_args = list()
) {
  control <- pv_validate_control(control)
  if (is.null(fit_function)) {
    pv_abort("`fit_function` is required in this package stage; live brms/cmdstanr integration is added later.")
  }
  if (is.null(draws_function)) {
    pv_abort("`draws_function` is required when using an injected `fit_function`.")
  }
  extract_log_lik <- pv_assert_scalar_logical(extract_log_lik, "extract_log_lik")

  prepared <- pv_prepare_stack_data(
    data = data,
    formula = formula,
    pv_cols = pv_cols,
    weight_col = weight_col
  )
  fit_args <- pv_stack_build_fit_args(
    prepared = prepared,
    family = family,
    prior = prior,
    chains = control$chains,
    iter = control$iter,
    warmup = control$warmup,
    cores = control$cores,
    seed = control$seed,
    backend = control$backend,
    cache_dir = cache_dir,
    cache_stem = cache_stem,
    additional_args = additional_args
  )
  fit <- do.call(fit_function, fit_args)
  draws <- as.matrix(draws_function(fit))
  map <- pv_stack_param_map(draws, param_map = param_map)
  stacked_draws <- map$draws_selected
  param_map <- list(
    fe_idx = map$fe_idx,
    vc_idx = map$vc_idx,
    fe_names = map$fe_names,
    vc_names = map$vc_names,
    original_fe_idx = map$original_fe_idx,
    original_vc_idx = map$original_vc_idx,
    original_dropped_idx = map$original_dropped_idx,
    dropped_names = map$dropped_names,
    map_source = map$map_source
  )
  diagnostics <- if (is.null(diagnose_function)) list() else diagnose_function(fit)
  log_lik_extracted <- NULL
  if (extract_log_lik) {
    if (is.null(log_lik_function)) {
      pv_abort("`log_lik_function` is required when `extract_log_lik = TRUE`.")
    }
    log_lik_extracted <- pv_validate_log_lik(
      as.matrix(log_lik_function(fit)),
      n_draws = nrow(stacked_draws),
      n_obs = nrow(prepared$data)
    )
  }
  log_lik <- if (isTRUE(control$keep_log_lik)) log_lik_extracted else NULL
  vc_policy <- pv_stack_vc_policy()
  prior_policy <- pv_stack_prior_policy(prior)
  warnings <- c(pv_stack_prior_warnings(prior), map$warnings)

  out <- list(
    stacked_draws = stacked_draws,
    diagnostics = diagnostics,
    log_lik = log_lik,
    psi_hat_fe = colMeans(stacked_draws[, param_map$fe_idx, drop = FALSE]),
    param_map = param_map,
    formula = prepared$formula,
    formula_string = prepared$formula_string,
    weight_summary = prepared$weight_summary,
    meta = list(
      topology = "single_long_fit",
      engine_id = "single_long_fit",
      fit_engine = "injected_fit_function",
      n_fits = 1L,
      long_data_rows = nrow(prepared$data),
      long_data_hash = prepared$weight_summary$long_data_hash,
      log_lik_extracted = !is.null(log_lik_extracted),
      log_lik_retained = !is.null(log_lik),
      vc_policy = vc_policy,
      vc_policy_id = vc_policy$policy_id,
      vc_confirmatory_reporting_allowed = vc_policy$confirmatory_reporting_allowed,
      reportable_parameter_scope = vc_policy$reportable_parameter_scope,
      prior_policy = prior_policy,
      prior_diagnostic = prior_policy,
      param_map_source = param_map$map_source,
      dropped_draw_columns = param_map$dropped_names
    ),
    fit = if (isTRUE(control$keep_backend_fit)) fit else NULL,
    prepared_data = if (isTRUE(control$keep_data)) prepared$data else NULL,
    control = control,
    schema_version = "0.1.0",
    provenance = list(
      function_name = "pv_stack_fit",
      package = "pvstackr",
      topology = "single_long_fit",
      method_stage = "stack_fit",
      pv_cols = prepared$pv_cols,
      weight_col = weight_col,
      long_data_hash = prepared$weight_summary$long_data_hash,
      long_data_hash_columns = prepared$weight_summary$long_data_hash_columns,
      formula_string = prepared$formula_string,
      backend = control$backend,
      backend_fit_retained = isTRUE(control$keep_backend_fit),
      log_lik_retained = isTRUE(control$keep_log_lik),
      schema_version = "0.1.0"
    ),
    warnings = warnings
  )
  class(out) <- c("pvstackr_stack_fit", "list")
  out
}
