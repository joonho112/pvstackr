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

pv_hash_stack_data <- function(data, columns = pv_stack_reserved_cols()) {
  pv_hash_payload(data[columns])
}

pv_stack_materialized_design <- function(model_bundle) {
  model_matrix <- model_bundle$model_matrix
  column_names <- sprintf("pvstackrMM%03d", seq_len(ncol(model_matrix)))
  reportable_fe_names <- paste0("b_", colnames(model_matrix))
  reportable_fe_names[reportable_fe_names == "b_(Intercept)"] <- "b_Intercept"
  list(
    column_names = column_names,
    backend_fe_names = paste0("b_", column_names),
    reportable_fe_names = reportable_fe_names,
    offset_name = if (is.null(model_bundle$offset)) NULL else "pvstackrOffset"
  )
}

pv_prepare_stack_data <- function(
  data,
  formula,
  pv_cols,
  weight_col = NULL,
  model_bundle = NULL
) {
  if (!is.data.frame(data)) {
    pv_abort("`data` must be a data frame.")
  }
  materialized <- NULL
  if (!is.null(model_bundle)) {
    pv_binding_validate_model_bundle(model_bundle, data, formula)
    materialized <- pv_stack_materialized_design(model_bundle)
  }
  reserved_names <- c(
    pv_stack_reserved_cols(),
    materialized$column_names,
    materialized$offset_name
  )
  reserved <- intersect(reserved_names, names(data))
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
    if (!is.null(materialized)) {
      for (index in seq_along(materialized$column_names)) {
        part[[materialized$column_names[[index]]]] <-
          unname(model_bundle$model_matrix[, index])
      }
      if (!is.null(materialized$offset_name)) {
        part[[materialized$offset_name]] <- model_bundle$offset
      }
    }
    long_parts[[m]] <- part
  }
  long_data <- do.call(rbind, long_parts)
  rownames(long_data) <- NULL

  backend_rhs <- if (is.null(materialized)) {
    rhs
  } else {
    terms <- paste(materialized$column_names, collapse = " + ")
    if (!is.null(materialized$offset_name)) {
      terms <- paste0(terms, " + offset(", materialized$offset_name, ")")
    }
    paste("0 +", terms)
  }
  long_formula <- stats::as.formula(
    sprintf(".pvstackr_y | weights(.pvstackr_weight) ~ %s", backend_rhs),
    env = if (is.null(materialized)) environment(formula) else asNamespace("stats")
  )
  per_pv_weight_sum <- vapply(pv_cols, function(pv) {
    sum(long_data$.pvstackr_weight[long_data$.pvstackr_pv == pv])
  }, numeric(1))

  hash_columns <- c(
    pv_stack_reserved_cols(),
    materialized$column_names,
    materialized$offset_name
  )
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
      long_data_hash = pv_hash_stack_data(long_data, hash_columns),
      long_data_hash_columns = hash_columns,
      model_matrix_materialized = !is.null(materialized),
      model_matrix_bundle_hash = if (is.null(materialized)) NULL else model_bundle$bundle_hash,
      model_matrix_values_hash = if (is.null(materialized)) NULL else model_bundle$model_matrix_values_hash,
      offset_values_hash = if (is.null(materialized)) NULL else model_bundle$offset_values_hash
    ),
    draw_name_map = if (is.null(materialized)) NULL else stats::setNames(
      materialized$reportable_fe_names,
      materialized$backend_fe_names
    )
  )
}

pv_stack_cache_spec <- function(cache_dir = "cache", cache_stem = "pvstackr-stack",
                                package_managed = FALSE) {
  package_managed <- pv_assert_scalar_logical(package_managed, "package_managed")
  if (is.null(cache_dir)) {
    return(list(
      enabled = FALSE,
      package_managed = FALSE,
      file = NULL,
      file_refit = "never",
      cache_stem = NA_character_,
      directory_created = FALSE,
      writable = NA
    ))
  }
  cache_dir <- pv_assert_scalar_string(cache_dir, "cache_dir")
  cache_stem <- pv_assert_scalar_string(cache_stem, "cache_stem")
  if (!nzchar(cache_dir)) {
    pv_abort("`cache_dir` must be `NULL` or a non-empty directory path.")
  }
  if (!nzchar(cache_stem) || cache_stem %in% c(".", "..") ||
      !identical(basename(cache_stem), cache_stem) ||
      grepl("[/\\\\]", cache_stem)) {
    pv_abort("`cache_stem` must be a single safe file stem without path separators or traversal.")
  }

  directory_created <- FALSE
  writable <- NA
  resolved_dir <- path.expand(cache_dir)
  if (package_managed) {
    if (file.exists(resolved_dir) && !dir.exists(resolved_dir)) {
      pv_abort("`cache_dir` exists but is not a directory.")
    }
    if (!dir.exists(resolved_dir)) {
      created <- dir.create(
        resolved_dir,
        recursive = TRUE,
        showWarnings = FALSE,
        mode = "0700"
      )
      if (!isTRUE(created) || !dir.exists(resolved_dir)) {
        pv_abort("Could not create `cache_dir` for the bundled brms backend.")
      }
      directory_created <- TRUE
    }
    resolved_dir <- normalizePath(resolved_dir, winslash = "/", mustWork = TRUE)
    cache_file <- file.path(resolved_dir, cache_stem)
    cache_targets <- c(cache_file, paste0(cache_file, ".rds"))
    if (any(dir.exists(cache_targets))) {
      pv_abort("The bundled brms cache target collides with an existing directory.")
    }
    probe <- tempfile(".pvstackr-write-probe-", tmpdir = resolved_dir)
    writable <- isTRUE(suppressWarnings(file.create(probe)))
    if (file.exists(probe)) {
      unlink(probe, force = TRUE)
    }
    if (!writable) {
      pv_abort("`cache_dir` is not writable for the bundled brms backend.")
    }
  }

  list(
    enabled = TRUE,
    package_managed = package_managed,
    file = file.path(resolved_dir, cache_stem),
    file_refit = "on_change",
    cache_stem = cache_stem,
    directory_created = directory_created,
    writable = writable
  )
}

pv_stack_cache_provenance <- function(cache_spec) {
  list(
    enabled = isTRUE(cache_spec$enabled),
    policy = if (!isTRUE(cache_spec$enabled)) {
      "disabled"
    } else if (isTRUE(cache_spec$package_managed)) {
      "bundled_brms_managed"
    } else {
      "injected_adapter_managed"
    },
    cache_stem = cache_spec$cache_stem,
    file_refit = cache_spec$file_refit,
    directory_created = isTRUE(cache_spec$directory_created),
    writable = cache_spec$writable
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
  cache_spec = NULL,
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
  if (is.null(cache_spec)) {
    cache_spec <- pv_stack_cache_spec(
      cache_dir = cache_dir,
      cache_stem = cache_stem,
      package_managed = FALSE
    )
  }
  cache_required <- c(
    "enabled", "package_managed", "file", "file_refit", "cache_stem",
    "directory_created", "writable"
  )
  if (!is.list(cache_spec) || !all(cache_required %in% names(cache_spec))) {
    pv_abort("Internal cache specification is incomplete.")
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
      file = cache_spec$file,
      file_refit = cache_spec$file_refit
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

pv_stack_redact_dropped_names <- function(dropped_names) {
  dropped_names <- unname(as.character(dropped_names))
  if (length(dropped_names) == 0L) {
    return(character())
  }
  sampler_names <- c(
    "lp__", "accept_stat__", "stepsize__", "treedepth__",
    "n_leapfrog__", "divergent__", "energy__"
  )
  canonical_redaction <- grepl(
    "^<redacted_draw_column_[0-9]{3}>$",
    dropped_names
  )
  model_parameter <- grepl(
    "^(b_|sd_|cor_)[A-Za-z0-9_.\\[\\],:()\u007c-]+$",
    dropped_names,
    perl = TRUE
  ) | dropped_names %in% c("sigma", "tau")
  safe <- dropped_names %in% sampler_names | model_parameter |
    canonical_redaction
  redacted <- sprintf("<redacted_draw_column_%03d>", seq_along(dropped_names))
  dropped_names[!safe] <- redacted[!safe]
  dropped_names
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
  vc_names <- all_names[vc_idx]
  vc_names_valid <- vc_names %in% c("sigma", "tau") | grepl(
    "^(sd_|cor_)[A-Za-z0-9_.:()\u007c-]+$",
    vc_names,
    perl = TRUE
  )
  if (any(!vc_names_valid)) {
    pv_abort(
      "Variance-component selections may contain only `sd_`, `cor_`, `sigma`, or `tau` nuisance parameters; level-specific `r_` draws are not retained."
    )
  }
  dropped_idx <- setdiff(seq_along(all_names), keep_idx)
  dropped_names <- pv_stack_redact_dropped_names(all_names[dropped_idx])
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
  matrix(
    as.numeric(log_lik),
    nrow = n_draws,
    ncol = n_obs,
    dimnames = list(NULL, NULL)
  )
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

pv_stack_materialized_prior <- function(prior, draw_name_map) {
  if (is.null(draw_name_map) || is.null(prior)) {
    return(prior)
  }
  if (!is.data.frame(prior) ||
      !all(c("class", "coef") %in% names(prior)) ||
      !is.character(prior$class) || !is.character(prior$coef) ||
      anyNA(prior$class) || anyNA(prior$coef)) {
    pv_abort(
      paste0(
        "Materialized stack_direct requires NULL or a recognized invariant ",
        "prior table; opaque and coefficient-specific prior objects are unsupported."
      )
    )
  }
  scope_fields <- intersect(c("resp", "dpar", "nlpar", "group"), names(prior))
  scoped <- if (length(scope_fields) == 0L) {
    rep(FALSE, nrow(prior))
  } else {
    vapply(seq_len(nrow(prior)), function(index) {
      values <- unlist(prior[index, scope_fields, drop = FALSE], use.names = FALSE)
      any(is.na(values) | nzchar(as.character(values)))
    }, logical(1))
  }
  coefficient_specific <-
    prior$class == "Intercept" |
    nzchar(prior$coef) |
    scoped
  has_intercept <- "b_Intercept" %in% unname(draw_name_map)
  global_b_changes_scope <- prior$class == "b" & !nzchar(prior$coef) &
    has_intercept
  allowed_class <- prior$class %in% c("b", "sigma")
  if (any(coefficient_specific | global_b_changes_scope | !allowed_class)) {
    pv_abort(
      paste0(
        "Coefficient-, intercept-, scoped-, or unrecognized priors cannot be ",
        "preserved exactly after stack_direct model-matrix materialization."
      )
    )
  }
  prior
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

pv_stack_direct_preflight_fields <- function() {
  c(
    "formula", "formula_string", "rhs_string", "fe_names", "target_hash",
    "target_source", "policy", "binding_proof"
  )
}

pv_stack_direct_preflight <- function(
  data,
  formula,
  target,
  family = NULL,
  return_model_bundle = FALSE
) {
  return_model_bundle <- pv_assert_scalar_logical(
    return_model_bundle,
    "return_model_bundle"
  )
  if (!is.data.frame(data)) {
    pv_abort("`data` must be a data frame.")
  }
  if (!is.list(target)) {
    validate_pvstackr_brr_target(target)
  }
  target_schema <- target$schema_version
  if (!identical(target_schema, "0.1.0") &&
      !identical(target_schema, "0.2.0")) {
    validate_pvstackr_brr_target(target)
  }
  unbound_target <- identical(target_schema, "0.1.0") ||
    !all(c("binding_manifest", "target_content") %in% names(target))
  if (unbound_target) {
    pv_binding_abort(
      "PV_BIND_E080",
      "Public stack_direct fitting requires a schema-0.2 target with owned binding manifest and target content.",
      "target_content"
    )
  }
  validate_pvstackr_brr_target(target)
  pv_binding_manifest_assert_reportable(target$binding_manifest)
  pv_binding_target_manifest_validate(
    target$target_content,
    target$binding_manifest
  )
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
  family_link <- pv_compatibility_family_link(family)
  model_bundle <- pv_binding_resolve_model_bundle(
    data = data,
    formula = formula
  )
  fe_names <- pv_normalize_fe_names(colnames(model_bundle$model_matrix))
  current_estimand <- pv_compatibility_brr_fay_estimand(
    fe_names = fe_names,
    interval_role = target$interval_role,
    coverage_claim_allowed = target$coverage_claim_allowed,
    estimand_id = target$binding_manifest$components$estimand$estimand_id,
    target_source = target$target_source,
    target_engine_id = target$target_content$target_engine_id,
    parameter_scope = target$binding_manifest$components$estimand$parameter_scope
  )
  estimand_fields <- c(
    "estimand_id", "target_source", "target_engine_id", "parameter_scope",
    "fe_names", "interval_role", "coverage_claim_allowed"
  )
  estimand_metadata <- current_estimand[estimand_fields]
  current_manifest <- pv_binding_manifest_build(
    data = data,
    formula = formula,
    pv_cols = target$pv_cols,
    weight_col = target$weight_col,
    rep_weight_cols = target$rep_weight_cols,
    fay_k = target$fay_k,
    id_cols = target$id_cols,
    family_link = family_link,
    estimand_contrast = NULL,
    estimand_metadata = estimand_metadata,
    model_bundle = model_bundle
  )
  current_manifest <- c(
    current_manifest,
    list(model_bundle_hash = model_bundle$bundle_hash)
  )
  # Migration identity is authenticated target provenance, not current-data
  # evidence. Reattach only the strict registered projection after every raw
  # component has been recomputed so migrated targets can produce the same
  # canonical manifest hash without reusing stored component hashes.
  if ("migration" %in% names(target$binding_manifest)) {
    current_manifest <- c(
      current_manifest,
      list(migration = target$binding_manifest$migration)
    )
  }
  current_manifest$manifest_hash <- pv_binding_hash_payload(
    pv_binding_manifest_hash_payload(current_manifest),
    "manifest"
  )
  pv_binding_manifest_validate(current_manifest)
  binding_proof <- pv_binding_proof_build(
    target_manifest = target$binding_manifest,
    current_manifest = current_manifest
  )
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
    ),
    binding_proof = binding_proof
  )
  if (!identical(names(out), pv_stack_direct_preflight_fields())) {
    pv_abort("Internal stack_direct preflight fields are not canonical.")
  }
  class(out) <- c("pvstackr_stack_direct_preflight", "list")
  if (isTRUE(return_model_bundle)) {
    return(list(preflight = out, model_bundle = model_bundle))
  }
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
  additional_args = list(),
  resolved_model_bundle = NULL,
  resolved_binding_manifest = NULL
) {
  control <- pv_validate_control(control)
  bundled_backend <- FALSE
  if (is.null(fit_function)) {
    if (identical(control$backend, "brms")) {
      bundled_backend <- TRUE
      fit_function <- pv_backend_brms_fit_function
      if (is.null(draws_function)) {
        draws_function <- pv_backend_brms_draws_function
      }
    } else {
      pv_abort("`fit_function` is required unless `pv_control(backend = \"brms\")` selects the bundled brms backend.")
    }
  }
  if (is.null(draws_function)) {
    pv_abort("`draws_function` is required when using an injected `fit_function`.")
  }
  extract_log_lik <- pv_assert_scalar_logical(extract_log_lik, "extract_log_lik")

  engine_spec <- if (bundled_backend) {
    pv_backend_brms_engine_spec()
  } else {
    pv_backend_injected_engine_spec(control$backend)
  }
  cache_spec <- pv_stack_cache_spec(
    cache_dir = cache_dir,
    cache_stem = cache_stem,
    package_managed = bundled_backend
  )
  cache_provenance <- pv_stack_cache_provenance(cache_spec)

  prepared <- pv_prepare_stack_data(
    data = data,
    formula = formula,
    pv_cols = pv_cols,
    weight_col = weight_col,
    model_bundle = resolved_model_bundle
  )
  if (!is.null(resolved_binding_manifest)) {
    pv_binding_manifest_validate(resolved_binding_manifest)
    if (is.null(resolved_model_bundle) ||
        !identical(
          resolved_model_bundle$bundle_hash,
          resolved_binding_manifest$model_bundle_hash
        )) {
      pv_abort(
        "Authenticated long-data provenance requires the exact resolved model bundle."
      )
    }
    prepared$weight_summary$long_data_hash <-
      pv_binding_stack_long_data_hash(resolved_binding_manifest)
  }
  backend_prior <- pv_stack_materialized_prior(
    prior,
    prepared$draw_name_map
  )
  fit_args <- pv_stack_build_fit_args(
    prepared = prepared,
    family = family,
    prior = backend_prior,
    chains = control$chains,
    iter = control$iter,
    warmup = control$warmup,
    cores = control$cores,
    seed = control$seed,
    backend = if (bundled_backend) engine_spec$resolved_backend else control$backend,
    cache_dir = cache_dir,
    cache_stem = cache_stem,
    cache_spec = cache_spec,
    additional_args = additional_args
  )
  fit <- do.call(fit_function, fit_args)
  draws <- as.matrix(draws_function(fit))
  if (!is.null(prepared$draw_name_map)) {
    backend_names <- names(prepared$draw_name_map)
    reportable_names <- unname(prepared$draw_name_map)
    if (!is.null(param_map)) {
      map <- pv_stack_param_map(draws, param_map = param_map)
      if (length(map$fe_idx) != length(reportable_names)) {
        pv_abort("Explicit param_map fixed effects must match the complete bound design width and order.")
      }
      nonselected <- setdiff(seq_len(ncol(draws)), map$original_fe_idx)
      if (any(colnames(draws)[nonselected] %in% c(backend_names, reportable_names))) {
        pv_abort("Explicit param_map leaves conflicting materialized or reportable fixed-effect columns unselected.")
      }
      colnames(map$draws_selected)[map$fe_idx] <- reportable_names
      map$fe_names <- reportable_names
    } else {
      present_backend <- backend_names %in% colnames(draws)
      present_reportable <- reportable_names %in% colnames(draws)
      if (any(present_backend)) {
        if (!all(present_backend)) {
          pv_abort("Materialized backend draws contain only part of the fixed-effect design block.")
        }
        conflicts <- reportable_names %in% setdiff(colnames(draws), backend_names)
        if (any(conflicts)) {
          pv_abort("Materialized backend draw names conflict with reportable fixed-effect names.")
        }
        colnames(draws)[match(backend_names, colnames(draws))] <- reportable_names
      } else if (!all(present_reportable)) {
        pv_abort("Backend draws do not expose the complete bound fixed-effect design block.")
      }
      map <- pv_stack_param_map(draws, param_map = NULL)
    }
  } else {
    map <- pv_stack_param_map(draws, param_map = param_map)
  }
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
  diagnostics <- pv_stack_sampler_diagnostics(
    fit = fit,
    bundled_backend = bundled_backend,
    diagnose_function = diagnose_function,
    chains = control$chains,
    iter = control$iter,
    warmup = control$warmup
  )
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
  retained_formula <- if (!is.null(prepared$draw_name_map)) {
    stats::as.formula(prepared$formula_string, env = asNamespace("stats"))
  } else {
    stats::as.formula(prepared$formula_string, env = baseenv())
  }

  out <- list(
    stacked_draws = stacked_draws,
    diagnostics = diagnostics,
    log_lik = log_lik,
    psi_hat_fe = colMeans(stacked_draws[, param_map$fe_idx, drop = FALSE]),
    param_map = param_map,
    formula = retained_formula,
    formula_string = prepared$formula_string,
    weight_summary = prepared$weight_summary,
    meta = list(
      topology = "single_long_fit",
      engine_id = engine_spec$engine_id,
      fit_engine = engine_spec$adapter_id,
      adapter_source = engine_spec$adapter_source,
      resolved_backend = engine_spec$resolved_backend,
      backend_selection_reason = engine_spec$selection_reason,
      cache = cache_provenance,
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
      dropped_draw_columns = param_map$dropped_names,
      model_matrix_materialized =
        isTRUE(prepared$weight_summary$model_matrix_materialized),
      model_matrix_bundle_hash =
        prepared$weight_summary$model_matrix_bundle_hash,
      model_matrix_values_hash =
        prepared$weight_summary$model_matrix_values_hash,
      offset_values_hash = prepared$weight_summary$offset_values_hash
    ),
    fit = if (isTRUE(control$keep_backend_fit)) fit else NULL,
    prepared_data = if (isTRUE(control$keep_data)) prepared$data else NULL,
    control = control,
    schema_version = "0.2.0",
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
      model_matrix_materialized =
        isTRUE(prepared$weight_summary$model_matrix_materialized),
      model_matrix_bundle_hash =
        prepared$weight_summary$model_matrix_bundle_hash,
      model_matrix_values_hash =
        prepared$weight_summary$model_matrix_values_hash,
      offset_values_hash = prepared$weight_summary$offset_values_hash,
      backend = control$backend,
      engine = engine_spec,
      cache = cache_provenance,
      backend_fit_retained = isTRUE(control$keep_backend_fit),
      log_lik_retained = isTRUE(control$keep_log_lik),
      schema_version = "0.2.0"
    ),
    warnings = warnings
  )
  class(out) <- c("pvstackr_stack_fit", "list")
  requested_return_draws <- control$return_draws
  out$control$return_draws <- TRUE
  pv_stack_fit_draw_projection(out, requested_return_draws)
}
