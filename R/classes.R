pv_schema_version <- function() {
  "0.1.0"
}

pv_required_fields <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    pv_abort(sprintf("%s is missing required field(s): %s.", label, paste(missing, collapse = ", ")))
  }
  invisible(x)
}

pv_provenance <- function(function_name, ..., schema_version = pv_schema_version()) {
  function_name <- pv_assert_scalar_string(function_name, "function_name")
  schema_version <- pv_validate_schema_version(schema_version)
  extra <- list(...)
  if (length(extra) > 0L && (is.null(names(extra)) || any(!nzchar(names(extra))))) {
    pv_abort("Provenance extras must be named.")
  }
  c(
    list(
      function_name = function_name,
      package = "pvstackr",
      schema_version = schema_version
    ),
    extra
  )
}

pv_validate_named_list_field <- function(x, name) {
  if (!is.list(x)) {
    pv_abort(sprintf("`%s` must be a list.", name))
  }
  if (length(x) > 0L && (is.null(names(x)) || any(!nzchar(names(x))))) {
    pv_abort(sprintf("`%s` must be a named list when non-empty.", name))
  }
  x
}

pv_validate_character_field <- function(x, name, allow_empty = TRUE) {
  if (!is.character(x) || any(is.na(x))) {
    pv_abort(sprintf("`%s` must be a character vector with no missing values.", name))
  }
  if (!allow_empty && length(x) == 0L) {
    pv_abort(sprintf("`%s` must not be empty.", name))
  }
  if (any(!nzchar(x))) {
    pv_abort(sprintf("`%s` must not contain empty values.", name))
  }
  x
}

pv_validate_schema_version <- function(x) {
  pv_assert_scalar_string(x, "schema_version")
}

pv_validate_design_formula <- function(formula) {
  rhs <- pv_formula_rhs_checked(formula)
  if (pv_formula_has_weights_call(rhs)) {
    pv_abort("Do not embed `weights()` in `formula`; pass weight columns explicitly.")
  }
  formula
}

pv_design_hashes <- function(data, pv_cols, weight_col, rep_weight_cols, id_cols, fay_k) {
  row_support <- if (length(id_cols) > 0L) {
    data[id_cols]
  } else {
    list(n = nrow(data), row = seq_len(nrow(data)))
  }
  weight_payload <- list(
    weight_col = weight_col,
    weight = if (is.null(weight_col)) NULL else data[weight_col],
    rep_weight_cols = rep_weight_cols,
    rep_weights = if (length(rep_weight_cols) == 0L) NULL else data[rep_weight_cols],
    fay_k = fay_k
  )
  list(
    row_support_hash = pv_hash_payload(row_support),
    pv_value_hash = pv_hash_payload(data[pv_cols]),
    weight_design_hash = pv_hash_payload(weight_payload)
  )
}

pv_design_manifest <- function(data, formula, pv_cols, weight_col, rep_weight_cols, id_cols, fay_k) {
  list(
    n = nrow(data),
    columns = names(data),
    column_classes = vapply(data, function(x) paste(class(x), collapse = "/"), character(1)),
    formula_string = paste(deparse(formula, width.cutoff = 500L), collapse = ""),
    pv_cols = pv_cols,
    weight_col = weight_col,
    rep_weight_cols = rep_weight_cols,
    id_cols = id_cols,
    fay_k = fay_k
  )
}

new_pvstackr_design <- function(
  data,
  formula,
  pv_cols,
  weight_col = NULL,
  rep_weight_cols = NULL,
  fay_k = 0.5,
  id_cols = NULL,
  roles = list(),
  provenance = list(),
  warnings = character()
) {
  if (!is.data.frame(data)) {
    pv_abort("`data` must be a data frame.")
  }
  formula <- pv_validate_design_formula(formula)
  pv_cols <- pv_validate_pv_columns(data, pv_cols)

  if (!is.null(weight_col)) {
    weight_col <- pv_assert_scalar_string(weight_col, "weight_col")
    pv_validate_columns(data, weight_col, "weight_col")
    pv_validate_weight_vector(data[[weight_col]], weight_col, nrow(data))
  }

  if (is.null(rep_weight_cols)) {
    rep_weight_cols <- character()
  } else {
    rep_weight_cols <- pv_validate_columns(data, rep_weight_cols, "rep_weight_cols")
    if (length(rep_weight_cols) < 2L) {
      pv_abort("`rep_weight_cols`, when supplied, must contain at least two BRR replicate weights.")
    }
    for (col in rep_weight_cols) {
      pv_validate_weight_vector(data[[col]], col, nrow(data))
    }
  }
  if (!is.null(weight_col) && anyDuplicated(c(weight_col, rep_weight_cols))) {
    pv_abort("`weight_col` and `rep_weight_cols` must name distinct columns.")
  }

  fay_k <- pv_validate_fay_k(fay_k)
  id_cols <- pv_validate_id_columns(data, id_cols)
  roles <- pv_validate_named_list_field(roles, "roles")
  provenance <- pv_validate_named_list_field(provenance, "provenance")
  warnings <- pv_validate_character_field(warnings, "warnings")
  hashes <- pv_design_hashes(data, pv_cols, weight_col, rep_weight_cols, id_cols, fay_k)
  data_manifest <- pv_design_manifest(data, formula, pv_cols, weight_col, rep_weight_cols, id_cols, fay_k)
  design_hash <- pv_hash_payload(c(data_manifest, hashes))

  design <- list(
    data = data,
    data_manifest = data_manifest,
    formula = formula,
    formula_string = data_manifest$formula_string,
    pv_cols = pv_cols,
    weight_col = weight_col,
    rep_weight_cols = rep_weight_cols,
    fay_k = fay_k,
    id_cols = id_cols,
    M = length(pv_cols),
    R = length(rep_weight_cols),
    n = nrow(data),
    row_support_hash = hashes$row_support_hash,
    row_support = list(
      type = if (length(id_cols) > 0L) "id_cols" else "row_number",
      id_cols = id_cols,
      n = nrow(data),
      hash = hashes$row_support_hash
    ),
    pv_value_hash = hashes$pv_value_hash,
    weight_design_hash = hashes$weight_design_hash,
    design_hash = design_hash,
    roles = roles,
    created_at = as.character(Sys.time()),
    schema_version = pv_schema_version(),
    provenance = c(
      pv_provenance("new_pvstackr_design"),
      provenance
    ),
    warnings = warnings
  )
  class(design) <- c("pvstackr_design", "list")
  validate_pvstackr_design(design)
  design
}

validate_pvstackr_design <- function(design) {
  pv_assert_named_list(design, "design")
  required <- c(
    "data", "data_manifest", "formula", "formula_string", "pv_cols",
    "weight_col", "rep_weight_cols", "fay_k", "id_cols", "M", "R", "n",
    "row_support_hash", "row_support", "pv_value_hash", "weight_design_hash",
    "design_hash", "roles", "created_at", "schema_version", "provenance",
    "warnings"
  )
  pv_required_fields(design, required, "Design object")
  if (!is.data.frame(design$data)) {
    pv_abort("Design `data` must be a data frame.")
  }
  pv_validate_design_formula(design$formula)
  pv_cols <- pv_validate_pv_columns(design$data, design$pv_cols)
  if (!identical(design$M, length(pv_cols))) {
    pv_abort("Design `M` must equal the number of plausible-value columns.")
  }
  if (!is.null(design$weight_col)) {
    pv_assert_scalar_string(design$weight_col, "weight_col")
    pv_validate_columns(design$data, design$weight_col, "weight_col")
    pv_validate_weight_vector(design$data[[design$weight_col]], design$weight_col, nrow(design$data))
  }
  if (!is.character(design$rep_weight_cols)) {
    pv_abort("Design `rep_weight_cols` must be a character vector.")
  }
  if (length(design$rep_weight_cols) > 0L) {
    rep_weight_cols <- pv_validate_columns(design$data, design$rep_weight_cols, "rep_weight_cols")
    if (length(rep_weight_cols) < 2L) {
      pv_abort("Design `rep_weight_cols` must contain at least two BRR replicate weights when present.")
    }
    for (col in rep_weight_cols) {
      pv_validate_weight_vector(design$data[[col]], col, nrow(design$data))
    }
  }
  if (!is.null(design$weight_col) && anyDuplicated(c(design$weight_col, design$rep_weight_cols))) {
    pv_abort("Design weight columns must be distinct.")
  }
  if (!identical(design$R, length(design$rep_weight_cols))) {
    pv_abort("Design `R` must equal the number of replicate-weight columns.")
  }
  if (!identical(design$n, nrow(design$data))) {
    pv_abort("Design `n` must equal the number of data rows.")
  }
  pv_validate_fay_k(design$fay_k)
  id_cols <- pv_validate_id_columns(design$data, design$id_cols)
  if (!identical(design$id_cols, id_cols)) {
    pv_abort("Design `id_cols` must be normalized character columns.")
  }
  hashes <- pv_design_hashes(
    design$data,
    design$pv_cols,
    design$weight_col,
    design$rep_weight_cols,
    design$id_cols,
    design$fay_k
  )
  data_manifest <- pv_design_manifest(
    design$data,
    design$formula,
    design$pv_cols,
    design$weight_col,
    design$rep_weight_cols,
    design$id_cols,
    design$fay_k
  )
  design_hash <- pv_hash_payload(c(data_manifest, hashes))
  if (!identical(design$formula_string, data_manifest$formula_string)) {
    pv_abort("Design `formula_string` does not match `formula`.")
  }
  if (!identical(design$data_manifest, data_manifest)) {
    pv_abort("Design `data_manifest` does not match the current object fields.")
  }
  if (!identical(design$row_support_hash, hashes$row_support_hash) ||
      !identical(design$pv_value_hash, hashes$pv_value_hash) ||
      !identical(design$weight_design_hash, hashes$weight_design_hash) ||
      !identical(design$design_hash, design_hash)) {
    pv_abort("Design hashes do not match the current object fields.")
  }
  pv_assert_named_list(design$row_support, "row_support")
  row_support_required <- c("type", "id_cols", "n", "hash")
  missing_row_support <- setdiff(row_support_required, names(design$row_support))
  if (length(missing_row_support) > 0L) {
    pv_abort(sprintf("Design `row_support` is missing required field(s): %s.", paste(missing_row_support, collapse = ", ")))
  }
  expected_row_support <- list(
    type = if (length(design$id_cols) > 0L) "id_cols" else "row_number",
    id_cols = design$id_cols,
    n = nrow(design$data),
    hash = hashes$row_support_hash
  )
  if (!identical(design$row_support, expected_row_support)) {
    pv_abort("Design `row_support` metadata must match row-support fields and hash.")
  }
  pv_assert_scalar_string(design$created_at, "created_at")
  pv_validate_schema_version(design$schema_version)
  pv_validate_named_list_field(design$roles, "roles")
  pv_validate_named_list_field(design$provenance, "provenance")
  pv_validate_character_field(design$warnings, "warnings")
  invisible(design)
}

validate_pvstackr_brr_target <- function(target) {
  pv_validate_brr_target(target)
}

new_pvstackr_brr_target <- function(target) {
  pv_assert_named_list(target, "target")
  class(target) <- unique(c("pvstackr_brr_target", "list", class(target)))
  pv_validate_brr_target(target)
  target
}

validate_pvstackr_stack_fit <- function(stack_fit) {
  pv_assert_named_list(stack_fit, "stack_fit")
  required <- c(
    "stacked_draws", "param_map", "psi_hat_fe", "formula", "weight_summary",
    "diagnostics", "log_lik", "fit", "meta", "control", "schema_version",
    "provenance", "warnings"
  )
  missing <- setdiff(required, names(stack_fit))
  if (length(missing) > 0L) {
    pv_abort(sprintf("Stack-fit object is missing required field(s): %s.", paste(missing, collapse = ", ")))
  }
  draws <- ccc_as_draw_matrix(stack_fit$stacked_draws)
  pv_assert_named_list(stack_fit$param_map, "param_map")
  map_required <- c("fe_idx", "vc_idx", "fe_names", "vc_names")
  map_missing <- setdiff(map_required, names(stack_fit$param_map))
  if (length(map_missing) > 0L) {
    pv_abort(sprintf("Stack-fit `param_map` is missing required field(s): %s.", paste(map_missing, collapse = ", ")))
  }
  fe_idx <- stack_fit$param_map$fe_idx
  vc_idx <- stack_fit$param_map$vc_idx
  if (!is.integer(fe_idx) && !is.numeric(fe_idx)) {
    pv_abort("Stack-fit `param_map$fe_idx` must be integer indexes.")
  }
  if (!is.integer(vc_idx) && !is.numeric(vc_idx)) {
    pv_abort("Stack-fit `param_map$vc_idx` must be integer indexes.")
  }
  fe_idx <- as.integer(fe_idx)
  vc_idx <- as.integer(vc_idx)
  if (length(fe_idx) == 0L || any(fe_idx < 1L | fe_idx > ncol(draws))) {
    pv_abort("Stack-fit fixed-effect indexes are invalid.")
  }
  if (length(vc_idx) > 0L && any(vc_idx < 1L | vc_idx > ncol(draws))) {
    pv_abort("Stack-fit variance-component indexes are invalid.")
  }
  if (length(intersect(fe_idx, vc_idx)) > 0L) {
    pv_abort("Stack-fit fixed-effect and variance-component indexes must be disjoint.")
  }
  if (!identical(stack_fit$param_map$fe_names, colnames(draws)[fe_idx])) {
    pv_abort("Stack-fit fixed-effect names must match selected draw columns.")
  }
  if (!identical(stack_fit$param_map$vc_names, colnames(draws)[vc_idx])) {
    pv_abort("Stack-fit variance-component names must match selected draw columns.")
  }
  if (!is.null(stack_fit$param_map$original_fe_idx)) {
    original_fe_idx <- stack_fit$param_map$original_fe_idx
    if ((!is.integer(original_fe_idx) && !is.numeric(original_fe_idx)) ||
        length(original_fe_idx) != length(fe_idx) ||
        anyNA(original_fe_idx) ||
        any(original_fe_idx < 1L) ||
        anyDuplicated(as.integer(original_fe_idx))) {
      pv_abort("Stack-fit `param_map$original_fe_idx` must identify unique original draw columns.")
    }
  }
  if (!is.null(stack_fit$param_map$original_vc_idx)) {
    original_vc_idx <- stack_fit$param_map$original_vc_idx
    if ((!is.integer(original_vc_idx) && !is.numeric(original_vc_idx)) ||
        anyNA(original_vc_idx) ||
        any(original_vc_idx < 1L) ||
        anyDuplicated(as.integer(original_vc_idx))) {
      pv_abort("Stack-fit `param_map$original_vc_idx` must identify unique original draw columns.")
    }
  }
  if (!is.null(stack_fit$param_map$original_dropped_idx)) {
    original_dropped_idx <- stack_fit$param_map$original_dropped_idx
    if ((!is.integer(original_dropped_idx) && !is.numeric(original_dropped_idx)) ||
        anyNA(original_dropped_idx) ||
        any(original_dropped_idx < 1L) ||
        anyDuplicated(as.integer(original_dropped_idx))) {
      pv_abort("Stack-fit `param_map$original_dropped_idx` must identify unique original draw columns.")
    }
  }
  if (!is.null(stack_fit$param_map$dropped_names)) {
    pv_validate_character_field(stack_fit$param_map$dropped_names, "param_map$dropped_names")
  }
  if (!is.null(stack_fit$param_map$map_source)) {
    map_source <- pv_assert_scalar_string(stack_fit$param_map$map_source, "param_map$map_source")
    if (!map_source %in% c("auto_regex", "explicit")) {
      pv_abort("Stack-fit `param_map$map_source` must be `auto_regex` or `explicit`.")
    }
  }
  if (!is.numeric(stack_fit$psi_hat_fe) || any(!is.finite(stack_fit$psi_hat_fe)) ||
      !identical(names(stack_fit$psi_hat_fe), stack_fit$param_map$fe_names)) {
    pv_abort("Stack-fit `psi_hat_fe` must be finite, named, and aligned to fixed effects.")
  }
  psi_hat_expected <- colMeans(draws[, fe_idx, drop = FALSE])
  if (max(abs(stack_fit$psi_hat_fe - psi_hat_expected)) > 1e-10) {
    pv_abort("Stack-fit `psi_hat_fe` must equal the fixed-effect draw means.")
  }
  if (!inherits(stack_fit$formula, "formula")) {
    pv_abort("Stack-fit `formula` must be a formula.")
  }
  pv_assert_named_list(stack_fit$weight_summary, "weight_summary")
  ws_required <- c("M", "n_original", "n_long", "total_long_weight", "long_data_hash")
  ws_missing <- setdiff(ws_required, names(stack_fit$weight_summary))
  if (length(ws_missing) > 0L) {
    pv_abort(sprintf("Stack-fit `weight_summary` is missing required field(s): %s.", paste(ws_missing, collapse = ", ")))
  }
  if (stack_fit$weight_summary$n_long != stack_fit$weight_summary$n_original * stack_fit$weight_summary$M) {
    pv_abort("Stack-fit long-data dimensions are inconsistent with `M` and original row count.")
  }
  if (abs(stack_fit$weight_summary$total_long_weight - stack_fit$weight_summary$n_original) > 1e-8) {
    pv_abort("Stack-fit total long weight must equal the original row count.")
  }
  if (!is.null(stack_fit$log_lik)) {
    pv_validate_log_lik(stack_fit$log_lik, nrow(draws), stack_fit$weight_summary$n_long)
  }
  if (!is.null(stack_fit$meta$long_data_hash) &&
      !identical(stack_fit$meta$long_data_hash, stack_fit$weight_summary$long_data_hash)) {
    pv_abort("Stack-fit `meta$long_data_hash` must match `weight_summary$long_data_hash`.")
  }
  pv_validate_named_list_field(stack_fit$diagnostics, "diagnostics")
  pv_validate_named_list_field(stack_fit$meta, "meta")
  for (prior_field in c("prior_policy", "prior_diagnostic")) {
    prior_policy <- stack_fit$meta[[prior_field]]
    if (is.null(prior_policy)) {
      next
    }
    pv_assert_named_list(prior_policy, paste0("meta$", prior_field))
    legacy_fields <- intersect(names(prior_policy), c("non_flat_prior_warning", "warn_nonflat_prior"))
    if (length(legacy_fields) > 0L) {
      pv_abort(sprintf(
        "Stack-fit `meta$%s` uses legacy prior diagnostic field(s): %s. Use `explicit_prior_warning` and `warn_explicit_prior`.",
        prior_field,
        paste(legacy_fields, collapse = ", ")
      ))
    }
    if (isTRUE(prior_policy$explicit_prior_warning) ||
        isTRUE(prior_policy$warn_explicit_prior)) {
      reason_code <- prior_policy$reason_code
      if (!identical(reason_code, "explicit_prior_warning")) {
        pv_abort(sprintf(
          "Stack-fit `meta$%s` explicit-prior diagnostics must use `reason_code = \"explicit_prior_warning\"`.",
          prior_field
        ))
      }
    }
  }
  pv_validate_control(stack_fit$control)
  pv_validate_schema_version(stack_fit$schema_version)
  pv_validate_named_list_field(stack_fit$provenance, "provenance")
  if (!is.null(stack_fit$provenance$long_data_hash) &&
      !identical(stack_fit$provenance$long_data_hash, stack_fit$weight_summary$long_data_hash)) {
    pv_abort("Stack-fit provenance hash must match `weight_summary$long_data_hash`.")
  }
  if (!is.null(stack_fit$provenance$schema_version) &&
      !identical(stack_fit$provenance$schema_version, stack_fit$schema_version)) {
    pv_abort("Stack-fit provenance schema version must match object `schema_version`.")
  }
  if (isTRUE(stack_fit$control$keep_backend_fit) && is.null(stack_fit$fit)) {
    pv_abort("Stack-fit retention policy requested backend fit retention, but `fit` is NULL.")
  }
  if (!isTRUE(stack_fit$control$keep_backend_fit) && !is.null(stack_fit$fit)) {
    pv_abort("Stack-fit retained `fit` despite `keep_backend_fit = FALSE`.")
  }
  if (isTRUE(stack_fit$control$keep_log_lik) && isTRUE(stack_fit$meta$log_lik_extracted) && is.null(stack_fit$log_lik)) {
    pv_abort("Stack-fit retention policy requested log-likelihood retention, but `log_lik` is NULL.")
  }
  if (!isTRUE(stack_fit$control$keep_log_lik) && !is.null(stack_fit$log_lik)) {
    pv_abort("Stack-fit retained `log_lik` despite `keep_log_lik = FALSE`.")
  }
  pv_validate_character_field(stack_fit$warnings, "warnings")
  invisible(stack_fit)
}

validate_pvstackr_ccc <- function(ccc) {
  pv_assert_named_list(ccc, "ccc")
  required <- c(
    "draws_calibrated", "draws_fe_cal", "A", "A_full", "psi_hat",
    "psi_raw", "psi_target", "Sigma_raw", "Sigma_target", "Sigma_cal_emp",
    "diagnostics", "flags", "param_map", "control", "target_source",
    "target_hash", "center", "policy", "ccc_status", "schema_version",
    "provenance", "warnings"
  )
  missing <- setdiff(required, names(ccc))
  if (length(missing) > 0L) {
    pv_abort(sprintf("CCC object is missing required field(s): %s.", paste(missing, collapse = ", ")))
  }
  draws_all <- ccc_as_draw_matrix(ccc$draws_calibrated)
  draws_fe <- ccc_as_draw_matrix(ccc$draws_fe_cal)
  pv_assert_named_list(ccc$param_map, "param_map")
  fe_idx <- as.integer(ccc$param_map$fe_idx)
  vc_idx <- as.integer(ccc$param_map$vc_idx)
  if (length(fe_idx) == 0L || any(fe_idx < 1L | fe_idx > ncol(draws_all))) {
    pv_abort("CCC fixed-effect indexes are invalid.")
  }
  if (length(vc_idx) > 0L && any(vc_idx < 1L | vc_idx > ncol(draws_all))) {
    pv_abort("CCC variance-component indexes are invalid.")
  }
  if (length(intersect(fe_idx, vc_idx)) > 0L) {
    pv_abort("CCC fixed-effect and variance-component indexes must be disjoint.")
  }
  if (!identical(draws_all[, fe_idx, drop = FALSE], draws_fe)) {
    pv_abort("CCC `draws_fe_cal` must equal the fixed-effect block of `draws_calibrated`.")
  }
  if (!is.numeric(ccc$psi_hat) || any(!is.finite(ccc$psi_hat)) ||
      !identical(names(ccc$psi_hat), colnames(draws_fe))) {
    pv_abort("CCC `psi_hat` must be finite, named, and aligned to fixed-effect draws.")
  }
  if (!is.numeric(ccc$psi_raw) || any(!is.finite(ccc$psi_raw)) ||
      !identical(names(ccc$psi_raw), colnames(draws_fe))) {
    pv_abort("CCC `psi_raw` must be finite, named, and aligned to fixed-effect draws.")
  }
  if (!is.numeric(ccc$psi_target) || any(!is.finite(ccc$psi_target)) ||
      !identical(names(ccc$psi_target), colnames(draws_fe))) {
    pv_abort("CCC `psi_target` must be finite, named, and aligned to fixed-effect draws.")
  }
  if (!is.matrix(ccc$A) || !is.numeric(ccc$A) || any(!is.finite(ccc$A)) ||
      !identical(dim(ccc$A), c(ncol(draws_fe), ncol(draws_fe)))) {
    pv_abort("CCC `A` must be a finite square matrix aligned to fixed-effect draws.")
  }
  if (!is.matrix(ccc$A_full) || !identical(dim(ccc$A_full), c(ncol(draws_all), ncol(draws_all)))) {
    pv_abort("CCC `A_full` dimensions must match all calibrated draw columns.")
  }
  if (max(abs(ccc$A_full[fe_idx, fe_idx, drop = FALSE] - ccc$A)) > 1e-12) {
    pv_abort("CCC `A_full` fixed-effect block must equal `A`.")
  }
  outside_idx <- setdiff(seq_len(ncol(draws_all)), fe_idx)
  if (length(outside_idx) > 0L) {
    identity_outside <- diag(ncol(draws_all))[outside_idx, outside_idx, drop = FALSE]
    if (max(abs(ccc$A_full[outside_idx, outside_idx, drop = FALSE] - identity_outside)) > 1e-12 ||
        max(abs(ccc$A_full[outside_idx, fe_idx, drop = FALSE])) > 1e-12 ||
        max(abs(ccc$A_full[fe_idx, outside_idx, drop = FALSE])) > 1e-12) {
      pv_abort("CCC `A_full` must be identity outside fixed-effect rows and columns.")
    }
  }
  pv_validate_target_matrix(ccc$Sigma_target, colnames(draws_fe))
  pv_assert_named_list(ccc$diagnostics, "diagnostics")
  diagnostic_required <- c(
    "delta_c_rel", "delta_c_max", "a_matrix_fro_rel", "kappa_A",
    "center_status", "center_separation", "conditioning_status",
    "conditioning"
  )
  diagnostic_missing <- setdiff(diagnostic_required, names(ccc$diagnostics))
  if (length(diagnostic_missing) > 0L) {
    pv_abort(sprintf("CCC diagnostics are missing required field(s): %s.", paste(diagnostic_missing, collapse = ", ")))
  }
  if (!is.numeric(ccc$diagnostics$delta_c_rel) || length(ccc$diagnostics$delta_c_rel) != 1L ||
      !is.finite(ccc$diagnostics$delta_c_rel) || ccc$diagnostics$delta_c_rel < 0) {
    pv_abort("CCC diagnostic `delta_c_rel` must be a non-negative finite scalar.")
  }
  if (!is.numeric(ccc$diagnostics$delta_c_max) || length(ccc$diagnostics$delta_c_max) != 1L ||
      !is.finite(ccc$diagnostics$delta_c_max) || ccc$diagnostics$delta_c_max < 0) {
    pv_abort("CCC diagnostic `delta_c_max` must be a non-negative finite scalar.")
  }
  if (!is.numeric(ccc$diagnostics$a_matrix_fro_rel) || length(ccc$diagnostics$a_matrix_fro_rel) != 1L ||
      !is.finite(ccc$diagnostics$a_matrix_fro_rel) || ccc$diagnostics$a_matrix_fro_rel < 0) {
    pv_abort("CCC diagnostic `a_matrix_fro_rel` must be a non-negative finite scalar.")
  }
  if (!is.numeric(ccc$diagnostics$kappa_A) || length(ccc$diagnostics$kappa_A) != 1L ||
      !is.finite(ccc$diagnostics$kappa_A) || ccc$diagnostics$kappa_A < 1) {
    pv_abort("CCC diagnostic `kappa_A` must be a finite scalar greater than or equal to 1.")
  }
  if (!identical(ccc$diagnostics$center_status, "ok") &&
      !identical(ccc$diagnostics$center_status, "warning") &&
      !identical(ccc$diagnostics$center_status, "blocked")) {
    pv_abort("CCC diagnostic `center_status` must be one of: ok, warning, blocked.")
  }
  if (!identical(ccc$diagnostics$conditioning_status, "ok") &&
      !identical(ccc$diagnostics$conditioning_status, "warning") &&
      !identical(ccc$diagnostics$conditioning_status, "blocked")) {
    pv_abort("CCC diagnostic `conditioning_status` must be one of: ok, warning, blocked.")
  }
  pv_assert_named_list(ccc$flags, "flags")
  pv_assert_named_list(ccc$control, "control")
  if (!identical(ccc$control$allow_target_nearpd, FALSE)) {
    pv_abort("CCC control `allow_target_nearpd` must be FALSE because target repair is not currently supported.")
  }
  ccc_center <- pv_validate_center(ccc$center)
  ccc_control_center <- pv_validate_center(ccc$control$center)
  if (!identical(ccc_center, ccc_control_center)) {
    pv_abort("CCC `center` and `control$center` must match.")
  }
  if (!identical(ccc$target_source, "external_brr_fay_rubin")) {
    pv_abort("CCC `target_source` must be `external_brr_fay_rubin`.")
  }
  if (!identical(ccc$ccc_status, "ok")) {
    pv_abort("CCC `ccc_status` must be `ok` for reportable calibration objects.")
  }
  if (!identical(ccc$policy$target_repair, "forbidden")) {
    pv_abort("CCC target repair policy must be forbidden.")
  }
  if (!identical(ccc$flags$target_repaired, FALSE) ||
      !identical(ccc$flags$nearpd_target, FALSE)) {
    pv_abort("CCC target repair and nearPD target flags must be FALSE when target repair is forbidden.")
  }
  if (!identical(ccc$flags$center_separation_warning, identical(ccc$diagnostics$center_status, "warning")) ||
      !identical(ccc$flags$center_separation_blocked, identical(ccc$diagnostics$center_status, "blocked"))) {
    pv_abort("CCC center separation flags must match diagnostic `center_status`.")
  }
  if (!identical(ccc$flags$conditioning_warning, identical(ccc$diagnostics$conditioning_status, "warning")) ||
      !identical(ccc$flags$conditioning_blocked, identical(ccc$diagnostics$conditioning_status, "blocked")) ||
      !identical(ccc$flags$kappa_a_warning, identical(ccc$diagnostics$conditioning_status, "warning")) ||
      !identical(ccc$flags$kappa_a_blocked, identical(ccc$diagnostics$conditioning_status, "blocked"))) {
    pv_abort("CCC conditioning flags must match diagnostic `conditioning_status`.")
  }
  pv_validate_schema_version(ccc$schema_version)
  pv_validate_named_list_field(ccc$provenance, "provenance")
  pv_validate_character_field(ccc$warnings, "warnings")
  invisible(ccc)
}

pv_validate_fit_reason_codes <- function(status, reason_codes) {
  reason_codes <- pv_validate_character_field(reason_codes, "reason_codes")
  if (identical(status, "ok") && length(reason_codes) > 0L) {
    pv_abort("OK fit objects must not include `reason_codes`.")
  }
  if (identical(status, "warning") && length(reason_codes) == 0L) {
    pv_abort("Warning fit objects must include at least one `reason_codes` value.")
  }
  if (identical(status, "blocked") && length(reason_codes) == 0L) {
    pv_abort("Blocked fit objects must include at least one `reason_codes` value.")
  }
  reason_codes
}

pv_validate_fit_warnings <- function(status, warnings) {
  warnings <- pv_validate_character_field(warnings, "warnings")
  if (identical(status, "ok") && length(warnings) > 0L) {
    pv_abort("OK fit objects must not include `warnings`.")
  }
  if (identical(status, "warning") && length(warnings) == 0L) {
    pv_abort("Warning fit objects must include at least one `warnings` value.")
  }
  warnings
}

pv_validate_method_implemented_status <- function(method, status) {
  invisible(TRUE)
}

validate_pvstackr_fit_target <- function(target) {
  if (inherits(target, "pvstackr_reference_pool")) {
    validate_pvstackr_reference_pool(target)
  } else {
    validate_pvstackr_brr_target(target)
  }
}

pv_validate_stack_direct_reportable_fit <- function(fit) {
  if (!identical(fit$method, "stack_direct")) {
    return(invisible(fit))
  }
  if (identical(fit$status, "blocked")) {
    if (is.data.frame(fit$estimates) && nrow(fit$estimates) > 0L) {
      pv_abort("Blocked stack_direct fit must not include reportable estimates.")
    }
    return(invisible(fit))
  }
  if (!fit$status %in% c("ok", "warning")) {
    return(invisible(fit))
  }

  required_non_null <- c("design", "target", "stack_fit", "ccc")
  missing <- required_non_null[vapply(required_non_null, function(name) is.null(fit[[name]]), logical(1))]
  if (length(missing) > 0L) {
    pv_abort(sprintf("Reportable stack_direct fit requires non-null component(s): %s.", paste(missing, collapse = ", ")))
  }
  if (!identical(fit$control$center, "target") ||
      !identical(fit$ccc$center, "target") ||
      !identical(fit$ccc$control$center, "target")) {
    pv_abort("Reportable stack_direct fit requires `control$center = \"target\"` and target-centered CCC output.")
  }
  pv_stack_direct_preflight(
    data = fit$design$data,
    formula = fit$design$formula,
    target = fit$target
  )
  if (!identical(fit$stack_fit$param_map$fe_names, fit$target$fe_names)) {
    pv_abort("Reportable stack_direct fit requires stack-fit fixed-effect names to match the target.")
  }
  center_status <- fit$ccc$diagnostics$center_status %||% "ok"
  if (identical(center_status, "blocked")) {
    pv_abort("Reportable stack_direct fit cannot include a blocked center separation diagnostic.")
  }
  if (identical(center_status, "warning") && identical(fit$status, "ok")) {
    pv_abort("Reportable stack_direct fit with warning-level center separation must use `status = \"warning\"`.")
  }
  center_reason_code <- fit$ccc$diagnostics$center_reason_code %||% NA_character_
  if (identical(center_status, "warning") &&
      !is.na(center_reason_code) &&
      !center_reason_code %in% fit$reason_codes) {
    pv_abort("Warning stack_direct fit must include the CCC center separation `reason_code`.")
  }
  conditioning_status <- fit$ccc$diagnostics$conditioning_status %||% "ok"
  if (identical(conditioning_status, "blocked")) {
    pv_abort("Reportable stack_direct fit cannot include a blocked CCC conditioning diagnostic.")
  }
  if (identical(conditioning_status, "warning") && identical(fit$status, "ok")) {
    pv_abort("Reportable stack_direct fit with warning-level CCC conditioning must use `status = \"warning\"`.")
  }
  conditioning_reason_code <- fit$ccc$diagnostics$conditioning_reason_code %||% NA_character_
  if (identical(conditioning_status, "warning") &&
      !is.na(conditioning_reason_code) &&
      !conditioning_reason_code %in% fit$reason_codes) {
    pv_abort("Warning stack_direct fit must include the CCC conditioning `reason_code`.")
  }
  if (!is.data.frame(fit$estimates) || nrow(fit$estimates) == 0L) {
    pv_abort("Reportable stack_direct fit requires non-empty estimates.")
  }
  if (!"term" %in% names(fit$estimates)) {
    pv_abort("Reportable stack_direct fit estimates must include a `term` column.")
  }
  if (!identical(as.character(fit$estimates$term), fit$target$fe_names)) {
    pv_abort("Reportable stack_direct fit estimate terms must align with target fixed-effect names.")
  }
  estimate_required <- c(
    "estimate", "se", "std.error", "df", "df_method", "conf_level",
    "conf_low", "conf_high", "conf.low", "conf.high", "interval_role",
    "coverage_claim_allowed", "parameter_scope", "target_source",
    "target_hash"
  )
  missing_estimate <- setdiff(estimate_required, names(fit$estimates))
  if (length(missing_estimate) > 0L) {
    pv_abort(sprintf("Reportable stack_direct fit estimates are missing required column(s): %s.", paste(missing_estimate, collapse = ", ")))
  }
  if (!identical(as.character(fit$estimates$parameter_scope), rep("fixed_effect", length(fit$target$fe_names)))) {
    pv_abort("Reportable stack_direct fit estimates must use `parameter_scope = \"fixed_effect\"`.")
  }
  if (!identical(as.character(fit$estimates$target_source), rep(fit$target$target_source, length(fit$target$fe_names))) ||
      !identical(as.character(fit$estimates$target_hash), rep(fit$target$target_hash, length(fit$target$fe_names)))) {
    pv_abort("Reportable stack_direct fit estimate target metadata must match the fit target.")
  }
  if (!identical(fit$estimates$std.error, fit$estimates$se)) {
    pv_abort("Reportable stack_direct fit estimate `std.error` must match `se`.")
  }
  if (!identical(fit$estimates$conf.low, fit$estimates$conf_low) ||
      !identical(fit$estimates$conf.high, fit$estimates$conf_high)) {
    pv_abort("Reportable stack_direct fit estimate interval aliases must match.")
  }
  expected_estimate <- unname(fit$target$beta[fit$target$fe_names])
  expected_se <- unname(sqrt(diag(fit$target$T_MI))[fit$target$fe_names])
  expected_df <- unname(fit$target$df[fit$target$fe_names])
  if (max(abs(fit$estimates$estimate - expected_estimate)) > 1e-8) {
    pv_abort("Reportable stack_direct fit estimates must match external target fixed-effect centers.")
  }
  if (max(abs(fit$estimates$se - expected_se)) > 1e-8) {
    pv_abort("Reportable stack_direct fit standard errors must match the target covariance.")
  }
  if (max(abs(fit$estimates$df - expected_df)) > 1e-8) {
    pv_abort("Reportable stack_direct fit degrees of freedom must match the target.")
  }
  if (!is.null(fit$draws) && !identical(fit$draws, fit$ccc$draws_fe_cal)) {
    pv_abort("Reportable stack_direct fit retained draws must match CCC fixed-effect calibrated draws.")
  }
  invisible(fit)
}

pv_validate_per_pv_reportable_fit <- function(fit) {
  if (!identical(fit$method, "per_pv") ||
      !fit$status %in% c("ok", "warning")) {
    return(invisible(fit))
  }
  if (is.null(fit$target) || !inherits(fit$target, "pvstackr_reference_pool")) {
    pv_abort("Reportable per_pv fit requires a non-null `pvstackr_reference_pool` target.")
  }
  if (!is.null(fit$stack_fit) || !is.null(fit$ccc)) {
    pv_abort("Reportable per_pv fit must not carry stack-fit or CCC components.")
  }
  if (!is.null(fit$draws)) {
    pv_abort("Reportable per_pv fit must keep per-PV draws nested in diagnostics, not top-level `draws`.")
  }
  validate_pvstackr_reference_pool(fit$target)
  if (!is.data.frame(fit$estimates) || nrow(fit$estimates) == 0L) {
    pv_abort("Reportable per_pv fit requires non-empty estimates.")
  }
  estimate_required <- c(
    "term", "estimate", "se", "std.error", "df", "df_method", "df_complete",
    "conf_level", "conf_low", "conf_high", "conf.low", "conf.high",
    "interval_role", "coverage_claim_allowed", "parameter_scope",
    "target_source", "target_hash", "pooling_source", "pooling_hash"
  )
  missing_estimate <- setdiff(estimate_required, names(fit$estimates))
  if (length(missing_estimate) > 0L) {
    pv_abort(sprintf("Reportable per_pv fit estimates are missing required column(s): %s.", paste(missing_estimate, collapse = ", ")))
  }
  if (!identical(as.character(fit$estimates$term), fit$target$fe_names)) {
    pv_abort("Reportable per_pv fit estimate terms must align with pooled fixed-effect names.")
  }
  if (!identical(as.character(fit$estimates$parameter_scope), rep("fixed_effect", length(fit$target$fe_names)))) {
    pv_abort("Reportable per_pv fit estimates must use `parameter_scope = \"fixed_effect\"`.")
  }
  if (!identical(as.character(fit$estimates$target_source), rep(fit$target$target_source, length(fit$target$fe_names))) ||
      !identical(as.character(fit$estimates$target_hash), rep(fit$target$target_hash, length(fit$target$fe_names))) ||
      !identical(as.character(fit$estimates$pooling_source), rep(fit$target$target_source, length(fit$target$fe_names))) ||
      !identical(as.character(fit$estimates$pooling_hash), rep(fit$target$target_hash, length(fit$target$fe_names)))) {
    pv_abort("Reportable per_pv fit estimate pooling metadata must match the reference pool.")
  }
  if (!identical(fit$estimates$std.error, fit$estimates$se) ||
      !identical(fit$estimates$conf.low, fit$estimates$conf_low) ||
      !identical(fit$estimates$conf.high, fit$estimates$conf_high)) {
    pv_abort("Reportable per_pv fit estimate aliases must match.")
  }
  expected_estimate <- unname(fit$target$beta[fit$target$fe_names])
  expected_se <- unname(fit$target$se[fit$target$fe_names])
  expected_df <- unname(fit$target$df[fit$target$fe_names])
  if (max(abs(fit$estimates$estimate - expected_estimate)) > 1e-8) {
    pv_abort("Reportable per_pv fit estimates must match Rubin pooled centers.")
  }
  if (max(abs(fit$estimates$se - expected_se)) > 1e-8) {
    pv_abort("Reportable per_pv fit standard errors must match Rubin pooled standard errors.")
  }
  if (max(abs(fit$estimates$df - expected_df)) > 1e-8) {
    pv_abort("Reportable per_pv fit degrees of freedom must match Rubin pooled degrees of freedom.")
  }
  invisible(fit)
}

pv_validate_stack_psis_reportable_fit <- function(fit) {
  if (!identical(fit$method, "stack_psis")) {
    return(invisible(fit))
  }
  if (identical(fit$status, "blocked")) {
    if (!is.list(fit$diagnostics) || !"psis" %in% names(fit$diagnostics)) {
      pv_abort("Blocked stack_psis fit requires PSIS diagnostics.")
    }
    if (is.data.frame(fit$estimates) && nrow(fit$estimates) > 0L) {
      pv_abort("Blocked stack_psis fit must not include reportable estimates.")
    }
    return(invisible(fit))
  }
  if (!fit$status %in% c("ok", "warning")) {
    return(invisible(fit))
  }
  if (!is.null(fit$ccc)) {
    pv_abort("Reportable stack_psis fit must not carry a CCC component.")
  }
  if (!is.null(fit$draws)) {
    pv_abort("Reportable stack_psis fit must keep weighted draws nested in diagnostics, not top-level `draws`.")
  }
  if (!is.list(fit$diagnostics) ||
      !all(c("psis", "pooling", "weighted") %in% names(fit$diagnostics))) {
    pv_abort("Reportable stack_psis fit requires PSIS, pooling, and weighted diagnostics.")
  }
  psis <- fit$diagnostics$psis
  pooling <- fit$diagnostics$pooling
  pv_assert_named_list(psis, "diagnostics$psis")
  pv_assert_named_list(pooling, "diagnostics$pooling")
  if (!all(c("pareto_k", "threshold", "pareto_k_max", "status") %in% names(psis))) {
    pv_abort("Reportable stack_psis fit PSIS diagnostics are incomplete.")
  }
  if (!is.numeric(psis$pareto_k) || length(psis$pareto_k) < 2L || any(is.na(psis$pareto_k))) {
    pv_abort("Reportable stack_psis fit Pareto-k diagnostics must be numeric and complete.")
  }
  if (identical(psis$status, "failed") && identical(fit$status, "ok")) {
    pv_abort("Reportable stack_psis fit with failed PSIS diagnostics must not use `status = \"ok\"`.")
  }
  if (!is.data.frame(fit$estimates) || nrow(fit$estimates) == 0L) {
    pv_abort("Reportable stack_psis fit requires non-empty estimates.")
  }
  estimate_required <- c(
    "term", "estimate", "se", "std.error", "df", "df_method", "df_complete",
    "conf_level", "conf_low", "conf_high", "conf.low", "conf.high",
    "interval_role", "coverage_claim_allowed", "parameter_scope",
    "target_source", "target_hash", "pooling_source", "pooling_hash",
    "psis_status", "pareto_k_max", "psis_k_threshold"
  )
  missing_estimate <- setdiff(estimate_required, names(fit$estimates))
  if (length(missing_estimate) > 0L) {
    pv_abort(sprintf("Reportable stack_psis fit estimates are missing required column(s): %s.", paste(missing_estimate, collapse = ", ")))
  }
  fe_names <- names(pooling$beta)
  if (is.null(fe_names) || length(fe_names) == 0L ||
      !identical(as.character(fit$estimates$term), fe_names)) {
    pv_abort("Reportable stack_psis fit estimate terms must align with pooled fixed-effect names.")
  }
  if (!identical(as.character(fit$estimates$parameter_scope), rep("fixed_effect", length(fe_names)))) {
    pv_abort("Reportable stack_psis fit estimates must use `parameter_scope = \"fixed_effect\"`.")
  }
  if (!identical(fit$estimates$std.error, fit$estimates$se) ||
      !identical(fit$estimates$conf.low, fit$estimates$conf_low) ||
      !identical(fit$estimates$conf.high, fit$estimates$conf_high)) {
    pv_abort("Reportable stack_psis fit estimate aliases must match.")
  }
  if (!identical(as.character(fit$estimates$pooling_hash), rep(pooling$pooling_hash, length(fe_names)))) {
    pv_abort("Reportable stack_psis fit estimate pooling metadata must match pooling diagnostics.")
  }
  if (max(abs(fit$estimates$estimate - unname(pooling$beta[fe_names]))) > 1e-8) {
    pv_abort("Reportable stack_psis fit estimates must match PSIS Rubin pooled centers.")
  }
  if (max(abs(fit$estimates$se - unname(sqrt(diag(pooling$T_MI))[fe_names]))) > 1e-8) {
    pv_abort("Reportable stack_psis fit standard errors must match PSIS Rubin pooled covariance.")
  }
  if (max(abs(fit$estimates$df - unname(pooling$df[fe_names]))) > 1e-8) {
    pv_abort("Reportable stack_psis fit degrees of freedom must match PSIS Rubin pooled degrees of freedom.")
  }
  invisible(fit)
}

new_pvstackr_fit <- function(
  method,
  design = NULL,
  target = NULL,
  stack_fit = NULL,
  ccc = NULL,
  estimates = data.frame(),
  draws = NULL,
  diagnostics = list(),
  status = "ok",
  control = NULL,
  reason_codes = character(),
  provenance = list(),
  warnings = character()
) {
  method <- pv_validate_method(method)
  control <- if (is.null(control)) pv_control(method = method) else pv_validate_control(control)
  if (!identical(control$method, method)) {
    pv_abort("Fit `control$method` must match `method`.")
  }
  status <- pv_assert_scalar_string(status, "status")
  if (!status %in% c("ok", "warning", "blocked")) {
    pv_abort("Fit `status` must be one of: ok, warning, blocked.")
  }
  pv_validate_method_implemented_status(method, status)
  reason_codes <- pv_validate_fit_reason_codes(status, reason_codes)
  diagnostics <- pv_validate_named_list_field(diagnostics, "diagnostics")
  provenance <- pv_validate_named_list_field(provenance, "provenance")
  warnings <- pv_validate_fit_warnings(status, warnings)
  if (!is.data.frame(estimates)) {
    pv_abort("Fit `estimates` must be a data frame.")
  }
  if (!is.null(draws)) {
    draws <- ccc_as_draw_matrix(draws)
    non_fe <- colnames(draws)[!grepl("^b_", colnames(draws))]
    if (length(non_fe) > 0L) {
      pv_abort("Fit `draws`, when retained, must contain reportable fixed-effect columns only.")
    }
  }

  fit <- list(
    method = method,
    design = design,
    target = target,
    stack_fit = stack_fit,
    ccc = ccc,
    estimates = estimates,
    draws = draws,
    diagnostics = diagnostics,
    status = status,
    control = control,
    reason_codes = reason_codes,
    schema_version = pv_schema_version(),
    provenance = c(
      pv_provenance("new_pvstackr_fit"),
      provenance
    ),
    warnings = warnings
  )
  class(fit) <- c("pvstackr_fit", "list")
  validate_pvstackr_fit(fit)
  fit
}

validate_pvstackr_fit <- function(fit) {
  pv_assert_named_list(fit, "fit")
  required <- c(
    "method", "design", "target", "stack_fit", "ccc", "estimates",
    "draws", "diagnostics", "status", "control", "reason_codes",
    "schema_version", "provenance", "warnings"
  )
  missing <- setdiff(required, names(fit))
  if (length(missing) > 0L) {
    pv_abort(sprintf("Fit object is missing required field(s): %s.", paste(missing, collapse = ", ")))
  }
  method <- pv_validate_method(fit$method)
  pv_validate_control(fit$control)
  if (!identical(fit$control$method, method)) {
    pv_abort("Fit `control$method` must match `method`.")
  }
  status <- pv_assert_scalar_string(fit$status, "status")
  if (!status %in% c("ok", "warning", "blocked")) {
    pv_abort("Fit `status` must be one of: ok, warning, blocked.")
  }
  pv_validate_method_implemented_status(method, status)
  reason_codes <- fit$reason_codes %||% character()
  reason_codes <- pv_validate_fit_reason_codes(status, reason_codes)
  if (!is.null(fit$design)) {
    validate_pvstackr_design(fit$design)
  }
  if (!is.null(fit$target)) {
    validate_pvstackr_fit_target(fit$target)
  }
  if (!is.null(fit$stack_fit)) {
    validate_pvstackr_stack_fit(fit$stack_fit)
  }
  if (!is.null(fit$ccc)) {
    validate_pvstackr_ccc(fit$ccc)
  }
  if (!is.null(fit$target) && !is.null(fit$ccc) &&
      !identical(fit$target$target_hash, fit$ccc$target_hash)) {
    pv_abort("Fit target and CCC object must carry the same `target_hash`.")
  }
  if (!is.data.frame(fit$estimates)) {
    pv_abort("Fit `estimates` must be a data frame.")
  }
  if (nrow(fit$estimates) > 0L && "term" %in% names(fit$estimates)) {
    non_fe_terms <- fit$estimates$term[!grepl("^b_", fit$estimates$term)]
    if (length(non_fe_terms) > 0L) {
      pv_abort("Fit `estimates` must contain reportable fixed-effect terms only.")
    }
  }
  if (!is.null(fit$draws)) {
    draws <- ccc_as_draw_matrix(fit$draws)
    non_fe <- colnames(draws)[!grepl("^b_", colnames(draws))]
    if (length(non_fe) > 0L) {
      pv_abort("Fit `draws`, when retained, must contain reportable fixed-effect columns only.")
    }
  }
  pv_validate_stack_direct_reportable_fit(fit)
  pv_validate_per_pv_reportable_fit(fit)
  pv_validate_stack_psis_reportable_fit(fit)
  pv_validate_named_list_field(fit$diagnostics, "diagnostics")
  pv_validate_schema_version(fit$schema_version)
  pv_validate_named_list_field(fit$provenance, "provenance")
  pv_validate_fit_warnings(status, fit$warnings)
  invisible(fit)
}
