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

pv_design_fields <- function() {
  c(
    "data", "data_manifest", "formula", "formula_string", "pv_cols",
    "weight_col", "rep_weight_cols", "fay_k", "id_cols", "M", "R", "n",
    "row_support_hash", "row_support", "pv_value_hash", "weight_design_hash",
    "design_hash", "roles", "created_at", "schema_version", "provenance",
    "warnings"
  )
}

pv_design_canonicalize_formula <- function(design) {
  validate_pvstackr_design(design)
  out <- design
  out$formula <- tryCatch(
    stats::as.formula(out$formula_string, env = baseenv()),
    error = function(error) NULL
  )
  if (is.null(out$formula)) {
    pv_abort("Fit design formula cannot be canonicalized for retention.")
  }
  validate_pvstackr_design(out)
  out
}

pv_design_data_free_snapshot <- function(design) {
  validate_pvstackr_design(design)
  formula <- tryCatch(
    stats::as.formula(design$formula_string, env = baseenv()),
    error = function(error) NULL
  )
  if (is.null(formula)) {
    pv_abort("A data-free design snapshot requires a canonical formula string.")
  }
  retained_fields <- c(
    "source", "target_hash", "target_manifest_hash", "target_content_hash",
    "binding_verification_policy", "pooling_hash"
  )
  retained <- design$provenance[
    retained_fields[retained_fields %in% names(design$provenance)]
  ]
  retained <- retained[vapply(
    retained,
    function(value) {
      is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
    },
    logical(1)
  )]
  roles <- list(
    outcome_placeholder = unname(as.character(
      design$roles$outcome_placeholder %||% "OUTCOME"
    )),
    method = unname(as.character(design$roles$method %||% "unknown"))
  )
  provenance <- c(
    list(
      function_name = "pv_design_data_free_snapshot",
      package = "pvstackr",
      schema_version = "0.2.0",
      retention_policy = "raw_data_redacted"
    ),
    retained
  )
  hashes <- list(
    row_support_hash = unname(as.character(design$row_support_hash)),
    pv_value_hash = unname(as.character(design$pv_value_hash)),
    weight_design_hash = unname(as.character(design$weight_design_hash))
  )
  formula_info <- pv_binding_formula_rhs_info(
    design$formula,
    data = design$data
  )
  bound_columns <- unique(c(
    design$id_cols,
    design$pv_cols,
    formula_info$data_variables,
    design$weight_col,
    design$rep_weight_cols
  ))
  bound_columns <- bound_columns[bound_columns %in% names(design$data)]
  data_manifest <- design$data_manifest
  data_manifest$columns <- bound_columns
  data_manifest$column_classes <- design$data_manifest$column_classes[
    bound_columns
  ]
  data_manifest$formula_string <- design$formula_string
  design_hash <- pv_hash_payload(c(data_manifest, hashes))
  out <- list(
    data = NULL,
    data_manifest = data_manifest,
    formula = formula,
    formula_string = unname(as.character(design$formula_string)),
    pv_cols = unname(as.character(design$pv_cols)),
    weight_col = if (is.null(design$weight_col)) NULL else
      unname(as.character(design$weight_col)),
    rep_weight_cols = unname(as.character(design$rep_weight_cols)),
    fay_k = as.numeric(design$fay_k),
    id_cols = unname(as.character(design$id_cols)),
    M = as.integer(design$M),
    R = as.integer(design$R),
    n = as.integer(design$n),
    row_support_hash = hashes$row_support_hash,
    row_support = list(
      type = unname(as.character(design$row_support$type)),
      id_cols = unname(as.character(design$row_support$id_cols)),
      n = as.integer(design$row_support$n),
      hash = hashes$row_support_hash
    ),
    pv_value_hash = hashes$pv_value_hash,
    weight_design_hash = hashes$weight_design_hash,
    design_hash = design_hash,
    roles = roles,
    created_at = "omitted_for_data_free_snapshot",
    schema_version = "0.2.0",
    provenance = provenance,
    warnings = unname(as.character(design$warnings))
  )
  if (!identical(names(out), pv_design_fields())) {
    pv_abort("Internal data-free design snapshot fields are not canonical.")
  }
  class(out) <- c("pvstackr_design", "list")
  validate_pvstackr_design(out)
  out
}

pv_design_target_data_free_snapshot <- function(target, binding_proof) {
  validate_pvstackr_brr_target(target)
  pv_binding_proof_validate(
    binding_proof,
    target_manifest = target$binding_manifest
  )
  components <- target$binding_manifest$components
  component_hashes <- target$binding_manifest$component_hashes
  id_cols <- if (is.null(target$id_cols)) character() else
    unname(as.character(target$id_cols))
  out <- list(
    data = NULL,
    data_manifest = list(
      snapshot_type = "target_bound",
      contract_id = target$binding_manifest$contract_id,
      manifest_hash = target$binding_manifest$manifest_hash,
      target_content_hash = target$target_content$target_content_hash,
      row_hash = component_hashes$row_hash,
      pv_hash = component_hashes$pv_hash,
      formula_hash = component_hashes$formula_hash,
      model_matrix_hash = component_hashes$model_matrix_hash,
      weight_hash = component_hashes$weight_hash,
      n = components$row$n,
      formula_string = target$formula_string,
      pv_cols = target$pv_cols,
      weight_col = target$weight_col,
      rep_weight_cols = target$rep_weight_cols,
      id_cols = id_cols,
      fay_k = target$fay_k
    ),
    formula = target$formula,
    formula_string = target$formula_string,
    pv_cols = target$pv_cols,
    weight_col = target$weight_col,
    rep_weight_cols = target$rep_weight_cols,
    fay_k = target$fay_k,
    id_cols = id_cols,
    M = target$M,
    R = target$R,
    n = components$row$n,
    row_support_hash = component_hashes$row_hash,
    row_support = list(
      type = if (length(id_cols) > 0L) "id_cols" else "row_number",
      id_cols = id_cols,
      n = components$row$n,
      hash = component_hashes$row_hash
    ),
    pv_value_hash = component_hashes$pv_hash,
    weight_design_hash = component_hashes$weight_hash,
    design_hash = target$binding_manifest$manifest_hash,
    roles = list(
      outcome_placeholder = "OUTCOME",
      method = "stack_direct"
    ),
    created_at = target$provenance$assembled_at,
    schema_version = "0.2.0",
    provenance = list(
      function_name = "pv_design_data_free_snapshot",
      package = "pvstackr",
      schema_version = "0.2.0",
      retention_policy = "raw_data_redacted",
      source = "pv_fit_direct",
      target_hash = target$target_hash,
      target_manifest_hash = binding_proof$target_manifest_hash,
      target_content_hash = target$target_content$target_content_hash,
      binding_verification_policy = binding_proof$verification_policy
    ),
    warnings = character()
  )
  if (!identical(names(out), pv_design_fields())) {
    pv_abort("Internal target-bound design snapshot fields are not canonical.")
  }
  class(out) <- c("pvstackr_design", "list")
  validate_pvstackr_design(out)
  out
}

pv_validate_design_target_data_free_snapshot <- function(design) {
  manifest_fields <- c(
    "snapshot_type", "contract_id", "manifest_hash", "target_content_hash",
    "row_hash", "pv_hash", "formula_hash", "model_matrix_hash",
    "weight_hash", "n", "formula_string", "pv_cols", "weight_col",
    "rep_weight_cols", "id_cols", "fay_k"
  )
  manifest <- design$data_manifest
  if (!is.list(manifest) || !identical(names(manifest), manifest_fields) ||
      !identical(attributes(manifest), list(names = manifest_fields)) ||
      !identical(manifest$snapshot_type, "target_bound") ||
      !identical(manifest$contract_id, "pvstackr_data_binding_v1") ||
      !identical(manifest$n, design$n) ||
      !identical(manifest$formula_string, design$formula_string) ||
      !identical(manifest$pv_cols, design$pv_cols) ||
      !identical(manifest$weight_col, design$weight_col) ||
      !identical(manifest$rep_weight_cols, design$rep_weight_cols) ||
      !identical(manifest$id_cols, design$id_cols) ||
      !identical(manifest$fay_k, design$fay_k)) {
    pv_abort("Target-bound data-free design manifest is not canonical.")
  }
  sha_fields <- c(
    "manifest_hash", "target_content_hash", "row_hash", "pv_hash",
    "formula_hash", "model_matrix_hash", "weight_hash"
  )
  sha_values <- unlist(manifest[sha_fields], use.names = FALSE)
  if (!is.character(sha_values) || anyNA(sha_values) ||
      any(!grepl("^sha256:[0-9a-f]{64}$", sha_values)) ||
      !identical(design$row_support_hash, manifest$row_hash) ||
      !identical(design$pv_value_hash, manifest$pv_hash) ||
      !identical(design$weight_design_hash, manifest$weight_hash) ||
      !identical(design$design_hash, manifest$manifest_hash)) {
    pv_abort("Target-bound data-free design hashes are malformed or inconsistent.")
  }
  if (!is.character(design$pv_cols) || length(design$pv_cols) < 1L ||
      anyNA(design$pv_cols) || anyDuplicated(design$pv_cols) ||
      !identical(design$M, as.integer(length(design$pv_cols))) ||
      !is.character(design$rep_weight_cols) ||
      length(design$rep_weight_cols) < 2L || anyNA(design$rep_weight_cols) ||
      anyDuplicated(design$rep_weight_cols) ||
      !identical(design$R, as.integer(length(design$rep_weight_cols))) ||
      !is.integer(design$n) || length(design$n) != 1L ||
      is.na(design$n) || design$n < 1L ||
      !is.character(design$id_cols) || anyNA(design$id_cols) ||
      anyDuplicated(design$id_cols)) {
    pv_abort("Target-bound data-free design role or dimension metadata is malformed.")
  }
  pv_validate_fay_k(design$fay_k)
  expected_row_support <- list(
    type = if (length(design$id_cols) > 0L) "id_cols" else "row_number",
    id_cols = design$id_cols,
    n = design$n,
    hash = design$row_support_hash
  )
  if (!identical(design$row_support, expected_row_support) ||
      !identical(
        design$roles,
        list(outcome_placeholder = "OUTCOME", method = "stack_direct")
      )) {
    pv_abort("Target-bound data-free design row support or roles are not canonical.")
  }
  provenance_fields <- c(
    "function_name", "package", "schema_version", "retention_policy",
    "source", "target_hash", "target_manifest_hash", "target_content_hash",
    "binding_verification_policy"
  )
  if (!is.list(design$provenance) ||
      !identical(names(design$provenance), provenance_fields) ||
      !identical(
        attributes(design$provenance),
        list(names = provenance_fields)
      ) ||
      !identical(design$provenance$function_name, "pv_design_data_free_snapshot") ||
      !identical(design$provenance$package, "pvstackr") ||
      !identical(design$provenance$schema_version, "0.2.0") ||
      !identical(design$provenance$retention_policy, "raw_data_redacted") ||
      !identical(design$provenance$source, "pv_fit_direct") ||
      !identical(
        design$provenance$target_manifest_hash,
        manifest$manifest_hash
      ) || !identical(
        design$provenance$target_content_hash,
        manifest$target_content_hash
      ) || any(!vapply(
        design$provenance,
        function(value) {
          is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
        },
        logical(1)
      ))) {
    pv_abort("Target-bound data-free design provenance is not canonical.")
  }
  pv_assert_scalar_string(design$created_at, "created_at")
  if (!identical(design$warnings, character())) {
    pv_abort("Target-bound data-free design warnings must be empty.")
  }
  invisible(design)
}

pv_validate_design_data_free_snapshot <- function(design) {
  fields <- pv_design_fields()
  root_attributes <- attributes(design)
  if (!is.list(design) || !identical(names(design), fields) ||
      !identical(names(root_attributes), c("names", "class")) ||
      !identical(root_attributes$names, fields) ||
      !identical(root_attributes$class, c("pvstackr_design", "list")) ||
      !is.null(design$data) || !identical(design$schema_version, "0.2.0")) {
    pv_abort("Data-free design snapshot must use the exact redacted schema-0.2 fields.")
  }
  canonical_formula <- tryCatch(
    stats::as.formula(design$formula_string, env = baseenv()),
    error = function(error) NULL
  )
  if (is.null(canonical_formula) || !identical(design$formula, canonical_formula) ||
      !identical(design$formula_string, pv_formula_string(design$formula))) {
    pv_abort("Data-free design snapshot formula must be the exact base-environment projection.")
  }
  if (identical(design$data_manifest$snapshot_type, "target_bound")) {
    return(pv_validate_design_target_data_free_snapshot(design))
  }
  manifest_fields <- c(
    "n", "columns", "column_classes", "formula_string", "pv_cols",
    "weight_col", "rep_weight_cols", "id_cols", "fay_k"
  )
  if (!is.list(design$data_manifest) ||
      !identical(names(design$data_manifest), manifest_fields) ||
      !identical(
        attributes(design$data_manifest),
        list(names = manifest_fields)
      ) ||
      !identical(design$data_manifest$n, design$n) ||
      !identical(design$data_manifest$formula_string, design$formula_string) ||
      !identical(design$data_manifest$pv_cols, design$pv_cols) ||
      !identical(design$data_manifest$weight_col, design$weight_col) ||
      !identical(design$data_manifest$rep_weight_cols, design$rep_weight_cols) ||
      !identical(design$data_manifest$id_cols, design$id_cols) ||
      !identical(design$data_manifest$fay_k, design$fay_k) ||
      !is.character(design$data_manifest$columns) ||
      !is.character(design$data_manifest$column_classes) ||
      !identical(
        names(design$data_manifest$column_classes),
        design$data_manifest$columns
      )) {
    pv_abort("Data-free design snapshot manifest is not canonical.")
  }
  if (!is.character(design$pv_cols) || length(design$pv_cols) < 1L ||
      anyNA(design$pv_cols) || any(!nzchar(design$pv_cols)) ||
      anyDuplicated(design$pv_cols) ||
      !identical(design$M, as.integer(length(design$pv_cols))) ||
      !is.character(design$rep_weight_cols) || anyNA(design$rep_weight_cols) ||
      anyDuplicated(design$rep_weight_cols) ||
      !identical(design$R, as.integer(length(design$rep_weight_cols))) ||
      !is.integer(design$n) || length(design$n) != 1L || is.na(design$n) ||
      design$n < 1L ||
      !(is.null(design$weight_col) ||
        (is.character(design$weight_col) && length(design$weight_col) == 1L &&
          !is.na(design$weight_col) && nzchar(design$weight_col))) ||
      !is.character(design$id_cols) || anyNA(design$id_cols) ||
      anyDuplicated(design$id_cols)) {
    pv_abort("Data-free design snapshot role and dimension metadata is malformed.")
  }
  pv_validate_fay_k(design$fay_k)
  hash_fields <- c(
    "row_support_hash", "pv_value_hash", "weight_design_hash", "design_hash"
  )
  hashes <- unlist(design[hash_fields], use.names = FALSE)
  if (!is.character(hashes) || anyNA(hashes) ||
      any(!grepl("^[0-9a-f]{8}$", hashes))) {
    pv_abort("Data-free design snapshot hashes must be canonical legacy digests.")
  }
  expected_design_hash <- pv_hash_payload(c(
    design$data_manifest,
    list(
      row_support_hash = design$row_support_hash,
      pv_value_hash = design$pv_value_hash,
      weight_design_hash = design$weight_design_hash
    )
  ))
  expected_row_support <- list(
    type = if (length(design$id_cols) > 0L) "id_cols" else "row_number",
    id_cols = design$id_cols,
    n = design$n,
    hash = design$row_support_hash
  )
  if (!identical(design$design_hash, expected_design_hash) ||
      !identical(design$row_support, expected_row_support)) {
    pv_abort("Data-free design snapshot hashes or row support are stale.")
  }
  if (!identical(
    design$roles,
    list(
      outcome_placeholder = "OUTCOME",
      method = design$roles$method
    )
  ) || !is.character(design$roles$method) ||
      length(design$roles$method) != 1L || is.na(design$roles$method) ||
      !design$roles$method %in% c("stack_direct", "per_pv", "stack_psis")) {
    pv_abort("Data-free design snapshot roles are not canonical.")
  }
  provenance_fields <- c(
    "function_name", "package", "schema_version", "retention_policy",
    "source", "target_hash", "target_manifest_hash", "target_content_hash",
    "binding_verification_policy", "pooling_hash"
  )
  if (!is.list(design$provenance) ||
      any(!names(design$provenance) %in% provenance_fields) ||
      !identical(
        attributes(design$provenance),
        list(names = names(design$provenance))
      ) ||
      !identical(
        names(design$provenance)[seq_len(4L)],
        provenance_fields[seq_len(4L)]
      ) ||
      !identical(design$provenance$function_name, "pv_design_data_free_snapshot") ||
      !identical(design$provenance$package, "pvstackr") ||
      !identical(design$provenance$schema_version, "0.2.0") ||
      !identical(design$provenance$retention_policy, "raw_data_redacted") ||
      any(!vapply(
        design$provenance,
        function(value) {
          is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
        },
        logical(1)
      ))) {
    pv_abort("Data-free design snapshot provenance is not canonical or scalar-only.")
  }
  pv_assert_scalar_string(design$created_at, "created_at")
  pv_validate_character_field(design$warnings, "warnings")
  invisible(design)
}

pv_validate_fit_data_retention_control <- function(control) {
  pv_validate_control(control)
  if (isTRUE(control$keep_backend_fit) && !isTRUE(control$keep_data)) {
    pv_abort(
      "`keep_backend_fit = TRUE` requires `keep_data = TRUE` because an opaque backend object may retain analysis data."
    )
  }
  invisible(control)
}

pv_stack_fit_draw_projection <- function(stack_fit, return_draws) {
  if (is.null(stack_fit)) {
    return(NULL)
  }
  validate_pvstackr_stack_fit(stack_fit)
  return_draws <- pv_assert_scalar_logical(return_draws, "return_draws")
  out <- stack_fit
  if (isTRUE(return_draws)) {
    if (is.null(out$stacked_draws)) {
      pv_abort("Requested retained stack draws are not available.")
    }
    out$control$return_draws <- TRUE
  } else {
    out["stacked_draws"] <- list(NULL)
    out$control$return_draws <- FALSE
  }
  validate_pvstackr_stack_fit(out)
  out
}

pv_ccc_draw_projection <- function(ccc, return_draws) {
  if (is.null(ccc)) {
    return(NULL)
  }
  validate_pvstackr_ccc(ccc)
  return_draws <- pv_assert_scalar_logical(return_draws, "return_draws")
  out <- ccc
  if (isTRUE(return_draws)) {
    if (is.null(out$draws_calibrated) || is.null(out$draws_fe_cal)) {
      pv_abort("Requested retained CCC draws are not available.")
    }
    out$provenance$draws_retained <- TRUE
  } else {
    out["draws_calibrated"] <- list(NULL)
    out["draws_fe_cal"] <- list(NULL)
    out$provenance$draws_retained <- FALSE
  }
  validate_pvstackr_ccc(out)
  out
}

pv_stack_fit_composite_projection <- function(
  stack_fit,
  control,
  canonicalize_formula = FALSE
) {
  if (is.null(stack_fit)) {
    return(NULL)
  }
  validate_pvstackr_stack_fit(stack_fit)
  pv_validate_fit_data_retention_control(control)
  canonicalize_formula <- pv_assert_scalar_logical(
    canonicalize_formula,
    "canonicalize_formula"
  )
  out <- stack_fit
  out["prepared_data"] <- list(NULL)
  out$control$keep_data <- FALSE
  if (isTRUE(control$keep_backend_fit)) {
    if (is.null(out$fit)) {
      pv_abort(
        "`keep_backend_fit = TRUE` requested an opaque backend object that is not available."
      )
    }
    out$control$keep_backend_fit <- TRUE
  } else {
    out["fit"] <- list(NULL)
    out$control$keep_backend_fit <- FALSE
  }
  if (isTRUE(control$keep_log_lik)) {
    out$control$keep_log_lik <- TRUE
  } else {
    out["log_lik"] <- list(NULL)
    out$control$keep_log_lik <- FALSE
    out$meta$log_lik_retained <- FALSE
    out$provenance$log_lik_retained <- FALSE
  }
  if (isTRUE(canonicalize_formula)) {
    out$formula <- tryCatch(
      stats::as.formula(out$formula_string, env = baseenv()),
      error = function(error) NULL
    )
    if (is.null(out$formula)) {
      pv_abort("Composite stack-fit formula cannot be canonicalized for retention.")
    }
  }
  pv_stack_fit_draw_projection(out, control$return_draws)
}

pv_validate_fit_data_retention <- function(fit) {
  control <- fit$control
  method <- fit$method
  pv_validate_fit_data_retention_control(control)
  diagnostic_fields <- if (identical(method, "stack_direct")) {
    c(
      "preflight", "sampler", "sampler_gate", "stack_fit",
      "stack_fit_warnings", "ccc", "redaction"
    )
  } else if (identical(method, "per_pv")) {
    c("reference", "pooling")
  } else if (identical(fit$status, "blocked")) {
    c("psis", "redaction")
  } else {
    c("psis", "pooling", "weighted")
  }
  direct_slim_fields <- c("preflight", "sampler", "sampler_gate", "redaction")
  diagnostic_names_valid <- identical(
    names(fit$diagnostics),
    diagnostic_fields[diagnostic_fields %in% names(fit$diagnostics)]
  ) || (identical(method, "stack_direct") && identical(
    names(fit$diagnostics),
    direct_slim_fields[direct_slim_fields %in% names(fit$diagnostics)]
  )) || length(fit$diagnostics) == 0L
  if (!is.list(fit$diagnostics) || !diagnostic_names_valid ||
      any(!names(fit$diagnostics) %in% union(diagnostic_fields, direct_slim_fields)) ||
      !identical(
        attributes(fit$diagnostics),
        if (length(fit$diagnostics) == 0L) NULL else
          list(names = names(fit$diagnostics))
      )) {
    pv_abort("Fit diagnostics must use the exact method-specific retention envelope.")
  }

  if (!is.null(fit$design)) {
    if (isTRUE(control$keep_data)) {
      if (!is.data.frame(fit$design$data) ||
          identical(fit$design$schema_version, "0.2.0")) {
        pv_abort("`keep_data = TRUE` requires a retained-data design when a design is present.")
      }
    } else if (!is.null(fit$design$data) ||
        !identical(fit$design$schema_version, "0.2.0")) {
      pv_abort("`keep_data = FALSE` requires a data-free design snapshot.")
    }
    if (!identical(fit$design$roles$method, method)) {
      pv_abort("Fit design retention role must match the fitted method.")
    }
    expected_source <- c(
      stack_direct = "pv_fit_direct",
      per_pv = "pv_fit_reference",
      stack_psis = "pv_fit_stack_psis"
    )[[method]]
    if (!isTRUE(control$keep_data) &&
        !identical(fit$design$provenance$source, expected_source)) {
      pv_abort("Data-free design snapshot source must match the fitted method.")
    }
  }

  if (identical(method, "per_pv") && length(fit$diagnostics) > 0L) {
    reference_fields <- c(
      "source", "topology", "M", "pv_cols", "fe_names", "draw_counts",
      "dropped_draw_columns", "map_sources", "backend_fits", "per_pv_draws"
    )
    pooling_fields <- c(
      "beta", "U_bar", "B", "T_MI", "lambda", "df", "df_classic",
      "df_method", "df_complete", "target_hash", "target_source"
    )
    if (!identical(names(fit$diagnostics$reference), reference_fields) ||
        !identical(names(fit$diagnostics$pooling), pooling_fields) ||
        !identical(
          attributes(fit$diagnostics$reference),
          list(names = reference_fields)
        ) || !identical(
          attributes(fit$diagnostics$pooling),
          list(names = pooling_fields)
        )) {
      pv_abort("Per-PV diagnostics use a noncanonical retention schema.")
    }
    backend_fits <- fit$diagnostics$reference$backend_fits
    per_pv_draws <- fit$diagnostics$reference$per_pv_draws
    source <- fit$diagnostics$reference$source
    if (!is.character(source) || length(source) != 1L || is.na(source) ||
        !source %in% c("draws", "injected_fit")) {
      pv_abort("Per-PV diagnostic source is not canonical.")
    }
    if (!isTRUE(control$keep_backend_fit) && !is.null(backend_fits)) {
      pv_abort("Per-PV fit retained backend objects without authorization.")
    }
    if (isTRUE(control$keep_backend_fit) &&
        identical(source, "injected_fit") && is.null(backend_fits)) {
      pv_abort("Per-PV fit omitted requested injected backend objects.")
    }
    if (identical(source, "injected_fit") && is.null(fit$design)) {
      pv_abort("Injected per-PV fits must retain a design or data-free snapshot.")
    }
    if (!identical(!is.null(per_pv_draws), control$return_draws)) {
      pv_abort("Per-PV draw retention must match `control$return_draws`.")
    }
    if (!is.null(per_pv_draws)) {
      pv_cols <- fit$diagnostics$reference$pv_cols
      fe_names <- fit$diagnostics$reference$fe_names
      if (!is.list(per_pv_draws) ||
          !identical(names(per_pv_draws), pv_cols) ||
          !identical(attributes(per_pv_draws), list(names = pv_cols))) {
        pv_abort("Retained per-PV draws must use the exact PV list envelope.")
      }
      for (pv in pv_cols) {
        draws <- per_pv_draws[[pv]]
        draw_attributes <- attributes(draws)
        if (!is.matrix(draws) || !is.numeric(draws) || any(!is.finite(draws)) ||
            !identical(colnames(draws), fe_names) ||
            !identical(names(draw_attributes), c("dim", "dimnames")) ||
            !is.null(rownames(draws))) {
          pv_abort("Retained per-PV draws must be canonical fixed-effect-only matrices.")
        }
        expected_beta <- unname(fit$target$per_pv$beta[pv, fe_names])
        expected_U <- fit$target$per_pv$U[[pv]][fe_names, fe_names, drop = FALSE]
        if (!identical(
              nrow(draws),
              as.integer(fit$target$per_pv$draw_counts[[pv]])
            ) ||
            max(abs(colMeans(draws) - expected_beta)) > 1e-10 ||
            max(abs(stats::cov(draws) - expected_U)) > 1e-10) {
          pv_abort("Retained per-PV draws must reproduce their pooled mean/covariance summaries.")
        }
      }
    }
    if (!is.null(fit$design) && !isTRUE(control$keep_data) &&
        !identical(
          fit$design$provenance$target_hash,
          fit$target$target_hash
        )) {
      pv_abort("Per-PV data-free design target hash is not linked to its pool.")
    }
    expected_reference <- list(
      source = source,
      topology = "one_fit_per_plausible_value",
      M = fit$target$M,
      pv_cols = fit$target$pv_cols,
      fe_names = fit$target$fe_names,
      draw_counts = fit$target$per_pv$draw_counts,
      dropped_draw_columns = fit$target$per_pv$dropped_draw_columns,
      map_sources = fit$target$per_pv$map_sources,
      backend_fits = backend_fits,
      per_pv_draws = per_pv_draws
    )
    expected_pooling <- list(
      beta = fit$target$beta,
      U_bar = fit$target$U_bar,
      B = fit$target$B,
      T_MI = fit$target$T_MI,
      lambda = fit$target$lambda,
      df = fit$target$df,
      df_classic = fit$target$df_classic,
      df_method = fit$target$df_method,
      df_complete = fit$target$df_complete,
      target_hash = fit$target$target_hash,
      target_source = fit$target$target_source
    )
    if (!identical(fit$diagnostics$reference, expected_reference) ||
        !identical(fit$diagnostics$pooling, expected_pooling)) {
      pv_abort("Per-PV diagnostics must exactly mirror the validated reference pool and authorized retained payloads.")
    }
  }

  if (identical(method, "stack_psis") && !is.null(fit$stack_fit)) {
    if (!is.null(fit$stack_fit$prepared_data) ||
        !is.null(fit$stack_fit$stacked_draws) ||
        !identical(fit$stack_fit$control$keep_data, FALSE) ||
        !identical(fit$stack_fit$control$return_draws, FALSE) ||
        !identical(
          fit$stack_fit$control$keep_backend_fit,
          control$keep_backend_fit
        ) || !identical(
          fit$stack_fit$control$keep_log_lik,
          control$keep_log_lik
        )) {
      pv_abort("Composite stack_psis fit has inconsistent nested data retention.")
    }
    canonical_formula <- tryCatch(
      stats::as.formula(fit$stack_fit$formula_string, env = baseenv()),
      error = function(error) NULL
    )
    if (is.null(canonical_formula) ||
        !identical(fit$stack_fit$formula, canonical_formula)) {
      pv_abort("Composite stack_psis formula must not retain a caller environment.")
    }
  }
  if (identical(method, "stack_psis") &&
      identical(fit$provenance$stacked_source, "injected_fit") &&
      !identical(fit$status, "blocked") && is.null(fit$design)) {
    pv_abort("Injected stack_psis fits must retain a design or data-free snapshot.")
  }
  if (identical(method, "stack_psis") && !is.null(fit$design) &&
      !isTRUE(control$keep_data)) {
    pooling_hash <- fit$provenance$pooling_hash
    if (!identical(fit$design$provenance$pooling_hash, pooling_hash) ||
        !identical(fit$diagnostics$pooling$pooling_hash, pooling_hash)) {
      pv_abort("stack_psis data-free design pooling hash is not canonically linked.")
    }
  }
  if (identical(method, "stack_psis") &&
      !identical(fit$status, "blocked") &&
      is.list(fit$diagnostics$weighted)) {
    weighted_draws <- fit$diagnostics$weighted$proposal_draws
    weights <- fit$diagnostics$weighted$weights
    if (!identical(!is.null(weighted_draws), control$return_draws) ||
        !identical(!is.null(weights), control$return_draws)) {
      pv_abort("stack_psis retained draw/weight payload must match `control$return_draws`.")
    }
    if (!is.null(weighted_draws)) {
      draw_attributes <- attributes(weighted_draws)
      fe_names <- colnames(fit$diagnostics$weighted$beta)
      if (!is.matrix(weighted_draws) || !is.numeric(weighted_draws) ||
          any(!is.finite(weighted_draws)) ||
          !identical(colnames(weighted_draws), fe_names) ||
          !identical(names(draw_attributes), c("dim", "dimnames")) ||
          !is.null(rownames(weighted_draws))) {
        pv_abort("stack_psis retained source draws must be canonical fixed-effect-only matrices.")
      }
    }
    if (!is.null(weights)) {
      weight_attributes <- attributes(weights)
      pv_cols <- rownames(fit$diagnostics$weighted$beta)
      if (!is.matrix(weights) || !is.numeric(weights) ||
          any(!is.finite(weights)) || any(weights < 0) ||
          !identical(colnames(weights), pv_cols) ||
          !identical(names(weight_attributes), c("dim", "dimnames")) ||
          !is.null(rownames(weights)) ||
          any(abs(colSums(weights) - 1) > 1e-12)) {
        pv_abort("stack_psis retained weights must use the canonical PV-aligned matrix envelope.")
      }
      if (nrow(weights) != nrow(weighted_draws)) {
        pv_abort("stack_psis retained source draws and weights must align by draw count.")
      }
      for (pv in pv_cols) {
        recomputed <- pv_weighted_mean_cov(weighted_draws, weights[, pv])
        expected_beta <- unname(fit$diagnostics$weighted$beta[pv, ])
        expected_U <- fit$diagnostics$weighted$U[[pv]]
        if (max(abs(recomputed$mean - expected_beta)) > 1e-10 ||
            max(abs(recomputed$cov - expected_U)) > 1e-10) {
          pv_abort("stack_psis retained proposal/weight payload must reproduce weighted summaries.")
        }
      }
    }
  }
  invisible(fit)
}

validate_pvstackr_design <- function(design) {
  pv_assert_named_list(design, "design")
  required <- pv_design_fields()
  pv_required_fields(design, required, "Design object")
  if (identical(design$schema_version, "0.2.0") && is.null(design$data)) {
    return(pv_validate_design_data_free_snapshot(design))
  }
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

pv_stack_fit_assert_canonical_owned_value <- function(value, path) {
  if (is.null(value)) {
    return(invisible(value))
  }
  if (is.environment(value) || is.function(value) || is.data.frame(value) ||
      is.raw(value)) {
    pv_abort(sprintf("Current stack-fit owned payload `%s` contains an unsupported retained object.", path))
  }
  if (inherits(value, "formula")) {
    formula_attributes <- attributes(value)
    if (!identical(names(formula_attributes), c("class", ".Environment")) ||
        !identical(formula_attributes$class, "formula") ||
        !(identical(formula_attributes$.Environment, baseenv()) ||
          identical(formula_attributes$.Environment, asNamespace("stats")))) {
      pv_abort(sprintf("Current stack-fit formula `%s` is not canonical.", path))
    }
    return(invisible(value))
  }
  if (is.matrix(value)) {
    matrix_attributes <- attributes(value)
    if (!identical(names(matrix_attributes), c("dim", "dimnames")) ||
        !is.list(matrix_attributes$dimnames) ||
        length(matrix_attributes$dimnames) != 2L ||
        !is.null(attributes(matrix_attributes$dimnames))) {
      pv_abort(sprintf("Current stack-fit matrix `%s` has noncanonical attributes.", path))
    }
    return(invisible(value))
  }
  if (is.list(value)) {
    value_names <- names(value)
    exact_list_schemas <- list(
      "stack_fit$diagnostics$sampler" = c(
        pv_sampler_diagnostic_required_fields(),
        "diagnostic_reason_codes"
      ),
      "stack_fit$meta$cache" = c(
        "enabled", "policy", "cache_stem", "file_refit",
        "directory_created", "writable"
      ),
      "stack_fit$meta$vc_policy" = c(
        "policy_id", "status", "calibration_status", "validation_status",
        "reporting_status", "confirmatory_reporting_allowed",
        "pipeline_consumption_allowed", "reportable_parameter_scope",
        "reportable_parameter_regex"
      ),
      "stack_fit$meta$prior_policy" = c(
        "identity_scope", "explicit_prior_supplied", "explicit_prior",
        "flatness", "explicit_prior_warning", "warn_explicit_prior",
        "reason_code", "warning"
      ),
      "stack_fit$meta$prior_diagnostic" = c(
        "identity_scope", "explicit_prior_supplied", "explicit_prior",
        "flatness", "explicit_prior_warning", "warn_explicit_prior",
        "reason_code", "warning"
      ),
      "stack_fit$provenance$engine" = c(
        "adapter_source", "adapter_id", "requested_backend",
        "resolved_backend", "engine_id", "selection_policy",
        "selection_reason", "package_versions", "cmdstan_state"
      ),
      "stack_fit$provenance$cache" = c(
        "enabled", "policy", "cache_stem", "file_refit",
        "directory_created", "writable"
      )
    )
    variable_list_schemas <- list(
      "stack_fit$diagnostics" = list(
        "sampler",
        c("custom_diagnostic_reason_codes", "sampler"),
        c("custom_sampler_override_ignored", "sampler"),
        c(
          "custom_diagnostic_reason_codes",
          "custom_sampler_override_ignored",
          "sampler"
        )
      ),
      "stack_fit$provenance$engine$package_versions" = list(
        NULL,
        c("brms", "posterior", "cmdstanr", "cmdstan", "rstan")
      ),
      "stack_fit$provenance$engine$cmdstan_state" = list(
        c(
          "namespace_available", "cmdstan_configured", "toolchain_checked",
          "state_reason"
        ),
        c(
          "namespace_available", "package_version", "cmdstan_configured",
          "cmdstan_version", "cmdstan_path_basename", "state_reason",
          "toolchain_checked"
        )
      )
    )
    root_owned_list_paths <- c(
      "stack_fit$diagnostics", "stack_fit$param_map",
      "stack_fit$weight_summary", "stack_fit$meta", "stack_fit$provenance"
    )
    nested_list_paths <- c(
      names(exact_list_schemas), names(variable_list_schemas)
    )
    if (!path %in% c(root_owned_list_paths, nested_list_paths)) {
      pv_abort(sprintf(
        "Current stack-fit owned payload `%s` contains an unsupported nested list.",
        path
      ))
    }
    expected_fields <- exact_list_schemas[[path]]
    if (!is.null(expected_fields) && !identical(value_names, expected_fields)) {
      pv_abort(sprintf(
        "Current stack-fit list `%s` does not use its exact canonical fields.",
        path
      ))
    }
    allowed_fields <- variable_list_schemas[[path]]
    if (!is.null(allowed_fields) &&
        !any(vapply(allowed_fields, identical, logical(1), value_names))) {
      pv_abort(sprintf(
        "Current stack-fit list `%s` does not use an allowed canonical schema.",
        path
      ))
    }
    expected_attributes <- if (length(value) == 0L) {
      NULL
    } else {
      list(names = value_names)
    }
    if (length(value) > 0L &&
        (is.null(value_names) || anyNA(value_names) ||
          any(!nzchar(value_names)) || anyDuplicated(value_names)) ||
        !identical(attributes(value), expected_attributes)) {
      pv_abort(sprintf("Current stack-fit list `%s` has noncanonical fields or attributes.", path))
    }
    for (index in seq_along(value)) {
      pv_stack_fit_assert_canonical_owned_value(
        value[[index]],
        paste0(path, "$", value_names[[index]])
      )
    }
    return(invisible(value))
  }
  value_attributes <- attributes(value)
  named_vector_paths <- c(
    "stack_fit$psi_hat_fe",
    "stack_fit$weight_summary$per_pv_weight_sum"
  )
  expected_attributes <- if (path %in% named_vector_paths) {
    list(names = names(value))
  } else {
    NULL
  }
  if (!identical(value_attributes, expected_attributes)) {
    pv_abort(sprintf("Current stack-fit scalar/vector `%s` has noncanonical attributes.", path))
  }
  invisible(value)
}

pv_validate_current_stack_fit_semantics <- function(stack_fit) {
  ws <- stack_fit$weight_summary
  meta <- stack_fit$meta
  provenance <- stack_fit$provenance
  engine <- provenance$engine

  if (!identical(meta$vc_policy, pv_stack_vc_policy()) ||
      !identical(meta$vc_policy_id, meta$vc_policy$policy_id) ||
      !identical(
        meta$vc_confirmatory_reporting_allowed,
        meta$vc_policy$confirmatory_reporting_allowed
      ) ||
      !identical(
        meta$reportable_parameter_scope,
        meta$vc_policy$reportable_parameter_scope
      )) {
    pv_abort("Current stack-fit variance-component policy is not canonical.")
  }
  canonical_prior_policies <- list(
    pv_stack_prior_policy(NULL),
    pv_stack_prior_policy(TRUE)
  )
  prior_is_canonical <- any(vapply(
    canonical_prior_policies,
    identical,
    logical(1),
    meta$prior_policy
  ))
  if (!prior_is_canonical ||
      !identical(meta$prior_diagnostic, meta$prior_policy)) {
    pv_abort("Current stack-fit prior policy and diagnostic tuple are not canonical.")
  }

  dropped_names <- stack_fit$param_map$dropped_names
  dropped_idx <- stack_fit$param_map$original_dropped_idx
  if (!is.character(dropped_names) || anyNA(dropped_names) ||
      any(!nzchar(dropped_names)) || anyDuplicated(dropped_names) ||
      !identical(dropped_names, pv_stack_redact_dropped_names(dropped_names)) ||
      !is.integer(dropped_idx) || anyNA(dropped_idx) ||
      length(dropped_names) != length(dropped_idx) ||
      !identical(meta$dropped_draw_columns, dropped_names) ||
      !identical(meta$param_map_source, stack_fit$param_map$map_source)) {
    pv_abort("Current stack-fit dropped-column metadata is not canonical.")
  }
  pv_cols <- provenance$pv_cols
  if (!is.character(pv_cols) || length(pv_cols) < 1L || anyNA(pv_cols) ||
      any(!nzchar(pv_cols)) || anyDuplicated(pv_cols) ||
      !identical(ws$M, as.integer(length(pv_cols))) ||
      !identical(names(ws$per_pv_weight_sum), pv_cols)) {
    pv_abort("Current stack-fit PV metadata is not canonical.")
  }
  materialized <- ws$model_matrix_materialized
  if (!identical(materialized, TRUE) && !identical(materialized, FALSE)) {
    pv_abort("Current stack-fit model-matrix materialization flag is not canonical.")
  }
  expected_hash_columns <- c(
    ".pvstackr_y", ".pvstackr_pv", ".pvstackr_row", ".pvstackr_weight",
    if (materialized) {
      unique(regmatches(
        stack_fit$formula_string,
        gregexpr(
          "pvstackr(?:MM[0-9]{3}|Offset)",
          stack_fit$formula_string,
          perl = TRUE
        )
      )[[1L]])
    } else {
      character()
    }
  )
  if (!identical(ws$long_data_hash_columns, expected_hash_columns) ||
      !identical(provenance$long_data_hash_columns, expected_hash_columns)) {
    pv_abort("Current stack-fit long-data hash columns are not canonical.")
  }

  if (!identical(ws$topology, "single_long_fit") ||
      !is.integer(ws$n_original) || length(ws$n_original) != 1L ||
      is.na(ws$n_original) || ws$n_original < 1L ||
      !is.integer(ws$n_long) || length(ws$n_long) != 1L ||
      is.na(ws$n_long) || ws$n_long < 1L ||
      !is.numeric(ws$fractional_weight) ||
      length(ws$fractional_weight) != 1L ||
      !is.finite(ws$fractional_weight) ||
      !isTRUE(all.equal(ws$fractional_weight, 1 / ws$M, tolerance = 1e-12)) ||
      !is.character(ws$weight_source) || length(ws$weight_source) != 1L ||
      is.na(ws$weight_source) || !nzchar(ws$weight_source) ||
      !identical(
        ws$weight_source,
        if (is.null(ws$weight_col)) "constant_fractional" else ws$weight_col
      ) ||
      !identical(ws$weight_col, provenance$weight_col) ||
      !is.numeric(ws$mean_long_weight) ||
      length(ws$mean_long_weight) != 1L ||
      !is.finite(ws$mean_long_weight) ||
      !isTRUE(all.equal(ws$mean_long_weight, 1 / ws$M, tolerance = 1e-12)) ||
      !is.numeric(ws$total_long_weight) ||
      length(ws$total_long_weight) != 1L ||
      !is.finite(ws$total_long_weight) ||
      !is.numeric(ws$per_pv_weight_sum) || anyNA(ws$per_pv_weight_sum) ||
      length(ws$per_pv_weight_sum) != ws$M ||
      any(!is.finite(ws$per_pv_weight_sum)) ||
      any(abs(ws$per_pv_weight_sum - ws$n_original / ws$M) > 1e-8) ||
      !is.character(ws$long_data_hash) ||
      length(ws$long_data_hash) != 1L || is.na(ws$long_data_hash) ||
      !grepl("^[0-9a-f]{8}$", ws$long_data_hash)) {
    pv_abort("Current stack-fit weight summary is not canonical.")
  }
  sha_fields <- c(
    "model_matrix_bundle_hash", "model_matrix_values_hash",
    "offset_values_hash"
  )
  if (materialized) {
    sha_values <- unlist(ws[sha_fields], use.names = FALSE)
    if (!is.character(sha_values) || length(sha_values) != length(sha_fields) ||
        anyNA(sha_values) || any(!grepl("^sha256:[0-9a-f]{64}$", sha_values))) {
      pv_abort("Current materialized stack-fit hashes are not canonical SHA-256 values.")
    }
  } else if (any(!vapply(ws[sha_fields], is.null, logical(1)))) {
    pv_abort("Current non-materialized stack-fit must not retain model-matrix hashes.")
  }

  mirror_fields <- c(
    "long_data_hash", "model_matrix_materialized",
    "model_matrix_bundle_hash", "model_matrix_values_hash",
    "offset_values_hash"
  )
  if (!all(vapply(
    mirror_fields,
    function(field) {
      identical(ws[[field]], meta[[field]]) &&
        identical(ws[[field]], provenance[[field]])
    },
    logical(1)
  )) ||
      !identical(meta$topology, "single_long_fit") ||
      !identical(meta$n_fits, 1L) ||
      !identical(meta$long_data_rows, ws$n_long) ||
      !is.logical(meta$log_lik_extracted) ||
      length(meta$log_lik_extracted) != 1L ||
      is.na(meta$log_lik_extracted) ||
      !is.logical(meta$log_lik_retained) ||
      length(meta$log_lik_retained) != 1L ||
      is.na(meta$log_lik_retained) ||
      !identical(meta$log_lik_retained, !is.null(stack_fit$log_lik)) ||
      !identical(provenance$function_name, "pv_stack_fit") ||
      !identical(provenance$package, "pvstackr") ||
      !identical(provenance$topology, "single_long_fit") ||
      !identical(provenance$method_stage, "stack_fit") ||
      !identical(provenance$long_data_hash, ws$long_data_hash) ||
      !identical(provenance$formula_string, stack_fit$formula_string) ||
      !identical(
        stack_fit$formula_string,
        paste(deparse(stack_fit$formula, width.cutoff = 500L), collapse = "")
      ) ||
      !is.logical(provenance$backend_fit_retained) ||
      length(provenance$backend_fit_retained) != 1L ||
      is.na(provenance$backend_fit_retained) ||
      !is.logical(provenance$log_lik_retained) ||
      length(provenance$log_lik_retained) != 1L ||
      is.na(provenance$log_lik_retained) ||
      !identical(provenance$backend_fit_retained, !is.null(stack_fit$fit)) ||
      !identical(provenance$log_lik_retained, !is.null(stack_fit$log_lik)) ||
      !identical(provenance$schema_version, "0.2.0")) {
    pv_abort("Current stack-fit metadata and provenance mirrors are not canonical.")
  }

  if (identical(engine$adapter_source, "injected")) {
    expected_state <- list(
      namespace_available = NA,
      cmdstan_configured = NA,
      toolchain_checked = FALSE,
      state_reason = "not_evaluated_for_injected_adapter"
    )
    if (!identical(engine$package_versions, list()) ||
        !identical(engine$cmdstan_state, expected_state)) {
      pv_abort("Injected stack-fit engine provenance must use the exact redacted tuple.")
    }
  } else if (identical(engine$adapter_source, "bundled")) {
    versions <- engine$package_versions
    valid_version <- function(value) {
      is.character(value) && length(value) == 1L &&
        (is.na(value) || grepl("^[0-9]+([.-][0-9A-Za-z]+)*$", value))
    }
    if (any(!vapply(versions, valid_version, logical(1)))) {
      pv_abort("Bundled stack-fit package versions are not canonical scalars.")
    }
  }
  invisible(stack_fit)
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
  if (identical(stack_fit$schema_version, "0.2.0")) {
    current_fields <- c(
      "stacked_draws", "diagnostics", "log_lik", "psi_hat_fe",
      "param_map", "formula", "formula_string", "weight_summary", "meta",
      "fit", "prepared_data", "control", "schema_version", "provenance",
      "warnings"
    )
    root_attributes <- attributes(stack_fit)
    if (!identical(names(stack_fit), current_fields) ||
        !identical(names(root_attributes), c("names", "class")) ||
        !identical(root_attributes$names, current_fields) ||
        !identical(
          root_attributes$class,
          c("pvstackr_stack_fit", "list")
        )) {
      pv_abort("Current stack-fit root fields, order, class, and attributes must be exact.")
    }
    nested_fields <- list(
      param_map = c(
        "fe_idx", "vc_idx", "fe_names", "vc_names", "original_fe_idx",
        "original_vc_idx", "original_dropped_idx", "dropped_names",
        "map_source"
      ),
      weight_summary = c(
        "topology", "M", "n_original", "n_long", "fractional_weight",
        "weight_source", "weight_col", "mean_long_weight",
        "total_long_weight", "per_pv_weight_sum", "long_data_hash",
        "long_data_hash_columns", "model_matrix_materialized",
        "model_matrix_bundle_hash", "model_matrix_values_hash",
        "offset_values_hash"
      ),
      meta = c(
        "topology", "engine_id", "fit_engine", "adapter_source",
        "resolved_backend", "backend_selection_reason", "cache", "n_fits",
        "long_data_rows", "long_data_hash", "log_lik_extracted",
        "log_lik_retained", "vc_policy", "vc_policy_id",
        "vc_confirmatory_reporting_allowed", "reportable_parameter_scope",
        "prior_policy", "prior_diagnostic", "param_map_source",
        "dropped_draw_columns", "model_matrix_materialized",
        "model_matrix_bundle_hash", "model_matrix_values_hash",
        "offset_values_hash"
      ),
      provenance = c(
        "function_name", "package", "topology", "method_stage", "pv_cols",
        "weight_col", "long_data_hash", "long_data_hash_columns",
        "formula_string", "model_matrix_materialized",
        "model_matrix_bundle_hash", "model_matrix_values_hash",
        "offset_values_hash", "backend", "engine", "cache",
        "backend_fit_retained", "log_lik_retained", "schema_version"
      )
    )
    for (field in names(nested_fields)) {
      value <- stack_fit[[field]]
      if (!is.list(value) || !identical(names(value), nested_fields[[field]]) ||
          !identical(attributes(value), list(names = nested_fields[[field]]))) {
        pv_abort(sprintf("Current stack-fit `%s` schema must be exact.", field))
      }
    }
    draws_retained <- !is.null(stack_fit$stacked_draws)
    if (!identical(stack_fit$control$return_draws, draws_retained)) {
      pv_abort("Current stack-fit draw retention must match `control$return_draws`.")
    }
    draws_attributes <- attributes(stack_fit$stacked_draws)
    if ((draws_retained &&
          (!identical(names(draws_attributes), c("dim", "dimnames")) ||
            !is.null(dimnames(stack_fit$stacked_draws)[[1L]]))) ||
        !identical(stack_fit$warnings, unname(as.character(stack_fit$warnings))) ||
        !is.null(attributes(stack_fit$warnings))) {
      pv_abort("Current stack-fit draws and warnings must use canonical bare attributes.")
    }
    if (!is.null(stack_fit$log_lik)) {
      log_lik_attributes <- attributes(stack_fit$log_lik)
      if (!identical(names(log_lik_attributes), c("dim", "dimnames")) ||
          !identical(dimnames(stack_fit$log_lik), list(NULL, NULL))) {
        pv_abort("Current stack-fit log-likelihood matrix must use canonical bare attributes.")
      }
    }
    owned_fields <- setdiff(
      current_fields,
      c("fit", "prepared_data", "control")
    )
    for (field in owned_fields) {
      pv_stack_fit_assert_canonical_owned_value(
        stack_fit[[field]],
        paste0("stack_fit$", field)
      )
    }
    pv_validate_current_stack_fit_semantics(stack_fit)
  }
  draws_retained <- !is.null(stack_fit$stacked_draws)
  draws <- if (draws_retained) {
    ccc_as_draw_matrix(stack_fit$stacked_draws)
  } else {
    NULL
  }
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
  draw_names <- if (draws_retained) {
    colnames(draws)
  } else {
    c(stack_fit$param_map$fe_names, stack_fit$param_map$vc_names)
  }
  n_params <- length(draw_names)
  if (length(fe_idx) == 0L || any(fe_idx < 1L | fe_idx > n_params)) {
    pv_abort("Stack-fit fixed-effect indexes are invalid.")
  }
  if (length(vc_idx) > 0L && any(vc_idx < 1L | vc_idx > n_params)) {
    pv_abort("Stack-fit variance-component indexes are invalid.")
  }
  if (length(intersect(fe_idx, vc_idx)) > 0L) {
    pv_abort("Stack-fit fixed-effect and variance-component indexes must be disjoint.")
  }
  if (!identical(stack_fit$param_map$fe_names, draw_names[fe_idx])) {
    pv_abort("Stack-fit fixed-effect names must match selected draw columns.")
  }
  if (!identical(stack_fit$param_map$vc_names, draw_names[vc_idx])) {
    pv_abort("Stack-fit variance-component names must match selected draw columns.")
  }
  if (!draws_retained &&
      (!identical(fe_idx, seq_along(stack_fit$param_map$fe_names)) ||
        !identical(
          vc_idx,
          if (length(stack_fit$param_map$vc_names) > 0L) {
            seq.int(length(fe_idx) + 1L, n_params)
          } else {
            integer()
          }
        ))) {
    pv_abort("Draw-free stack-fit parameter indexes must use the canonical selected-column order.")
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
  if (draws_retained) {
    psi_hat_expected <- colMeans(draws[, fe_idx, drop = FALSE])
    if (max(abs(stack_fit$psi_hat_fe - psi_hat_expected)) > 1e-10) {
      pv_abort("Stack-fit `psi_hat_fe` must equal the fixed-effect draw means.")
    }
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
    if (draws_retained) {
      pv_validate_log_lik(
        stack_fit$log_lik,
        nrow(draws),
        stack_fit$weight_summary$n_long
      )
    } else if (!is.matrix(stack_fit$log_lik) ||
        !is.numeric(stack_fit$log_lik) ||
        any(!is.finite(stack_fit$log_lik)) ||
        nrow(stack_fit$log_lik) < 1L ||
        ncol(stack_fit$log_lik) != stack_fit$weight_summary$n_long) {
      pv_abort("Draw-free stack-fit log-likelihood dimensions are invalid.")
    }
  }
  if (!is.null(stack_fit$meta$long_data_hash) &&
      !identical(stack_fit$meta$long_data_hash, stack_fit$weight_summary$long_data_hash)) {
    pv_abort("Stack-fit `meta$long_data_hash` must match `weight_summary$long_data_hash`.")
  }
  pv_validate_named_list_field(stack_fit$diagnostics, "diagnostics")
  diagnostic_fields <- c(
    "custom_diagnostic_reason_codes", "custom_sampler_override_ignored",
    "sampler"
  )
  if (any(!names(stack_fit$diagnostics) %in% diagnostic_fields) ||
      !identical(
        names(stack_fit$diagnostics),
        diagnostic_fields[diagnostic_fields %in% names(stack_fit$diagnostics)]
      ) || !"sampler" %in% names(stack_fit$diagnostics)) {
    pv_abort("Stack-fit diagnostics must use the strict normalized retention envelope.")
  }
  if (!is.null(stack_fit$diagnostics$custom_diagnostic_reason_codes) &&
      !identical(
        stack_fit$diagnostics$custom_diagnostic_reason_codes,
        "diagnostic_extraction_failed"
      )) {
    pv_abort("Stack-fit custom diagnostic reason is not canonical.")
  }
  if (!is.null(stack_fit$diagnostics$custom_sampler_override_ignored) &&
      !identical(
        stack_fit$diagnostics$custom_sampler_override_ignored,
        TRUE
      )) {
    pv_abort("Stack-fit custom sampler override flag is not canonical.")
  }
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
  if (identical(stack_fit$schema_version, "0.2.0")) {
    if (is.null(stack_fit$diagnostics$sampler)) {
      pv_abort("Current stack-fit objects require normalized sampler diagnostics.")
    }
    pv_validate_sampler_diagnostics(stack_fit$diagnostics$sampler)
  } else if (!is.null(stack_fit$diagnostics$sampler)) {
    pv_validate_sampler_diagnostics(stack_fit$diagnostics$sampler)
  }
  pv_validate_named_list_field(stack_fit$provenance, "provenance")
  engine <- stack_fit$provenance$engine
  cache <- stack_fit$provenance$cache
  legacy_provenance <- identical(stack_fit$schema_version, "0.1.0") &&
    is.null(engine) && is.null(cache) &&
    is.null(stack_fit$meta$adapter_source) &&
    is.null(stack_fit$meta$resolved_backend) &&
    is.null(stack_fit$meta$backend_selection_reason) &&
    is.null(stack_fit$meta$cache)
  if (!legacy_provenance && (is.null(engine) || is.null(cache))) {
    pv_abort("Current stack-fit objects require both backend and cache provenance.")
  }
  if (!is.null(engine)) {
    pv_assert_named_list(engine, "provenance$engine")
    engine_required <- c(
      "adapter_source", "adapter_id", "requested_backend",
      "resolved_backend", "engine_id", "selection_policy",
      "selection_reason", "package_versions", "cmdstan_state"
    )
    engine_missing <- setdiff(engine_required, names(engine))
    if (length(engine_missing) > 0L) {
      pv_abort(sprintf(
        "Stack-fit backend provenance is missing field(s): %s.",
        paste(engine_missing, collapse = ", ")
      ))
    }
    for (field in c(
      "adapter_source", "adapter_id", "requested_backend",
      "resolved_backend", "engine_id", "selection_policy",
      "selection_reason"
    )) {
      pv_assert_scalar_string(engine[[field]], paste0("provenance$engine$", field))
    }
    pv_assert_named_list(engine$package_versions, "provenance$engine$package_versions")
    pv_assert_named_list(engine$cmdstan_state, "provenance$engine$cmdstan_state")
    if (!engine$adapter_source %in% c("bundled", "injected")) {
      pv_abort("Stack-fit backend provenance has an invalid `adapter_source`.")
    }
    if (!identical(stack_fit$meta$engine_id, engine$engine_id) ||
        !identical(stack_fit$meta$fit_engine, engine$adapter_id) ||
        !identical(stack_fit$meta$adapter_source, engine$adapter_source) ||
        !identical(stack_fit$meta$resolved_backend, engine$resolved_backend) ||
        !identical(stack_fit$meta$backend_selection_reason, engine$selection_reason)) {
      pv_abort("Stack-fit backend metadata must match backend provenance.")
    }
    if (!identical(engine$requested_backend, stack_fit$control$backend) ||
        !identical(engine$requested_backend, stack_fit$provenance$backend)) {
      pv_abort("Stack-fit requested backend must match control and provenance metadata.")
    }
    if (identical(engine$adapter_source, "bundled") &&
        (!identical(engine$adapter_id, "bundled_brms") ||
         !identical(engine$requested_backend, "brms") ||
         !engine$resolved_backend %in% c("cmdstanr", "rstan") ||
         !identical(engine$engine_id, paste0("bundled_brms_", engine$resolved_backend)) ||
         !identical(
           engine$selection_policy,
           "cmdstanr_when_namespace_and_cmdstan_configured_else_rstan"
         ))) {
      pv_abort("Bundled brms engine provenance is inconsistent.")
    }
    if (identical(engine$adapter_source, "bundled")) {
      state_required <- c(
        "namespace_available", "package_version", "cmdstan_configured",
        "cmdstan_version", "cmdstan_path_basename", "state_reason",
        "toolchain_checked"
      )
      if (!all(state_required %in% names(engine$cmdstan_state))) {
        pv_abort("Bundled brms CmdStan provenance is incomplete.")
      }
      state <- engine$cmdstan_state
      if (!is.logical(state$namespace_available) ||
          length(state$namespace_available) != 1L ||
          is.na(state$namespace_available) ||
          !is.logical(state$cmdstan_configured) ||
          length(state$cmdstan_configured) != 1L ||
          is.na(state$cmdstan_configured) ||
          !identical(state$toolchain_checked, FALSE)) {
        pv_abort("Bundled brms CmdStan provenance has invalid state fields.")
      }
      if (identical(engine$resolved_backend, "cmdstanr") &&
          (!isTRUE(state$namespace_available) ||
           !isTRUE(state$cmdstan_configured) ||
           !identical(engine$selection_reason, "configured_cmdstan_selected") ||
           !is.character(state$package_version) ||
           length(state$package_version) != 1L || is.na(state$package_version) ||
           !nzchar(state$package_version) ||
           !is.character(state$cmdstan_version) ||
           length(state$cmdstan_version) != 1L || is.na(state$cmdstan_version) ||
           !nzchar(state$cmdstan_version) ||
           !is.character(state$cmdstan_path_basename) ||
           length(state$cmdstan_path_basename) != 1L ||
           is.na(state$cmdstan_path_basename) ||
           !nzchar(state$cmdstan_path_basename) ||
           !identical(basename(state$cmdstan_path_basename),
                      state$cmdstan_path_basename) ||
           !identical(
             state$state_reason,
             "cmdstanr_namespace_and_cmdstan_configured"
           ))) {
        pv_abort("Bundled cmdstanr provenance requires a configured CmdStan state.")
      }
      if (identical(engine$resolved_backend, "rstan") &&
          (isTRUE(state$namespace_available) ||
           isTRUE(state$cmdstan_configured) ||
           !identical(
             engine$selection_reason,
             "cmdstanr_namespace_absent_rstan_selected"
           ) ||
           !identical(state$package_version, NA_character_) ||
           !identical(state$cmdstan_version, NA_character_) ||
           !identical(state$cmdstan_path_basename, NA_character_) ||
           !identical(state$state_reason, "cmdstanr_namespace_unavailable"))) {
        pv_abort("Bundled rstan provenance requires cmdstanr namespace absence.")
      }
    }
    if (identical(engine$adapter_source, "injected") &&
        (!identical(engine$adapter_id, "injected_fit_function") ||
         !identical(engine$resolved_backend, "injected") ||
         !identical(engine$engine_id, "injected_fit_function") ||
         !identical(engine$selection_policy, "caller_supplied_fit_function") ||
         !identical(engine$selection_reason, "fit_function_supplied") ||
         !identical(engine$cmdstan_state$toolchain_checked, FALSE) ||
         !identical(
           engine$cmdstan_state$state_reason,
           "not_evaluated_for_injected_adapter"
         ))) {
      pv_abort("Injected engine provenance is inconsistent.")
    }
  }
  if (!is.null(cache)) {
    pv_assert_named_list(cache, "provenance$cache")
    cache_required <- c(
      "enabled", "policy", "cache_stem", "file_refit",
      "directory_created", "writable"
    )
    cache_missing <- setdiff(cache_required, names(cache))
    if (length(cache_missing) > 0L) {
      pv_abort(sprintf(
        "Stack-fit cache provenance is missing field(s): %s.",
        paste(cache_missing, collapse = ", ")
      ))
    }
    if (!is.logical(cache$enabled) || length(cache$enabled) != 1L ||
        is.na(cache$enabled)) {
      pv_abort("Stack-fit cache provenance is invalid.")
    }
    pv_assert_scalar_string(cache$policy, "provenance$cache$policy")
    if (!cache$policy %in% c(
      "disabled", "bundled_brms_managed", "injected_adapter_managed"
    )) {
      pv_abort("Stack-fit cache provenance is invalid.")
    }
    if (!identical(stack_fit$meta$cache, cache)) {
      pv_abort("Stack-fit cache metadata must match cache provenance.")
    }
    if (!cache$enabled &&
        (!identical(cache$policy, "disabled") ||
         !identical(cache$file_refit, "never") ||
         !identical(cache$cache_stem, NA_character_) ||
         !identical(cache$directory_created, FALSE) ||
         !is.logical(cache$writable) || length(cache$writable) != 1L ||
         !is.na(cache$writable))) {
      pv_abort("Disabled stack-fit cache provenance has an inconsistent tuple.")
    }
    if (isTRUE(cache$enabled)) {
      pv_assert_scalar_string(cache$cache_stem, "provenance$cache$cache_stem")
      if (!nzchar(cache$cache_stem) || cache$cache_stem %in% c(".", "..") ||
          !identical(basename(cache$cache_stem), cache$cache_stem) ||
          grepl("[/\\\\]", cache$cache_stem) ||
          !identical(cache$file_refit, "on_change") ||
          !is.logical(cache$directory_created) ||
          length(cache$directory_created) != 1L ||
          is.na(cache$directory_created)) {
        pv_abort("Enabled stack-fit cache provenance has an inconsistent tuple.")
      }
    }
    if (isTRUE(cache$enabled) &&
        identical(engine$adapter_source, "bundled") &&
        (!identical(cache$policy, "bundled_brms_managed") ||
         !identical(cache$writable, TRUE))) {
      pv_abort("Bundled brms cache provenance must use the package-managed policy.")
    }
    if (isTRUE(cache$enabled) &&
        identical(engine$adapter_source, "injected") &&
        (!identical(cache$policy, "injected_adapter_managed") ||
         !identical(cache$directory_created, FALSE) ||
         !is.logical(cache$writable) || length(cache$writable) != 1L ||
         !is.na(cache$writable))) {
      pv_abort("Injected adapter cache provenance must use the adapter-managed policy.")
    }
  }
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
  if (identical(stack_fit$schema_version, "0.2.0")) {
    expected_warnings <- c(
      stack_fit$meta$prior_policy$warning,
      pv_stack_param_drop_warning(
        stack_fit$param_map$dropped_names,
        stack_fit$param_map$map_source
      )
    )
    if (!identical(stack_fit$warnings, expected_warnings)) {
      pv_abort("Current stack-fit warnings must be exactly reconstructed from canonical diagnostics.")
    }
  }
  invisible(stack_fit)
}

validate_pvstackr_ccc <- function(ccc) {
  pv_assert_named_list(ccc, "ccc")
  schema_version <- ccc$schema_version
  if (identical(schema_version, "0.2.0")) {
    fields <- pv_ccc_schema02_fields()
    root_attributes <- attributes(ccc)
    if (!identical(names(ccc), fields) ||
        !identical(names(root_attributes), c("names", "class")) ||
        !identical(root_attributes$names, fields) ||
        !identical(root_attributes$class, c("pvstackr_ccc", "list"))) {
      pv_abort("Schema-0.2 CCC root fields, order, and class must be exact.")
    }
    provenance_fields <- c("function_name", "package", "draws_retained")
    if (!is.list(ccc$provenance) ||
        !identical(names(ccc$provenance), provenance_fields) ||
        !identical(
          attributes(ccc$provenance),
          list(names = provenance_fields)
        ) ||
        !identical(ccc$provenance$function_name, "ccc_calibrate") ||
        !identical(ccc$provenance$package, "pvstackr") ||
        (!identical(ccc$provenance$draws_retained, TRUE) &&
          !identical(ccc$provenance$draws_retained, FALSE))) {
      pv_abort("Schema-0.2 CCC provenance and draw-retention flag must be exact.")
    }
    draws_retained <- ccc$provenance$draws_retained
    if (!identical(!is.null(ccc$draws_calibrated), draws_retained) ||
        !identical(!is.null(ccc$draws_fe_cal), draws_retained)) {
      pv_abort("Schema-0.2 CCC draw payloads must match its retention flag.")
    }
    pv_binding_proof_validate(ccc$binding_proof)
    if (!is.character(ccc$target_hash) || length(ccc$target_hash) != 1L ||
        is.na(ccc$target_hash) ||
        !grepl("^sha256:[0-9a-f]{64}$", ccc$target_hash)) {
      pv_binding_abort(
        "PV_BIND_E090",
        "Schema-0.2 CCC target_hash must be a canonical target-content SHA-256 link.",
        "target_content",
        observed_hash = ccc$target_hash
      )
    }
  } else if (identical(schema_version, "0.1.0")) {
    required <- setdiff(pv_ccc_schema01_fields(), "Sigma_cal_emp_raw")
    missing <- setdiff(required, names(ccc))
    if (length(missing) > 0L) {
      pv_abort(sprintf("CCC object is missing required field(s): %s.", paste(missing, collapse = ", ")))
    }
    if ("binding_proof" %in% names(ccc)) {
      pv_abort("Schema-0.1 diagnostic CCC objects cannot carry a binding proof.")
    }
    draws_retained <- !is.null(ccc$draws_calibrated) &&
      !is.null(ccc$draws_fe_cal)
  } else {
    pv_abort("CCC `schema_version` must be 0.1.0 or 0.2.0.")
  }
  pv_assert_named_list(ccc$param_map, "param_map")
  fe_idx <- as.integer(ccc$param_map$fe_idx)
  vc_idx <- as.integer(ccc$param_map$vc_idx)
  fe_names <- ccc$param_map$fe_names
  vc_names <- ccc$param_map$vc_names
  if (!is.character(fe_names) || length(fe_names) < 1L || anyNA(fe_names) ||
      any(!nzchar(fe_names)) || anyDuplicated(fe_names) ||
      !is.character(vc_names) || anyNA(vc_names) ||
      any(!nzchar(vc_names)) || anyDuplicated(vc_names) ||
      length(intersect(fe_names, vc_names)) > 0L) {
    pv_abort("CCC parameter names must be unique canonical fixed-effect and nuisance labels.")
  }
  n_all <- length(fe_names) + length(vc_names)
  if (length(fe_idx) == 0L || any(fe_idx < 1L | fe_idx > n_all)) {
    pv_abort("CCC fixed-effect indexes are invalid.")
  }
  if (length(vc_idx) > 0L && any(vc_idx < 1L | vc_idx > n_all)) {
    pv_abort("CCC variance-component indexes are invalid.")
  }
  if (length(intersect(fe_idx, vc_idx)) > 0L) {
    pv_abort("CCC fixed-effect and variance-component indexes must be disjoint.")
  }
  if (length(fe_idx) != length(fe_names) ||
      length(vc_idx) != length(vc_names) ||
      !identical(sort(c(fe_idx, vc_idx)), seq_len(n_all))) {
    pv_abort("CCC parameter indexes must cover every selected draw column exactly once.")
  }
  all_names <- character(n_all)
  all_names[fe_idx] <- fe_names
  all_names[vc_idx] <- vc_names
  if (draws_retained) {
    draws_all <- ccc_as_draw_matrix(ccc$draws_calibrated)
    draws_fe <- ccc_as_draw_matrix(ccc$draws_fe_cal)
    if (!identical(colnames(draws_all), all_names) ||
        !identical(colnames(draws_fe), fe_names) ||
        !identical(draws_all[, fe_idx, drop = FALSE], draws_fe)) {
      pv_abort("CCC retained draw matrices must align exactly with the canonical parameter map.")
    }
  }
  if (!is.numeric(ccc$psi_hat) || any(!is.finite(ccc$psi_hat)) ||
      !identical(names(ccc$psi_hat), fe_names)) {
    pv_abort("CCC `psi_hat` must be finite, named, and aligned to fixed-effect draws.")
  }
  if (!is.numeric(ccc$psi_raw) || any(!is.finite(ccc$psi_raw)) ||
      !identical(names(ccc$psi_raw), fe_names)) {
    pv_abort("CCC `psi_raw` must be finite, named, and aligned to fixed-effect draws.")
  }
  if (!is.numeric(ccc$psi_target) || any(!is.finite(ccc$psi_target)) ||
      !identical(names(ccc$psi_target), fe_names)) {
    pv_abort("CCC `psi_target` must be finite, named, and aligned to fixed-effect draws.")
  }
  if (!is.matrix(ccc$A) || !is.numeric(ccc$A) || any(!is.finite(ccc$A)) ||
      !identical(dim(ccc$A), c(length(fe_names), length(fe_names)))) {
    pv_abort("CCC `A` must be a finite square matrix aligned to fixed-effect draws.")
  }
  if (!is.matrix(ccc$A_full) || !identical(dim(ccc$A_full), c(n_all, n_all))) {
    pv_abort("CCC `A_full` dimensions must match all calibrated draw columns.")
  }
  if (max(abs(ccc$A_full[fe_idx, fe_idx, drop = FALSE] - ccc$A)) > 1e-12) {
    pv_abort("CCC `A_full` fixed-effect block must equal `A`.")
  }
  outside_idx <- setdiff(seq_len(n_all), fe_idx)
  if (length(outside_idx) > 0L) {
    identity_outside <- diag(n_all)[outside_idx, outside_idx, drop = FALSE]
    if (max(abs(ccc$A_full[outside_idx, outside_idx, drop = FALSE] - identity_outside)) > 1e-12 ||
        max(abs(ccc$A_full[outside_idx, fe_idx, drop = FALSE])) > 1e-12 ||
        max(abs(ccc$A_full[fe_idx, outside_idx, drop = FALSE])) > 1e-12) {
      pv_abort("CCC `A_full` must be identity outside fixed-effect rows and columns.")
    }
  }
  pv_validate_target_matrix(ccc$Sigma_target, fe_names)
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

pv_blocked_forbidden_diagnostic_fields <- function() {
  c(
    "fit", "backend_fit", "backend_fits", "data", "prepared_data",
    "draws", "stacked_draws", "draws_calibrated", "draws_fe_cal",
    "proposal_draws", "per_pv_draws", "weights", "weighted", "log_lik",
    "pooling", "beta", "U", "U_bar", "B", "T_MI", "lambda", "se",
    "std.error", "df", "df_complete", "A", "A_full", "psi_hat",
    "psi_raw", "psi_target", "Sigma_raw", "Sigma_target",
    "Sigma_cal_emp", "Sigma_cal_emp_raw", "param_map"
  )
}

pv_validate_blocked_diagnostic_tree <- function(x, path = "diagnostics") {
  if (is.environment(x) || is.function(x) ||
      typeof(x) %in% c("externalptr", "weakref")) {
    pv_abort(sprintf("Blocked fit `%s` contains a nonportable retained object.", path))
  }
  if (is.data.frame(x) || is.matrix(x) || is.array(x)) {
    pv_abort(sprintf("Blocked fit `%s` contains a forbidden tabular or array payload.", path))
  }
  if (is.raw(x) || is.complex(x)) {
    pv_abort(sprintf("Blocked fit `%s` contains a forbidden binary or complex payload.", path))
  }
  if (is.list(x)) {
    fields <- names(x)
    expected_list_attributes <- if (identical(path, "diagnostics$preflight")) {
      list(
        names = fields,
        class = c("pvstackr_stack_direct_preflight_snapshot", "list")
      )
    } else if (length(x) == 0L) {
      NULL
    } else {
      list(names = fields)
    }
    if (length(x) > 0L &&
        (is.null(fields) || anyNA(fields) || any(!nzchar(fields)) ||
          anyDuplicated(fields) ||
          !identical(attributes(x), expected_list_attributes))) {
      pv_abort(sprintf("Blocked fit `%s` must use an exact named-list envelope.", path))
    }
    forbidden <- intersect(fields, pv_blocked_forbidden_diagnostic_fields())
    if (length(forbidden) > 0L) {
      pv_abort(sprintf(
        "Blocked fit diagnostics retain forbidden result field(s): %s.",
        paste(forbidden, collapse = ", ")
      ))
    }
    for (field in fields) {
      pv_validate_blocked_diagnostic_tree(
        x[[field]],
        paste0(path, "$", field)
      )
    }
    return(invisible(x))
  }

  value_attributes <- attributes(x)
  if (!is.null(value_attributes)) {
    allowed_named_atomic_paths <- c(
      "diagnostics$psis$pareto_k",
      "diagnostics$psis$weight_ess_iid",
      "diagnostics$psis$weight_ess_fraction",
      "diagnostics$psis$max_normalized_weight"
    )
    named_atomic <- path %in% allowed_named_atomic_paths &&
      identical(names(value_attributes), "names") &&
      is.character(value_attributes$names) &&
      length(value_attributes$names) == length(x) &&
      !anyNA(value_attributes$names) &&
      all(nzchar(value_attributes$names)) &&
      !anyDuplicated(value_attributes$names)
    if (!named_atomic) {
      pv_abort(sprintf("Blocked fit `%s` has noncanonical or hidden attributes.", path))
    }
  }
  if ((is.numeric(x) || is.logical(x)) && length(x) > 1L &&
      !path %in% c(
        "diagnostics$psis$pareto_k",
        "diagnostics$psis$weight_ess_iid",
        "diagnostics$psis$weight_ess_fraction",
        "diagnostics$psis$max_normalized_weight"
      )) {
    pv_abort(sprintf("Blocked fit `%s` contains a forbidden numeric result vector.", path))
  }
  invisible(x)
}

pv_validate_blocked_fit_redaction <- function(fit) {
  if (!identical(fit$status, "blocked")) {
    return(invisible(fit))
  }
  if (!identical(fit$estimates, data.frame())) {
    pv_abort("Every blocked fit must use the canonical empty estimates table.")
  }
  if (identical(fit$method, "per_pv")) {
    pv_abort("Current per_pv fits have no typed blocked-object path and cannot be relabeled blocked.")
  }
  forbidden_components <- c("design", "stack_fit", "ccc", "draws")
  retained <- forbidden_components[!vapply(
    fit[forbidden_components],
    is.null,
    logical(1)
  )]
  if (length(retained) > 0L) {
    pv_abort(sprintf(
      "Every blocked fit must redact component(s): %s.",
      paste(retained, collapse = ", ")
    ))
  }
  diagnostic_names <- names(fit$diagnostics) %||% character()
  diagnostic_variants <- switch(
    fit$method,
    stack_direct = list(
      character(),
      c("preflight", "sampler", "sampler_gate", "redaction"),
      c("preflight", "sampler", "sampler_gate", "ccc", "redaction")
    ),
    per_pv = list(character()),
    stack_psis = list(c("psis", "redaction"))
  )
  if (!any(vapply(
    diagnostic_variants,
    identical,
    logical(1),
    diagnostic_names
  ))) {
    pv_abort("Blocked fit diagnostics do not match an exact current slim variant.")
  }
  if (identical(fit$method, "stack_direct") && length(diagnostic_names) == 0L) {
    if (!identical(fit$reason_codes, "preflight_failed") ||
        length(fit$warnings) != 0L ||
        !identical(fit$provenance, pv_provenance("new_pvstackr_fit")) ||
        !is.null(fit$target)) {
      pv_abort("Hollow blocked stack_direct fits are restricted to the canonical preflight-failure record.")
    }
  }
  retention_flags <- c(
    "return_draws", "keep_data", "keep_backend_fit", "keep_log_lik"
  )
  if (any(!vapply(
    fit$control[retention_flags],
    identical,
    logical(1),
    FALSE
  ))) {
    pv_abort("Every blocked fit must record all heavy-retention controls as effective FALSE.")
  }

  independent_target <- identical(
    fit$provenance$independent_target_retained,
    TRUE
  )
  if (!is.null(fit$target)) {
    if (!identical(fit$method, "stack_direct") ||
        !independent_target ||
        !inherits(fit$target, "pvstackr_brr_target")) {
      pv_abort("A blocked fit may retain only an independently valid external stack_direct target.")
    }
  } else if (independent_target) {
    pv_abort("Blocked-fit independent-target provenance requires the retained external target.")
  }
  pv_validate_blocked_diagnostic_tree(fit$diagnostics)
  invisible(fit)
}

pv_blocked_scalar_number <- function(x, path, finite = TRUE) {
  valid <- is.numeric(x) && length(x) == 1L && !is.na(x) &&
    is.null(attributes(x)) && x >= 0
  if (isTRUE(finite)) {
    valid <- valid && is.finite(x)
  }
  if (!valid) {
    pv_abort(sprintf("Blocked CCC diagnostic `%s` must be a canonical nonnegative scalar.", path))
  }
  unname(as.numeric(x))
}

pv_validate_blocked_ccc_diagnostics <- function(diagnostics, target_hash) {
  root_fields <- c(
    "source", "target_hash", "prior", "center", "conditioning", "residual"
  )
  nested_fields <- list(
    prior = c("explicit_warning", "reason_code"),
    center = c(
      "status", "reason_code", "delta_c_rel", "delta_c_max",
      "warn_threshold", "block_threshold"
    ),
    conditioning = c(
      "status", "reason_code", "band", "kappa_A", "a_matrix_fro_rel",
      "warn_threshold", "block_threshold"
    ),
    residual = c("rho1", "rho2", "empirical_fro_rel")
  )
  if (!is.list(diagnostics) ||
      !identical(names(diagnostics), root_fields) ||
      !identical(attributes(diagnostics), list(names = root_fields)) ||
      !identical(diagnostics$source, "ccc_reportability_gate") ||
      !identical(diagnostics$target_hash, target_hash)) {
    pv_abort("Blocked stack_direct CCC diagnostics must use the exact slim root schema.")
  }
  for (field in names(nested_fields)) {
    value <- diagnostics[[field]]
    if (!is.list(value) ||
        !identical(names(value), nested_fields[[field]]) ||
        !identical(attributes(value), list(names = nested_fields[[field]]))) {
      pv_abort(sprintf("Blocked stack_direct CCC `%s` diagnostics are noncanonical.", field))
    }
  }

  prior <- diagnostics$prior
  if (!is.logical(prior$explicit_warning) ||
      length(prior$explicit_warning) != 1L ||
      is.na(prior$explicit_warning) ||
      !is.null(attributes(prior$explicit_warning)) ||
      !identical(
        prior$reason_code,
        if (isTRUE(prior$explicit_warning)) "explicit_prior_warning" else NA_character_
      )) {
    pv_abort("Blocked stack_direct CCC prior diagnostic is noncanonical.")
  }

  center <- diagnostics$center
  center_thresholds <- ccc_center_thresholds()
  delta_c_rel <- pv_blocked_scalar_number(center$delta_c_rel, "center$delta_c_rel")
  delta_c_max <- pv_blocked_scalar_number(center$delta_c_max, "center$delta_c_max")
  if (delta_c_max >= center_thresholds$block) {
    expected_center <- c(status = "blocked", reason_code = "center_separation_red")
  } else if (delta_c_max >= center_thresholds$warn) {
    expected_center <- c(status = "warning", reason_code = "center_separation_yellow")
  } else {
    expected_center <- c(status = "ok", reason_code = NA_character_)
  }
  if (!identical(center$status, unname(expected_center[["status"]])) ||
      !identical(center$reason_code, unname(expected_center[["reason_code"]])) ||
      !identical(center$warn_threshold, center_thresholds$warn) ||
      !identical(center$block_threshold, center_thresholds$block)) {
    pv_abort("Blocked stack_direct center diagnostic does not reproduce its gate status.")
  }

  conditioning <- diagnostics$conditioning
  conditioning_thresholds <- ccc_conditioning_thresholds()$kappa_A
  kappa_A <- pv_blocked_scalar_number(
    conditioning$kappa_A,
    "conditioning$kappa_A",
    finite = FALSE
  )
  pv_blocked_scalar_number(
    conditioning$a_matrix_fro_rel,
    "conditioning$a_matrix_fro_rel"
  )
  expected_band <- ccc_conditioning_band(kappa_A, conditioning_thresholds)
  expected_conditioning <- c(
    green = "ok", yellow = "warning", red = "blocked"
  )[[expected_band]]
  expected_conditioning_reason <- c(
    green = NA_character_,
    yellow = "ccc_conditioning_yellow",
    red = "ccc_conditioning_red"
  )[[expected_band]]
  if (!identical(conditioning$status, expected_conditioning) ||
      !identical(conditioning$reason_code, expected_conditioning_reason) ||
      !identical(conditioning$band, expected_band) ||
      !identical(conditioning$warn_threshold, conditioning_thresholds$warn) ||
      !identical(conditioning$block_threshold, conditioning_thresholds$block)) {
    pv_abort("Blocked stack_direct conditioning diagnostic does not reproduce its gate status.")
  }
  for (field in names(diagnostics$residual)) {
    pv_blocked_scalar_number(
      diagnostics$residual[[field]],
      paste0("residual$", field)
    )
  }
  if (!"blocked" %in% c(center$status, conditioning$status)) {
    pv_abort("Blocked stack_direct CCC diagnostics require at least one blocked CCC gate.")
  }
  invisible(diagnostics)
}

pv_validate_stack_direct_current_status <- function(fit) {
  if (is.null(fit$stack_fit) ||
      !identical(fit$stack_fit$schema_version, "0.2.0")) {
    return(NULL)
  }
  if (is.null(fit$ccc) ||
      is.null(fit$diagnostics$sampler) ||
      is.null(fit$diagnostics$sampler_gate) ||
      !identical(
        fit$diagnostics$sampler,
        fit$stack_fit$diagnostics$sampler
      )) {
    pv_abort("Current stack_direct fit requires aligned top-level sampler diagnostics, gate, stack fit, and CCC output.")
  }
  if (!identical(fit$diagnostics$stack_fit, fit$stack_fit$diagnostics) ||
      !identical(fit$diagnostics$stack_fit_warnings, fit$stack_fit$warnings) ||
      !identical(fit$diagnostics$ccc, fit$ccc$diagnostics)) {
    pv_abort("Current stack_direct nested diagnostics must exactly mirror their validated stack-fit and CCC sources.")
  }
  pv_validate_sampler_diagnostics(fit$diagnostics$sampler)
  pv_validate_sampler_gate(fit$diagnostics$sampler_gate)
  expected_gate <- pv_sampler_gate(
    fit$diagnostics$sampler,
    expected_chains = fit$stack_fit$control$chains,
    expected_post_warmup_draws_per_chain =
      fit$stack_fit$control$iter - fit$stack_fit$control$warmup
  )
  if (!identical(fit$diagnostics$sampler_gate, expected_gate)) {
    pv_abort("Stored sampler gate must equal the gate recomputed from diagnostics and control.")
  }
  if (identical(expected_gate$status, "blocked")) {
    pv_abort("A current stack_direct fit carrying raw stack-fit or CCC payloads cannot have a blocked sampler gate.")
  }
  expected_status <- pv_fit_direct_status(
    fit$stack_fit,
    fit$ccc,
    expected_gate
  )
  observed_status <- list(
    status = fit$status,
    reason_codes = fit$reason_codes,
    warnings = fit$warnings
  )
  if (!identical(observed_status, expected_status)) {
    pv_abort("Current stack_direct status, reason_codes, and warnings must exactly equal the tuple recomputed from sampler, prior, and CCC diagnostics.")
  }
  expected_gate
}

pv_validate_stack_direct_binding_links <- function(fit) {
  required_non_null <- c("design", "target", "stack_fit", "ccc")
  missing <- required_non_null[vapply(
    required_non_null,
    function(name) is.null(fit[[name]]),
    logical(1)
  )]
  if (length(missing) > 0L) {
    pv_abort(sprintf(
      "Reportable stack_direct fit requires non-null component(s): %s.",
      paste(missing, collapse = ", ")
    ))
  }

  stored_preflight <- fit$diagnostics$preflight
  if (!is.list(stored_preflight) ||
      !identical(names(stored_preflight), pv_stack_direct_preflight_fields()) ||
      !identical(class(stored_preflight), c("pvstackr_stack_direct_preflight", "list")) ||
      is.null(stored_preflight$binding_proof)) {
    pv_abort("Current stack_direct fit must retain its exact canonical preflight and binding proof.")
  }
  pv_binding_proof_validate(
    stored_preflight$binding_proof,
    target_manifest = fit$target$binding_manifest
  )
  if (isTRUE(fit$control$keep_data)) {
    if (!is.data.frame(fit$design$data) ||
        identical(fit$design$schema_version, "0.2.0")) {
      pv_abort("`keep_data = TRUE` requires the canonical retained-data design.")
    }
    pv_binding_retained_raw_inputs_validate(
      data = fit$design$data,
      formula = fit$design$formula,
      target = fit$target
    )
  } else {
    expected_design <- pv_design_target_data_free_snapshot(
      fit$target,
      stored_preflight$binding_proof
    )
    if (!identical(fit$design, expected_design)) {
      pv_abort("`keep_data = FALSE` requires the exact target-linked data-free design snapshot.")
    }
  }
  expected_stack_control <- fit$control
  expected_stack_control$keep_data <- FALSE
  expected_stack_control$return_draws <- FALSE
  if (!is.null(fit$stack_fit$prepared_data) ||
      !identical(fit$stack_fit$control, expected_stack_control)) {
    pv_abort("Composite stack_direct fits must keep raw data only in the documented design location.")
  }
  expected_preflight <- list(
    formula = fit$design$formula,
    formula_string = pv_formula_string(fit$design$formula),
    rhs_string = pv_deparse_expr(pv_formula_rhs_checked(fit$design$formula)),
    fe_names = fit$target$fe_names,
    target_hash = fit$target$target_hash,
    target_source = fit$target$target_source,
    policy = list(
      fixed_effects_only = fit$target$policy$fixed_effects_only,
      target_repair = fit$target$policy$target_repair
    ),
    binding_proof = stored_preflight$binding_proof
  )
  class(expected_preflight) <- c(
    "pvstackr_stack_direct_preflight", "list"
  )
  if (!identical(stored_preflight, expected_preflight)) {
    pv_abort("Current stack_direct preflight must match the retained authenticated design projection.")
  }
  proof <- stored_preflight$binding_proof
  target_matrix <- fit$target$binding_manifest$components$model_matrix
  bundle_hashes <- c(
    fit$stack_fit$weight_summary$model_matrix_bundle_hash,
    fit$stack_fit$meta$model_matrix_bundle_hash,
    fit$stack_fit$provenance$model_matrix_bundle_hash
  )
  matrix_value_hashes <- c(
    fit$stack_fit$weight_summary$model_matrix_values_hash,
    fit$stack_fit$meta$model_matrix_values_hash,
    fit$stack_fit$provenance$model_matrix_values_hash
  )
  offset_value_hashes <- c(
    fit$stack_fit$weight_summary$offset_values_hash,
    fit$stack_fit$meta$offset_values_hash,
    fit$stack_fit$provenance$offset_values_hash
  )
  internal_columns <- sprintf(
    "pvstackrMM%03d",
    seq_len(target_matrix$column_count)
  )
  expected_hash_columns <- c(
    pv_stack_reserved_cols(),
    internal_columns,
    if (isTRUE(target_matrix$has_offset)) "pvstackrOffset" else NULL
  )
  expected_backend_rhs <- paste(
    "0 +", paste(internal_columns, collapse = " + ")
  )
  if (isTRUE(target_matrix$has_offset)) {
    expected_backend_rhs <- paste0(
      expected_backend_rhs, " + offset(pvstackrOffset)"
    )
  }
  expected_backend_formula <- paste0(
    ".pvstackr_y | weights(.pvstackr_weight) ~ ", expected_backend_rhs
  )
  canonical_backend_formula <- stats::as.formula(
    expected_backend_formula,
    env = asNamespace("stats")
  )
  expected_bundle_hash <- fit$target$binding_manifest$model_bundle_hash
  long_data_hashes <- c(
    fit$stack_fit$weight_summary$long_data_hash,
    fit$stack_fit$meta$long_data_hash,
    fit$stack_fit$provenance$long_data_hash
  )
  expected_long_data_hash <- pv_binding_stack_long_data_hash(
    fit$target$binding_manifest
  )
  if (!identical(fit$stack_fit$weight_summary$model_matrix_materialized, TRUE) ||
      !identical(fit$stack_fit$meta$model_matrix_materialized, TRUE) ||
      !identical(fit$stack_fit$provenance$model_matrix_materialized, TRUE) ||
      length(bundle_hashes) != 3L || anyNA(bundle_hashes) ||
      any(!grepl("^sha256:[0-9a-f]{64}$", bundle_hashes)) ||
      !identical(unname(bundle_hashes), rep(expected_bundle_hash, 3L)) ||
      !identical(
        unname(matrix_value_hashes),
        rep(target_matrix$values_hash, 3L)
      ) ||
      !identical(
        unname(offset_value_hashes),
        rep(target_matrix$offset_values_hash, 3L)
      ) ||
      !identical(
        fit$stack_fit$weight_summary$long_data_hash_columns,
        expected_hash_columns
      ) ||
      !identical(
        fit$stack_fit$provenance$long_data_hash_columns,
        expected_hash_columns
      ) ||
      length(long_data_hashes) != 3L || anyNA(long_data_hashes) ||
      any(!grepl("^[0-9a-f]{8}$", long_data_hashes)) ||
      !identical(
        unname(long_data_hashes),
        rep(expected_long_data_hash, 3L)
      ) ||
      !identical(fit$stack_fit$formula, canonical_backend_formula) ||
      !identical(fit$stack_fit$formula_string, expected_backend_formula) ||
      !identical(
        fit$stack_fit$provenance$formula_string,
        expected_backend_formula
      ) ||
      !identical(
        paste(deparse(fit$stack_fit$formula, width.cutoff = 500L), collapse = ""),
        expected_backend_formula
      )) {
    pv_abort("Current stack_direct backend must retain the exact canonical materialization attestation.")
  }
  if (!identical(fit$ccc$schema_version, "0.2.0") ||
      is.null(fit$ccc$binding_proof) ||
      !identical(fit$ccc$binding_proof, proof)) {
    pv_abort("Current stack_direct CCC must retain the exact verified preflight binding proof.")
  }
  if (!identical(fit$target$binding_manifest$manifest_hash, proof$target_manifest_hash) ||
      !identical(fit$target$target_content$target_content_hash, fit$target$target_hash) ||
      !identical(fit$ccc$target_hash, fit$target$target_content$target_content_hash)) {
    pv_abort("Current stack_direct target, content, manifest, and CCC hashes must be exactly cross-linked.")
  }
  expected_links <- list(
    target_hash = fit$target$target_hash,
    target_manifest_hash = proof$target_manifest_hash,
    target_content_hash = fit$target$target_content$target_content_hash,
    binding_verification_policy = proof$verification_policy
  )
  for (field in names(expected_links)) {
    if (!identical(fit$design$provenance[[field]], expected_links[[field]])) {
      pv_abort(sprintf(
        "Current stack_direct design provenance `%s` must match the verified target binding.",
        field
      ))
    }
    if (!identical(fit$provenance[[field]], expected_links[[field]])) {
      pv_abort(sprintf(
        "Current stack_direct fit provenance `%s` must match the verified target binding.",
        field
      ))
    }
  }
  invisible(proof)
}

pv_validate_stack_direct_slim_ccc_blocked_fit <- function(fit) {
  if (!identical(
        names(fit$diagnostics),
        c("preflight", "sampler", "sampler_gate", "ccc", "redaction")
      )) {
    pv_abort("CCC-blocked stack_direct fit requires the exact slim diagnostic variant.")
  }
  if (is.null(fit$target) ||
      !inherits(fit$target, "pvstackr_brr_target") ||
      !identical(fit$target$target_source, "external_brr_fay_rubin") ||
      !identical(fit$provenance$independent_target_retained, TRUE) ||
      !identical(fit$target, pv_fit_direct_independent_target(fit$target))) {
    pv_abort("CCC-blocked stack_direct fit must retain only its canonical independent external target.")
  }
  if (!identical(fit$control, pv_fit_direct_blocked_control(fit$control))) {
    pv_abort("CCC-blocked stack_direct control must use the exact fail-closed schema.")
  }

  preflight <- fit$diagnostics$preflight
  canonical_preflight <- pv_fit_direct_blocked_preflight(preflight)
  preflight_fields <- c(
    "formula_string", "rhs_string", "fe_names", "target_hash",
    "target_source", "policy", "binding_proof"
  )
  expected_policy <- list(
    fixed_effects_only = fit$target$policy$fixed_effects_only,
    target_repair = fit$target$policy$target_repair
  )
  pv_binding_proof_validate(
    preflight$binding_proof,
    target_manifest = fit$target$binding_manifest
  )
  if (!identical(preflight, canonical_preflight) ||
      !identical(names(preflight), preflight_fields) ||
      !identical(preflight$formula_string, fit$target$formula_string) ||
      !identical(preflight$rhs_string, fit$target$rhs_string) ||
      !identical(preflight$fe_names, fit$target$fe_names) ||
      !identical(preflight$target_hash, fit$target$target_hash) ||
      !identical(preflight$target_source, fit$target$target_source) ||
      !identical(preflight$policy, expected_policy)) {
    pv_abort("CCC-blocked stack_direct preflight must align with the retained target.")
  }

  sampler <- fit$diagnostics$sampler
  sampler_gate <- fit$diagnostics$sampler_gate
  pv_validate_sampler_diagnostics(sampler)
  pv_validate_sampler_gate(sampler_gate)
  expected_gate <- pv_sampler_gate(
    sampler,
    expected_chains = fit$control$chains,
    expected_post_warmup_draws_per_chain = fit$control$iter - fit$control$warmup
  )
  if (!identical(sampler_gate, expected_gate) ||
      identical(expected_gate$status, "blocked")) {
    pv_abort("CCC-blocked stack_direct fit requires an exact nonblocked sampler gate.")
  }
  pv_validate_blocked_ccc_diagnostics(
    fit$diagnostics$ccc,
    fit$target$target_hash
  )
  if (!identical(
        fit$diagnostics$redaction,
        pv_fit_direct_blocked_redaction("ccc_reportability_gate")
      )) {
    pv_abort("CCC-blocked stack_direct fit requires the exact redaction record.")
  }
  expected_status <- pv_fit_direct_slim_blocked_status(
    expected_gate,
    fit$diagnostics$ccc
  )
  observed_status <- list(
    status = fit$status,
    reason_codes = fit$reason_codes,
    warnings = fit$warnings
  )
  if (!identical(observed_status, expected_status)) {
    pv_abort("CCC-blocked stack_direct status must exactly reproduce slim gate evidence.")
  }

  expected_provenance <- c(
    pv_provenance("new_pvstackr_fit"),
    list(
      wrapper_function = "pv_fit_direct",
      target_hash = fit$target$target_hash,
      target_manifest_hash = fit$target$binding_manifest$manifest_hash,
      target_content_hash = fit$target$target_content$target_content_hash,
      binding_verification_policy = pv_binding_proof_policy(),
      independent_target_retained = TRUE,
      sampler_diagnostic_source = sampler$diagnostic_source,
      ccc_target_hash = fit$target$target_hash,
      reportability_policy = "generic_blocked_fail_closed"
    )
  )
  if (!identical(fit$provenance, expected_provenance)) {
    pv_abort("CCC-blocked stack_direct provenance must align with slim target and gate evidence.")
  }
  invisible(fit)
}

pv_validate_stack_direct_reportable_fit <- function(fit) {
  if (!identical(fit$method, "stack_direct")) {
    return(invisible(fit))
  }
  if (identical(fit$status, "blocked")) {
    if (identical(
          names(fit$diagnostics),
          c("preflight", "sampler", "sampler_gate", "ccc", "redaction")
        )) {
      return(pv_validate_stack_direct_slim_ccc_blocked_fit(fit))
    }
    if (!identical(fit$control$return_draws, FALSE)) {
      pv_abort("Blocked stack_direct fits must record effective `return_draws = FALSE`.")
    }
    if (!is.null(fit$draws) ||
        (!is.null(fit$stack_fit) && !is.null(fit$stack_fit$stacked_draws)) ||
        (!is.null(fit$ccc) &&
          (!is.null(fit$ccc$draws_calibrated) ||
            !is.null(fit$ccc$draws_fe_cal) ||
            !identical(fit$ccc$provenance$draws_retained, FALSE)))) {
      pv_abort("Blocked stack_direct fits must not retain any package-owned individual draws.")
    }
    if (is.data.frame(fit$estimates) && nrow(fit$estimates) > 0L) {
      pv_abort("Blocked stack_direct fit must not include reportable estimates.")
    }
    expected_gate <- NULL
    if (is.list(fit$diagnostics) &&
        !is.null(fit$diagnostics$sampler)) {
      pv_validate_sampler_diagnostics(fit$diagnostics$sampler)
      expected_gate <- pv_sampler_gate(
        fit$diagnostics$sampler,
        expected_chains = fit$control$chains,
        expected_post_warmup_draws_per_chain =
          fit$control$iter - fit$control$warmup
      )
    }
    if (is.null(expected_gate) ||
        !identical(expected_gate$status, "blocked")) {
      if (identical(fit$provenance$independent_target_retained, TRUE) &&
          (is.null(fit$stack_fit) || is.null(fit$ccc))) {
        pv_abort("A slim independently retained target is valid only when the recomputed sampler gate is blocked.")
      }
      if (!is.null(fit$stack_fit) &&
          identical(fit$stack_fit$schema_version, "0.2.0") &&
          !is.null(fit$ccc)) {
        pv_validate_stack_direct_binding_links(fit)
      }
      pv_validate_stack_direct_current_status(fit)
      return(invisible(fit))
    }
    fit_fields <- pv_fit_required_fields()
    if (!identical(names(fit), fit_fields) ||
        !identical(class(fit), c("pvstackr_fit", "list"))) {
      pv_abort("Sampler-blocked stack_direct fit must use the exact canonical fit schema.")
    }
    if (!is.null(fit$design) || !is.null(fit$stack_fit) ||
        !is.null(fit$ccc) || !is.null(fit$draws)) {
      pv_abort("Sampler-blocked stack_direct fit must redact design, stack-fit, CCC, and draws.")
    }
    if (is.null(fit$target) ||
        !inherits(fit$target, "pvstackr_brr_target") ||
        !identical(fit$target$target_source, "external_brr_fay_rubin") ||
        !identical(fit$provenance$independent_target_retained, TRUE)) {
      pv_abort("Sampler-blocked stack_direct fit must retain its independently valid external target with provenance.")
    }
    canonical_target <- pv_fit_direct_independent_target(fit$target)
    if (!identical(fit$target, canonical_target)) {
      pv_abort("Sampler-blocked stack_direct retained target must use the exact recursive canonical snapshot schema.")
    }
    canonical_control <- pv_fit_direct_blocked_control(fit$control)
    if (!identical(fit$control, canonical_control)) {
      pv_abort("Sampler-blocked stack_direct control must use the exact canonical schema.")
    }
    if (!identical(fit$estimates, data.frame())) {
      pv_abort("Sampler-blocked stack_direct estimates must be the canonical empty data frame.")
    }
    if (!identical(names(fit$diagnostics), c(
      "preflight", "sampler", "sampler_gate", "redaction"
    ))) {
      pv_abort("Sampler-blocked stack_direct fit requires only slim sampler diagnostics.")
    }
    preflight_fields <- c(
      "formula_string", "rhs_string", "fe_names", "target_hash",
      "target_source", "policy", "binding_proof"
    )
    canonical_preflight <- pv_fit_direct_blocked_preflight(
      fit$diagnostics$preflight
    )
    expected_preflight_policy <- list(
      fixed_effects_only = fit$target$policy$fixed_effects_only,
      target_repair = fit$target$policy$target_repair
    )
    pv_binding_proof_validate(
      fit$diagnostics$preflight$binding_proof,
      target_manifest = fit$target$binding_manifest
    )
    if (!identical(fit$diagnostics$preflight, canonical_preflight) ||
        !identical(names(fit$diagnostics$preflight), preflight_fields) ||
        !identical(
          fit$diagnostics$preflight$formula_string,
          fit$target$formula_string
        ) ||
        !identical(
          fit$diagnostics$preflight$rhs_string,
          fit$target$rhs_string
        ) ||
        !identical(
          fit$diagnostics$preflight$fe_names,
          fit$target$fe_names
        ) ||
        !identical(
          fit$diagnostics$preflight$target_hash,
          fit$target$target_hash
        ) ||
        !identical(
          fit$diagnostics$preflight$target_source,
          fit$target$target_source
        ) ||
        !identical(
          fit$diagnostics$preflight$policy,
          expected_preflight_policy
        )) {
      pv_abort("Sampler-blocked stack_direct fit requires a formula-free preflight snapshot aligned with the retained target.")
    }
    if (!identical(
          fit$diagnostics$redaction,
          pv_fit_direct_blocked_redaction("sampler_gate")
        )) {
      pv_abort("Sampler-blocked stack_direct fit requires the exact redaction record.")
    }
    pv_validate_sampler_diagnostics(fit$diagnostics$sampler)
    pv_validate_sampler_gate(fit$diagnostics$sampler_gate)
    if (!identical(fit$diagnostics$sampler_gate, expected_gate) ||
        !identical(fit$reason_codes, fit$diagnostics$sampler_gate$reason_codes) ||
        !identical(fit$warnings, fit$diagnostics$sampler_gate$warnings)) {
      pv_abort("Sampler-blocked stack_direct status must match its sampler gate.")
    }
    expected_provenance <- c(
      pv_provenance("new_pvstackr_fit"),
      list(
        wrapper_function = "pv_fit_direct",
        target_hash = fit$target$target_hash,
        target_manifest_hash = fit$target$binding_manifest$manifest_hash,
        target_content_hash = fit$target$target_content$target_content_hash,
        binding_verification_policy = pv_binding_proof_policy(),
        independent_target_retained = TRUE,
        sampler_diagnostic_source =
          fit$diagnostics$sampler$diagnostic_source
      )
    )
    if (!identical(fit$provenance, expected_provenance)) {
      pv_abort("Sampler-blocked stack_direct provenance must exactly align with the retained target and sampler source.")
    }
    return(invisible(fit))
  }
  if (!fit$status %in% c("ok", "warning")) {
    return(invisible(fit))
  }

  pv_validate_stack_direct_binding_links(fit)
  if (!identical(fit$control$center, "target") ||
      !identical(fit$ccc$center, "target") ||
      !identical(fit$ccc$control$center, "target")) {
    pv_abort("Reportable stack_direct fit requires `control$center = \"target\"` and target-centered CCC output.")
  }
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
  pv_validate_stack_direct_current_status(fit)
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
  if (!identical(!is.null(fit$draws), fit$control$return_draws)) {
    pv_abort("Reportable stack_direct top-level draw retention must match `control$return_draws`.")
  }
  if (!is.null(fit$stack_fit$stacked_draws) ||
      !identical(fit$stack_fit$control$return_draws, FALSE) ||
      !is.null(fit$ccc$draws_calibrated) ||
      !is.null(fit$ccc$draws_fe_cal) ||
      !identical(fit$ccc$provenance$draws_retained, FALSE)) {
    pv_abort("Reportable stack_direct fits must keep raw and duplicate calibrated draws out of nested components.")
  }
  if (!is.null(fit$draws)) {
    if (!identical(colnames(fit$draws), fit$target$fe_names) ||
        max(abs(colMeans(fit$draws) - fit$ccc$psi_hat)) > 1e-10 ||
        max(abs(stats::cov(fit$draws) - fit$ccc$Sigma_cal_emp_raw)) > 1e-10) {
      pv_abort("Reportable stack_direct retained draws must match canonical CCC fixed-effect summaries.")
    }
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

pv_stack_psis_psis_fields <- function() {
  c(
    "pareto_k", "threshold", "pareto_k_max", "status", "bad_pv_cols",
    "reason_code", "pv_cols", "source", "pareto_k_source",
    "weight_method", "producer", "producer_version", "normalization",
    "n_draws", "ess_definition", "weight_ess_iid",
    "weight_ess_fraction", "max_normalized_weight",
    "weight_diagnostic_authority",
    "fallback_requested", "fallback_effective"
  )
}

pv_validate_stack_psis_weight_diagnostics <- function(
  psis,
  pv_names,
  retained_weights = NULL
) {
  source_allowed <- c(
    "supplied_psis_weights", "injected_psis_function",
    "self_normalized_log_ratios"
  )
  method_allowed <- c(
    "unspecified_external", "self_normalized_raw_importance",
    "caller_declared_external_psis"
  )
  if (!is.character(psis$source) || length(psis$source) != 1L ||
      is.na(psis$source) || !psis$source %in% source_allowed ||
      !is.character(psis$pareto_k_source) ||
      length(psis$pareto_k_source) != 1L ||
      is.na(psis$pareto_k_source) ||
      !psis$pareto_k_source %in% c("supplied", "injected_function_output") ||
      !is.character(psis$weight_method) ||
      length(psis$weight_method) != 1L || is.na(psis$weight_method) ||
      !psis$weight_method %in% method_allowed ||
      !is.character(psis$producer) || length(psis$producer) != 1L ||
      !is.character(psis$producer_version) ||
      length(psis$producer_version) != 1L ||
      !identical(psis$normalization, "column_sum_one") ||
      !identical(psis$ess_definition, "kish_iid_normalized_weights_v1") ||
      !is.character(psis$weight_diagnostic_authority) ||
      length(psis$weight_diagnostic_authority) != 1L ||
      is.na(psis$weight_diagnostic_authority) ||
      !psis$weight_diagnostic_authority %in% c(
        "retained_weights_recomputed",
        "owned_stamp_bounded_projection"
      ) ||
      any(!vapply(
        psis[c(
          "source", "pareto_k_source", "weight_method", "producer",
          "producer_version", "normalization", "ess_definition",
          "weight_diagnostic_authority"
        )],
        function(value) is.null(attributes(value)),
        logical(1)
      ))) {
    pv_abort("stack_psis weight provenance and diagnostic definitions are not canonical.")
  }
  expected_pareto_source <- if (identical(
    psis$source,
    "injected_psis_function"
  )) "injected_function_output" else "supplied"
  if (!identical(psis$pareto_k_source, expected_pareto_source)) {
    pv_abort("stack_psis Pareto-k source must match the declared weight route.")
  }
  if (identical(psis$source, "self_normalized_log_ratios")) {
    if (!identical(psis$weight_method, "self_normalized_raw_importance") ||
        !identical(psis$producer, "pvstackr") ||
        !identical(psis$producer_version, "0.2.0")) {
      pv_abort("Self-normalized raw log ratios must retain package-owned unsmoothed provenance.")
    }
  } else if (identical(psis$weight_method, "caller_declared_external_psis")) {
    producer_values <- c(psis$producer, psis$producer_version)
    if (anyNA(producer_values) || any(!nzchar(trimws(producer_values))) ||
        !identical(producer_values, trimws(producer_values)) ||
        nchar(psis$producer, type = "bytes") > 128L ||
        nchar(psis$producer_version, type = "bytes") > 64L ||
        any(grepl("[[:cntrl:]]", producer_values))) {
      pv_abort("Caller-declared external PSIS weights require bounded producer and version provenance.")
    }
  } else if (!identical(psis$weight_method, "unspecified_external") ||
      !identical(psis$producer, NA_character_) ||
      !identical(psis$producer_version, NA_character_)) {
    pv_abort("Unverified external weights must not carry inferred PSIS producer provenance.")
  }

  if (!is.integer(psis$n_draws) || length(psis$n_draws) != 1L ||
      is.na(psis$n_draws) || psis$n_draws < 2L ||
      !is.null(attributes(psis$n_draws))) {
    pv_abort("stack_psis weight diagnostics require a bare positive integer draw count.")
  }
  vector_fields <- c(
    "weight_ess_iid", "weight_ess_fraction", "max_normalized_weight"
  )
  for (field in vector_fields) {
    value <- psis[[field]]
    if (!is.numeric(value) || length(value) != length(pv_names) ||
        any(!is.finite(value)) ||
        !identical(attributes(value), list(names = pv_names))) {
      pv_abort("stack_psis weight diagnostics must be bare named finite vectors aligned to every PV.")
    }
  }
  n_draws <- as.numeric(psis$n_draws)
  inverse_ess <- 1 / psis$weight_ess_iid
  max_weight <- psis$max_normalized_weight
  minimum_sum_squares <- max_weight^2 +
    (1 - max_weight)^2 / (n_draws - 1)
  if (any(psis$weight_ess_iid < 1 - 1e-12) ||
      any(psis$weight_ess_iid > n_draws + 1e-12) ||
      any(psis$weight_ess_fraction < (1 / n_draws) - 1e-12) ||
      any(psis$weight_ess_fraction > 1 + 1e-12) ||
      any(psis$max_normalized_weight < (1 / n_draws) - 1e-12) ||
      any(psis$max_normalized_weight > 1 + 1e-12) ||
      any(inverse_ess > max_weight + 1e-12) ||
      any(inverse_ess < minimum_sum_squares - 1e-12) ||
      !isTRUE(all.equal(
        unname(psis$weight_ess_fraction),
        unname(psis$weight_ess_iid / n_draws),
        tolerance = 1e-14,
        check.attributes = FALSE
      ))) {
    pv_abort("stack_psis weight ESS and maximum-weight diagnostics violate their canonical bounds or algebra.")
  }
  if (!is.null(retained_weights)) {
    if (!identical(
      psis$weight_diagnostic_authority,
      "retained_weights_recomputed"
    )) {
      pv_abort("Retained stack_psis weights require recomputed diagnostic authority.")
    }
    expected <- pv_stack_psis_weight_diagnostics(retained_weights, pv_names)
    observed <- psis[names(expected)]
    if (!identical(observed, expected)) {
      pv_abort("Retained stack_psis weights must exactly reproduce weight ESS and maximum-weight diagnostics.")
    }
  } else if (!identical(
    psis$weight_diagnostic_authority,
    "owned_stamp_bounded_projection"
  )) {
    pv_abort("Redacted stack_psis weights require bounded owned-stamp diagnostic authority.")
  }
  invisible(psis)
}

pv_validate_stack_psis_reportable_fit <- function(fit) {
  if (!identical(fit$method, "stack_psis")) {
    return(invisible(fit))
  }
  if (identical(fit$status, "blocked")) {
    if (!identical(fit$control$return_draws, FALSE)) {
      pv_abort("Blocked stack_psis fits must record effective `return_draws = FALSE`.")
    }
    if (!is.list(fit$diagnostics) || !"psis" %in% names(fit$diagnostics)) {
      pv_abort("Blocked stack_psis fit requires PSIS diagnostics.")
    }
    if (is.data.frame(fit$estimates) && nrow(fit$estimates) > 0L) {
      pv_abort("Blocked stack_psis fit must not include reportable estimates.")
    }
    if (!is.null(fit$design) || !is.null(fit$target) ||
        !is.null(fit$stack_fit) || !is.null(fit$draws) || !is.null(fit$ccc)) {
      pv_abort("Blocked stack_psis fit must redact data, targets, stack/backend fits, CCC payloads, and all draws.")
    }
    if (!setequal(names(fit$diagnostics), c("psis", "redaction")) ||
        length(fit$diagnostics) != 2L) {
      pv_abort("Blocked stack_psis fit diagnostics must contain only slim PSIS and redaction records.")
    }
    redaction <- fit$diagnostics$redaction
    if (!identical(redaction, pv_stack_psis_blocked_redaction())) {
      pv_abort("Blocked stack_psis fit requires an explicit withheld redaction record.")
    }
    psis <- fit$diagnostics$psis
    pv_assert_named_list(psis, "diagnostics$psis")
    psis_required <- pv_stack_psis_psis_fields()
    if (!identical(names(psis), psis_required) ||
        !identical(attributes(psis), list(names = psis_required))) {
      pv_abort("Blocked stack_psis fit PSIS diagnostics require the exact ordered current schema.")
    }
    threshold <- pv_validate_psis_k_threshold(psis$threshold)
    if (!identical(threshold, fit$control$psis_k_threshold)) {
      pv_abort("Blocked stack_psis fit PSIS threshold must match the fit control.")
    }
    if (!is.numeric(psis$pareto_k) || length(psis$pareto_k) < 2L ||
        is.null(names(psis$pareto_k)) || any(!nzchar(names(psis$pareto_k))) ||
        anyDuplicated(names(psis$pareto_k))) {
      pv_abort("Blocked stack_psis fit Pareto-k diagnostics must be a uniquely named numeric PV vector.")
    }
    if (!is.numeric(psis$pareto_k_max) || length(psis$pareto_k_max) != 1L ||
        !is.character(psis$status) || length(psis$status) != 1L || is.na(psis$status) ||
        !is.character(psis$bad_pv_cols) || anyNA(psis$bad_pv_cols) ||
        !is.character(psis$reason_code) || length(psis$reason_code) != 1L || is.na(psis$reason_code) ||
        !is.character(psis$fallback_requested) || length(psis$fallback_requested) != 1L ||
        is.na(psis$fallback_requested) ||
        !is.character(psis$fallback_effective) || length(psis$fallback_effective) != 1L ||
        is.na(psis$fallback_effective)) {
      pv_abort("Blocked stack_psis fit slim diagnostics have invalid field types or source provenance.")
    }
    if (!is.character(psis$pv_cols) || length(psis$pv_cols) < 2L ||
        anyNA(psis$pv_cols) || any(!nzchar(psis$pv_cols)) ||
        anyDuplicated(psis$pv_cols) || !identical(names(psis$pareto_k), psis$pv_cols)) {
      pv_abort("Blocked stack_psis fit Pareto-k diagnostics must align with the complete declared PV universe.")
    }
    pv_validate_stack_psis_weight_diagnostics(psis, psis$pv_cols)
    if (!psis$status %in% c(
      "failed", "not_evaluated", "unsmoothed", "provenance_incomplete"
    )) {
      pv_abort("Blocked stack_psis fit must record a failed, unevaluated, unsmoothed, or provenance-incomplete diagnostic status.")
    }
    expected_bad <- if (identical(psis$status, "not_evaluated")) {
      names(psis$pareto_k)[!is.finite(psis$pareto_k)]
    } else if (identical(psis$status, "failed")) {
      if (any(!is.finite(psis$pareto_k))) {
        pv_abort("Failed PSIS diagnostics must be finite; use `not_evaluated` for non-finite Pareto-k values.")
      }
      names(psis$pareto_k)[psis$pareto_k >= threshold]
    } else {
      if (any(!is.finite(psis$pareto_k)) || any(psis$pareto_k >= threshold)) {
        pv_abort("Weight-provenance PSIS blocks require complete Pareto-k values below threshold.")
      }
      character()
    }
    expected_reason <- if (identical(psis$status, "not_evaluated")) {
      "psis_k_not_evaluated"
    } else if (identical(psis$status, "failed")) {
      "psis_k_too_high"
    } else if (identical(psis$status, "unsmoothed")) {
      "psis_smoothing_not_applied"
    } else {
      "psis_weight_provenance_incomplete"
    }
    expected_max <- if (identical(psis$status, "not_evaluated")) NA_real_ else max(psis$pareto_k)
    if ((psis$status %in% c("failed", "not_evaluated") && length(expected_bad) == 0L) ||
        !identical(as.character(psis$bad_pv_cols), expected_bad) ||
        !identical(psis$reason_code, expected_reason) ||
        !identical(fit$reason_codes, expected_reason) ||
        !identical(as.numeric(psis$pareto_k_max), expected_max)) {
      pv_abort("Blocked stack_psis failure metadata are internally inconsistent.")
    }
    if (!psis$fallback_requested %in% c("block", "warn") ||
        !identical(psis$fallback_effective, "block")) {
      pv_abort("Blocked stack_psis fit must record an effective fail-closed fallback.")
    }
    expected_warning <- pv_stack_psis_blocked_warning(psis)
    if (!identical(fit$warnings, expected_warning)) {
      pv_abort("Blocked stack_psis warnings must exactly reproduce the slim PSIS gate evidence.")
    }
    expected_provenance <- c(
      pv_provenance("new_pvstackr_fit"),
      list(
        wrapper_function = "pv_fit_stack_psis",
        stacked_source = fit$provenance$stacked_source,
        psis_source = psis$source,
        reportability_policy = "immutable_psis_fail_closed"
      )
    )
    if (!is.character(fit$provenance$stacked_source) ||
        length(fit$provenance$stacked_source) != 1L ||
        !fit$provenance$stacked_source %in%
          c("stack_fit", "stacked_draws", "injected_fit") ||
        !identical(fit$provenance, expected_provenance)) {
      pv_abort("Blocked stack_psis provenance must exactly reproduce source and fail-closed policy.")
    }
    return(invisible(fit))
  }
  if (identical(fit$status, "warning")) {
    pv_abort("Warning-status stack_psis fits are legacy unsafe objects and are not reportable.")
  }
  if (!identical(fit$status, "ok")) {
    return(invisible(fit))
  }
  if (!is.null(fit$target)) {
    pv_abort("Reportable stack_psis fit must not carry a formal target object.")
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
  weighted <- fit$diagnostics$weighted
  pv_assert_named_list(psis, "diagnostics$psis")
  pv_assert_named_list(pooling, "diagnostics$pooling")
  pv_assert_named_list(weighted, "diagnostics$weighted")
  weighted_fields <- c(
    "beta", "U", "proposal_draws", "weights", "param_map"
  )
  if (!identical(names(weighted), weighted_fields) ||
      !identical(attributes(weighted), list(names = weighted_fields))) {
    pv_abort("Reportable stack_psis weighted diagnostics must use the exact retention schema.")
  }
  param_map_fields <- c(
    "fe_idx", "vc_idx", "fe_names", "vc_names", "dropped_names",
    "map_source"
  )
  if (!is.matrix(weighted$beta) || !is.numeric(weighted$beta) ||
      any(!is.finite(weighted$beta)) ||
      !identical(names(attributes(weighted$beta)), c("dim", "dimnames")) ||
      !is.list(weighted$U) ||
      !identical(names(weighted$U), rownames(weighted$beta)) ||
      !identical(attributes(weighted$U), list(names = names(weighted$U))) ||
      !is.list(weighted$param_map) ||
      !identical(names(weighted$param_map), param_map_fields) ||
      !identical(
        attributes(weighted$param_map),
        list(names = param_map_fields)
      )) {
    pv_abort("Reportable stack_psis weighted summaries and parameter map are not canonical.")
  }
  fe_names_weighted <- colnames(weighted$beta)
  if (is.null(fe_names_weighted) || length(fe_names_weighted) < 1L ||
      any(!nzchar(fe_names_weighted)) || anyDuplicated(fe_names_weighted)) {
    pv_abort("Reportable stack_psis weighted fixed-effect names are not canonical.")
  }
  for (pv in names(weighted$U)) {
    U_pv <- weighted$U[[pv]]
    if (!is.matrix(U_pv) || !is.numeric(U_pv) || any(!is.finite(U_pv)) ||
        !identical(dimnames(U_pv), list(fe_names_weighted, fe_names_weighted)) ||
        !identical(names(attributes(U_pv)), c("dim", "dimnames"))) {
      pv_abort("Reportable stack_psis weighted covariance summaries are not canonical.")
    }
  }
  if (!is.integer(weighted$param_map$fe_idx) ||
      !is.integer(weighted$param_map$vc_idx) ||
      !identical(
        weighted$param_map$fe_idx,
        seq_along(weighted$param_map$fe_names)
      ) ||
      !identical(weighted$param_map$fe_names, fe_names_weighted) ||
      !identical(weighted$param_map$vc_idx, integer()) ||
      !identical(weighted$param_map$vc_names, character()) ||
      !identical(weighted$param_map$dropped_names, character()) ||
      !identical(
        weighted$param_map$map_source,
        "fixed_effect_projection"
      ) ||
      any(!vapply(weighted$param_map, function(value) {
        is.null(attributes(value))
      }, logical(1)))) {
    pv_abort("Reportable stack_psis retained parameter-map values are not canonical.")
  }
  if (!identical(names(psis), pv_stack_psis_psis_fields()) ||
      !identical(attributes(psis), list(names = pv_stack_psis_psis_fields()))) {
    pv_abort("Reportable stack_psis fit PSIS diagnostics require the exact current schema.")
  }
  if (!is.numeric(psis$pareto_k) || length(psis$pareto_k) < 2L || any(!is.finite(psis$pareto_k))) {
    pv_abort("Reportable stack_psis fit Pareto-k diagnostics must be numeric and complete.")
  }
  threshold <- pv_validate_psis_k_threshold(psis$threshold)
  if (!identical(threshold, fit$control$psis_k_threshold)) {
    pv_abort("Reportable stack_psis fit PSIS threshold must match the selected fit control threshold.")
  }
  pv_names <- if (is.matrix(weighted$beta)) rownames(weighted$beta) else NULL
  if (is.null(pv_names) || length(pv_names) < 2L || any(!nzchar(pv_names)) ||
      anyDuplicated(pv_names) || !identical(names(psis$pareto_k), pv_names) ||
      !is.list(weighted$U) || !identical(names(weighted$U), pv_names)) {
    pv_abort("Reportable stack_psis fit requires Pareto-k diagnostics aligned to every weighted PV summary.")
  }
  if (!is.null(weighted$weights) &&
      (!is.matrix(weighted$weights) || !identical(colnames(weighted$weights), pv_names))) {
    pv_abort("Reportable stack_psis retained weights must align with every declared PV.")
  }
  if (!identical(psis$status, "ok") || any(psis$pareto_k >= threshold)) {
    pv_abort("Reportable stack_psis fit requires every Pareto-k to be finite and strictly below the immutable threshold.")
  }
  if (!identical(as.numeric(psis$pareto_k_max), max(psis$pareto_k))) {
    pv_abort("Reportable stack_psis fit Pareto-k maximum must match the complete PV diagnostics.")
  }
  canonical_pareto <- stats::setNames(
    as.numeric(psis$pareto_k),
    names(psis$pareto_k)
  )
  pv_validate_stack_psis_weight_diagnostics(
    psis,
    pv_names,
    retained_weights = weighted$weights
  )
  if (!identical(psis$weight_method, "caller_declared_external_psis") ||
      !psis$source %in% c("supplied_psis_weights", "injected_psis_function")) {
    pv_abort("Reportable stack_psis fits require caller-declared external PSIS-smoothed weights.")
  }
  if (!is.character(psis$fallback_requested) ||
      length(psis$fallback_requested) != 1L ||
      is.na(psis$fallback_requested) ||
      !psis$fallback_requested %in% c("block", "warn") ||
      !identical(psis$fallback_effective, psis$fallback_requested) ||
      any(!vapply(
        psis[c("fallback_requested", "fallback_effective")],
        function(value) is.null(attributes(value)),
        logical(1)
      ))) {
    pv_abort("Reportable stack_psis source and fallback diagnostics are not canonical.")
  }
  expected_psis <- c(
    pv_stack_psis_diagnostics(canonical_pareto, threshold),
    list(
      pv_cols = pv_names,
      source = psis$source,
      pareto_k_source = psis$pareto_k_source,
      weight_method = psis$weight_method,
      producer = psis$producer,
      producer_version = psis$producer_version,
      normalization = psis$normalization,
      n_draws = psis$n_draws,
      ess_definition = psis$ess_definition,
      weight_ess_iid = psis$weight_ess_iid,
      weight_ess_fraction = psis$weight_ess_fraction,
      max_normalized_weight = psis$max_normalized_weight,
      weight_diagnostic_authority = psis$weight_diagnostic_authority,
      fallback_requested = psis$fallback_requested,
      fallback_effective = psis$fallback_effective
    )
  )
  if (!identical(psis, expected_psis)) {
    pv_abort("Reportable stack_psis diagnostics must equal the canonical Pareto-k projection.")
  }
  recomputed_pool <- rubin_pool_matrix(
    beta = weighted$beta,
    U = weighted$U,
    orientation = "rows_pv",
    conf_level = fit$control$conf_level,
    df_method = pooling$df_method,
    df_complete = pooling$df_complete
  )
  pooling_hash <- pv_stack_psis_pool_hash(
    list(
      pv_cols = pv_names,
      fe_names = fe_names_weighted,
      beta = weighted$beta,
      U = weighted$U
    ),
    recomputed_pool,
    expected_psis
  )
  expected_pooling <- list(
    beta = recomputed_pool$beta,
    U_bar = recomputed_pool$U_bar,
    B = recomputed_pool$B,
    T_MI = recomputed_pool$T_MI,
    lambda = recomputed_pool$lambda,
    df = recomputed_pool$df,
    df_classic = recomputed_pool$df_classic,
    df_method = recomputed_pool$df_method,
    df_complete = recomputed_pool$df_complete,
    pooling_hash = pooling_hash,
    pooling_source = "stack_psis_rubin_pooling"
  )
  if (!identical(pooling, expected_pooling)) {
    pv_abort("Reportable stack_psis pooling diagnostics must exactly reproduce weighted summaries.")
  }
  if (!is.data.frame(fit$estimates) || nrow(fit$estimates) == 0L) {
    pv_abort("Reportable stack_psis fit requires non-empty estimates.")
  }
  estimate_required <- c(
    "term", "estimate", "se", "std.error", "df", "df_method", "df_complete",
    "conf_level", "conf_low", "conf_high", "conf.low", "conf.high",
    "interval_role", "coverage_claim_allowed", "parameter_scope",
    "target_source", "target_hash", "pooling_source", "pooling_hash",
    "psis_status", "pareto_k_max", "psis_k_threshold", "psis_source",
    "pareto_k_source", "weight_method", "psis_producer",
    "psis_producer_version"
  )
  missing_estimate <- setdiff(estimate_required, names(fit$estimates))
  if (length(missing_estimate) > 0L) {
    pv_abort(sprintf("Reportable stack_psis fit estimates are missing required column(s): %s.", paste(missing_estimate, collapse = ", ")))
  }
  expected_estimates <- pv_stack_psis_estimates(
    recomputed_pool,
    expected_psis,
    pooling_hash
  )
  if (!identical(fit$estimates, expected_estimates)) {
    pv_abort(
      paste(
        "Reportable stack_psis estimates must exactly equal the canonical",
        "pooled projection, including interval, coverage, target, and PSIS metadata."
      )
    )
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
  if (!identical(as.numeric(fit$estimates$psis_k_threshold), rep(threshold, length(fe_names))) ||
      !identical(as.character(fit$estimates$psis_status), rep("ok", length(fe_names))) ||
      !identical(as.numeric(fit$estimates$pareto_k_max), rep(max(psis$pareto_k), length(fe_names))) ||
      !identical(as.character(fit$estimates$psis_source), rep(psis$source, length(fe_names))) ||
      !identical(as.character(fit$estimates$pareto_k_source), rep(psis$pareto_k_source, length(fe_names))) ||
      !identical(as.character(fit$estimates$weight_method), rep(psis$weight_method, length(fe_names))) ||
      !identical(as.character(fit$estimates$psis_producer), rep(psis$producer, length(fe_names))) ||
      !identical(as.character(fit$estimates$psis_producer_version), rep(psis$producer_version, length(fe_names)))) {
    pv_abort("Reportable stack_psis fit estimate PSIS metadata must match the complete diagnostics.")
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

pv_fit_required_fields <- function() {
  c(
    "method", "design", "target", "stack_fit", "ccc", "estimates",
    "draws", "diagnostics", "status", "control", "reason_codes",
    "schema_version", "provenance", "validation", "warnings"
  )
}

pv_fit_validation_schema <- function() {
  "pvstackr_fit_validation_v1"
}

pv_fit_validation_policy <- function() {
  "deep_semantics_plus_owned_payload_sha256_v1"
}

pv_fit_validation_canonicalizer <- function() {
  "r_xdr_v2_opaque_backend_boundary_v1"
}

pv_fit_validation_stamp_sentinel <- function() {
  paste0("sha256:", strrep("0", 64L))
}

pv_fit_has_opaque_backend <- function(fit) {
  stack_backend <- is.list(fit$stack_fit) &&
    "fit" %in% names(fit$stack_fit) &&
    !is.null(fit$stack_fit$fit)
  reference_backend <- is.list(fit$diagnostics) &&
    is.list(fit$diagnostics$reference) &&
    "backend_fits" %in% names(fit$diagnostics$reference) &&
    !is.null(fit$diagnostics$reference$backend_fits)
  isTRUE(stack_backend) || isTRUE(reference_backend)
}

pv_fit_validation_record <- function(fit, stamp = pv_fit_validation_stamp_sentinel()) {
  list(
    schema_version = pv_fit_validation_schema(),
    policy_id = pv_fit_validation_policy(),
    canonicalizer_id = pv_fit_validation_canonicalizer(),
    fast_path_eligible = !pv_fit_has_opaque_backend(fit),
    stamp = stamp
  )
}

pv_validate_fit_root_envelope <- function(fit) {
  pv_assert_named_list(fit, "fit")
  required <- pv_fit_required_fields()
  root_attributes <- attributes(fit)
  if (!identical(names(fit), required) ||
      !identical(names(root_attributes), c("names", "class")) ||
      !identical(root_attributes$names, required) ||
      !identical(root_attributes$class, c("pvstackr_fit", "list"))) {
    pv_abort("Fit object fields, order, class, and root attributes must be exact.")
  }
  invisible(fit)
}

pv_validate_fit_validation_record <- function(fit) {
  validation <- fit$validation
  required <- c(
    "schema_version", "policy_id", "canonicalizer_id",
    "fast_path_eligible", "stamp"
  )
  if (!is.list(validation) ||
      !identical(names(validation), required) ||
      !identical(attributes(validation), list(names = required)) ||
      !identical(validation$schema_version, pv_fit_validation_schema()) ||
      !identical(validation$policy_id, pv_fit_validation_policy()) ||
      !identical(validation$canonicalizer_id, pv_fit_validation_canonicalizer()) ||
      !is.logical(validation$fast_path_eligible) ||
      length(validation$fast_path_eligible) != 1L ||
      is.na(validation$fast_path_eligible) ||
      !identical(validation$fast_path_eligible, !pv_fit_has_opaque_backend(fit)) ||
      !is.character(validation$stamp) ||
      length(validation$stamp) != 1L ||
      is.na(validation$stamp) ||
      !grepl("^sha256:[0-9a-f]{64}$", validation$stamp)) {
    pv_abort("Fit validation record must use the exact current schema, policy, eligibility, and SHA-256 stamp envelope.")
  }
  invisible(validation)
}

pv_fit_validation_projection <- function(fit) {
  projected <- fit
  projected$validation$stamp <- pv_fit_validation_stamp_sentinel()

  # Backend fit objects are intentionally opaque retention exceptions. Their
  # presence and authorization remain covered, but their implementation-owned
  # contents are outside the portable validation stamp.
  if (is.list(projected$stack_fit) &&
      "fit" %in% names(projected$stack_fit) &&
      !is.null(projected$stack_fit$fit)) {
    projected$stack_fit$fit <- "pvstackr::opaque_backend_fit"
  }
  if (is.list(projected$diagnostics) &&
      is.list(projected$diagnostics$reference) &&
      "backend_fits" %in% names(projected$diagnostics$reference) &&
      !is.null(projected$diagnostics$reference$backend_fits)) {
    backend_fits <- projected$diagnostics$reference$backend_fits
    projected$diagnostics$reference$backend_fits <- list(
      marker = "pvstackr::opaque_backend_fits",
      length = length(backend_fits),
      names = names(backend_fits)
    )
  }

  # Formula environments affect both evaluation and retention safety. Only the
  # two canonical package environments are portable. The serialized formula
  # uses baseenv() for byte stability, while an adjacent marker preserves the
  # exact base-versus-stats identity. A private environment fails closed even
  # in the cheap tier instead of being normalized away.
  normalize_formulas <- function(x) {
    if (inherits(x, "formula")) {
      formula_environment <- environment(x)
      environment_id <- if (identical(formula_environment, baseenv())) {
        "base"
      } else if (identical(formula_environment, asNamespace("stats"))) {
        "namespace:stats"
      } else {
        pv_abort("Fit validation requires every retained formula environment to be canonical and portable.")
      }
      environment(x) <- baseenv()
      return(list(formula = x, environment_id = environment_id))
    }
    if (is.list(x)) {
      for (i in seq_along(x)) {
        x[i] <- list(normalize_formulas(x[[i]]))
      }
    }
    x
  }
  normalize_formulas(projected)
}

pv_fit_validation_digest <- function(fit) {
  projected <- pv_fit_validation_projection(fit)
  payload <- tryCatch(
    # Serialization v2 is used deliberately: unlike v3, it does not encode
    # ALTREP implementation state, so value-identical fits hash identically
    # before and after a save/read round trip.
    serialize(projected, NULL, ascii = FALSE, xdr = TRUE, version = 2L),
    error = function(error) {
      pv_abort(sprintf("Fit validation payload could not be serialized: %s", conditionMessage(error)))
    }
  )
  bytes <- c(
    charToRaw("pvstackr-fit-validation-v1"),
    as.raw(0L),
    payload
  )
  paste0("sha256:", digest::digest(bytes, algo = "sha256", serialize = FALSE))
}

pv_fit_issue_validation_stamp <- function(fit) {
  fit$validation <- pv_fit_validation_record(fit)
  fit$validation$stamp <- pv_fit_validation_digest(fit)
  fit
}

pv_validate_fit_validation_stamp <- function(fit) {
  observed <- pv_fit_validation_digest(fit)
  if (!identical(fit$validation$stamp, observed)) {
    pv_abort("Fit validation stamp does not match the current package-owned payload.")
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
    validation = NULL,
    warnings = warnings
  )
  class(fit) <- c("pvstackr_fit", "list")
  fit$validation <- pv_fit_validation_record(fit)
  pv_validate_pvstackr_fit_deep(fit)
  fit <- pv_fit_issue_validation_stamp(fit)
  pv_validate_fit_validation_stamp(fit)
  fit
}

pv_validate_pvstackr_fit_deep <- function(fit) {
  pv_validate_fit_root_envelope(fit)
  pv_validate_fit_validation_record(fit)
  method <- pv_validate_method(fit$method)
  pv_validate_control(fit$control)
  if (!identical(fit$control$method, method)) {
    pv_abort("Fit `control$method` must match `method`.")
  }
  pv_validate_fit_data_retention(fit)
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
  estimate_attributes <- attributes(fit$estimates)
  if (any(!names(estimate_attributes) %in% c("names", "row.names", "class")) ||
      !identical(class(fit$estimates), "data.frame") ||
      !identical(
        rownames(fit$estimates),
        as.character(seq_len(nrow(fit$estimates)))
      )) {
    pv_abort("Fit `estimates` must use canonical data-frame attributes and row names.")
  }
  if (nrow(fit$estimates) > 0L && "term" %in% names(fit$estimates)) {
    non_fe_terms <- fit$estimates$term[!grepl("^b_", fit$estimates$term)]
    if (length(non_fe_terms) > 0L) {
      pv_abort("Fit `estimates` must contain reportable fixed-effect terms only.")
    }
  }
  if (!is.null(fit$draws)) {
    draws <- ccc_as_draw_matrix(fit$draws)
    draw_attributes <- attributes(fit$draws)
    if (!identical(names(draw_attributes), c("dim", "dimnames")) ||
        !is.null(dimnames(fit$draws)[[1L]])) {
      pv_abort("Fit `draws` must use canonical matrix attributes without row labels.")
    }
    non_fe <- colnames(draws)[!grepl("^b_", colnames(draws))]
    if (length(non_fe) > 0L) {
      pv_abort("Fit `draws`, when retained, must contain reportable fixed-effect columns only.")
    }
  }
  if (!is.null(attributes(fit$warnings)) ||
      !is.null(attributes(fit$reason_codes))) {
    pv_abort("Fit warnings and reason codes must be bare character vectors.")
  }
  provenance_allowed <- c(
    "function_name", "package", "schema_version", "wrapper_function",
    "source", "target_hash", "target_manifest_hash", "target_content_hash",
    "binding_verification_policy", "independent_target_retained",
    "sampler_diagnostic_source", "ccc_target_hash",
    "stack_fit_long_data_hash", "stacked_source", "psis_source",
    "reportability_policy", "pooling_hash"
  )
  if (!is.list(fit$provenance) ||
      any(!names(fit$provenance) %in% provenance_allowed) ||
      anyDuplicated(names(fit$provenance)) ||
      !identical(attributes(fit$provenance), list(names = names(fit$provenance))) ||
      any(!vapply(
        fit$provenance,
        function(value) {
          length(value) == 1L && !is.object(value) && is.null(attributes(value))
        },
        logical(1)
      ))) {
    pv_abort("Fit provenance must use the canonical scalar-only retention envelope.")
  }
  pv_validate_stack_direct_reportable_fit(fit)
  pv_validate_per_pv_reportable_fit(fit)
  pv_validate_stack_psis_reportable_fit(fit)
  pv_validate_blocked_fit_redaction(fit)
  pv_validate_named_list_field(fit$diagnostics, "diagnostics")
  pv_validate_schema_version(fit$schema_version)
  pv_validate_named_list_field(fit$provenance, "provenance")
  pv_validate_fit_warnings(status, fit$warnings)
  invisible(fit)
}

validate_pvstackr_fit <- function(fit, tier = c("deep", "cheap")) {
  tier <- match.arg(tier)
  pv_validate_fit_root_envelope(fit)
  pv_validate_fit_validation_record(fit)

  if (identical(tier, "cheap") && isTRUE(fit$validation$fast_path_eligible)) {
    pv_validate_fit_validation_stamp(fit)
    return(invisible(fit))
  }

  # Opaque backend payloads cannot participate in a portable content digest,
  # so a requested cheap validation deliberately falls back to the deep tier.
  pv_validate_pvstackr_fit_deep(fit)
  pv_validate_fit_validation_stamp(fit)
  invisible(fit)
}
