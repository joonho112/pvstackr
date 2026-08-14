pv_binding_error_registry <- function() {
  c(
    PV_BIND_E001 = "MANIFEST_MISSING",
    PV_BIND_E002 = "SCHEMA_UNSUPPORTED",
    PV_BIND_E003 = "CANONICALIZER_UNSUPPORTED",
    PV_BIND_E004 = "HASH_ALGORITHM_UNSUPPORTED",
    PV_BIND_E005 = "MANIFEST_HASH_STALE",
    PV_BIND_E010 = "ROW_COUNT_MISMATCH",
    PV_BIND_E011 = "ROW_MEMBERSHIP_MISMATCH",
    PV_BIND_E012 = "ROW_ORDER_MISMATCH",
    PV_BIND_E020 = "PV_NAMES_OR_ORDER_MISMATCH",
    PV_BIND_E021 = "PV_VALUES_MISMATCH",
    PV_BIND_E030 = "PREDICTOR_SCHEMA_MISMATCH",
    PV_BIND_E031 = "PREDICTOR_VALUES_MISMATCH",
    PV_BIND_E040 = "FORMULA_MISMATCH",
    PV_BIND_E041 = "MODEL_MATRIX_SCHEMA_MISMATCH",
    PV_BIND_E042 = "MODEL_MATRIX_VALUES_MISMATCH",
    PV_BIND_E043 = "FACTOR_CONTRAST_MISMATCH",
    PV_BIND_E050 = "WEIGHT_NAMES_OR_ORDER_MISMATCH",
    PV_BIND_E051 = "WEIGHT_VALUES_MISMATCH",
    PV_BIND_E052 = "WEIGHT_POLICY_MISMATCH",
    PV_BIND_E053 = "REPLICATE_DESIGN_MISMATCH",
    PV_BIND_E060 = "FAMILY_MISMATCH",
    PV_BIND_E061 = "LINK_MISMATCH",
    PV_BIND_E070 = "ESTIMAND_MISMATCH",
    PV_BIND_E071 = "COEFFICIENT_ORDER_MISMATCH",
    PV_BIND_E072 = "ESTIMAND_CONTRAST_MISMATCH",
    PV_BIND_E080 = "LEGACY_NOT_REVALIDATED",
    PV_BIND_E081 = "BINDING_RECOMPUTE_FAILED",
    PV_BIND_E090 = "TARGET_CONTENT_HASH_STALE"
  )
}

pv_binding_hash_prefix_hex_length <- function() {
  12L
}

pv_binding_sanitize_hash_metadata <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.character(x) && length(x) == 1L && !is.na(x)) {
    if (grepl("^sha256:[0-9a-f]{64}$", x)) {
      return(substr(x, 1L, 7L + pv_binding_hash_prefix_hex_length()))
    }
    if (grepl(sprintf("^sha256:[0-9a-f]{%d}$", pv_binding_hash_prefix_hex_length()), x)) {
      return(x)
    }
    if (identical(x, "malformed")) {
      return(x)
    }
  }
  "malformed"
}

pv_binding_abort <- function(
  code,
  detail,
  component = NULL,
  expected_hash = NULL,
  observed_hash = NULL,
  all_codes = NULL
) {
  registry <- pv_binding_error_registry()
  if (!is.character(code) || length(code) != 1L || is.na(code) || !code %in% names(registry)) {
    stop("pvstackr: internal binding error used an unknown error code.", call. = FALSE)
  }
  if (!is.character(detail) || length(detail) != 1L || is.na(detail) || !nzchar(detail)) {
    stop("pvstackr: internal binding error requires a non-empty detail.", call. = FALSE)
  }
  scalar_or_null <- function(x) {
    is.null(x) || (is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x))
  }
  expected_hash <- pv_binding_sanitize_hash_metadata(expected_hash)
  observed_hash <- pv_binding_sanitize_hash_metadata(observed_hash)
  if (!scalar_or_null(component) || !scalar_or_null(expected_hash) || !scalar_or_null(observed_hash)) {
    stop("pvstackr: internal binding error metadata must be NULL or a non-empty string.", call. = FALSE)
  }
  if (!is.null(all_codes)) {
    registry_order <- names(registry)
    if (!is.character(all_codes) || length(all_codes) < 1L || anyNA(all_codes) ||
        anyDuplicated(all_codes) || any(!all_codes %in% registry_order) ||
        !identical(all_codes, registry_order[registry_order %in% all_codes]) ||
        !identical(code, all_codes[[1L]])) {
      stop(
        "pvstackr: internal binding error `all_codes` must be a registry-ordered vector beginning with `code`.",
        call. = FALSE
      )
    }
  }

  condition_name <- unname(registry[[code]])
  condition_class <- paste0("pvstackr_binding_", tolower(condition_name))
  message <- sprintf("pvstackr [%s %s]: %s", code, condition_name, detail)
  condition <- structure(
    list(
      message = message,
      call = NULL,
      code = code,
      condition_name = condition_name,
      component = component,
      expected_hash = expected_hash,
      observed_hash = observed_hash,
      all_codes = all_codes
    ),
    class = c(
      condition_class,
      "pvstackr_binding_error",
      "pvstackr_error",
      "error",
      "condition"
    )
  )
  stop(condition)
}

pv_binding_hash_domains <- function() {
  c(
    "row", "pv", "predictor", "formula", "model_matrix", "weight",
    "factor_contrast", "estimand_contrast", "family_link", "estimand",
    "target_content", "manifest"
  )
}

pv_binding_formula_rhs_info <- function(formula, data = NULL) {
  if (!inherits(formula, "formula") || length(formula) != 3L ||
      !identical(formula[[2L]], as.name("OUTCOME"))) {
    pv_binding_abort(
      "PV_BIND_E081",
      "A binding formula must be two-sided and use OUTCOME as its left-hand placeholder.",
      "formula"
    )
  }
  if (!is.null(data)) {
    if (!is.data.frame(data) || is.null(names(data)) || anyNA(names(data)) ||
        any(!nzchar(names(data))) || anyDuplicated(names(data))) {
      pv_binding_abort(
        "PV_BIND_E081",
        "Formula validation requires a data frame with unique non-empty column names.",
        "formula"
      )
    }
  }

  rhs <- formula[[3L]]
  data_variables <- unique(all.vars(rhs))
  if (!is.null(data)) {
    missing_variables <- setdiff(data_variables, names(data))
    if (length(missing_variables) > 0L) {
      pv_binding_abort(
        "PV_BIND_E081",
        sprintf(
          "Every RHS variable must be an explicit data column; missing: %s.",
          paste(missing_variables, collapse = ", ")
        ),
        "formula"
      )
    }
  }

  list(
    rhs = rhs,
    data_variables = data_variables
  )
}

pv_binding_language_ast <- function(expr) {
  if (is.null(expr)) {
    return(list(node_type = "null"))
  }
  if (is.symbol(expr)) {
    return(list(
      node_type = "symbol",
      name = as.character(expr)
    ))
  }
  if (is.call(expr)) {
    call_parts <- as.list(expr)
    arguments <- call_parts[-1L]
    argument_names <- names(call_parts)[-1L]
    if (is.null(argument_names)) {
      argument_names <- rep("", length(arguments))
    } else {
      argument_names[is.na(argument_names)] <- ""
    }
    return(list(
      node_type = "call",
      head = pv_binding_language_ast(call_parts[[1L]]),
      argument_names = unname(argument_names),
      arguments = unname(lapply(arguments, pv_binding_language_ast))
    ))
  }
  if (is.atomic(expr) && length(expr) == 1L &&
      typeof(expr) %in% c("logical", "integer", "double", "character")) {
    attributes(expr) <- NULL
    return(list(
      node_type = "constant",
      constant_type = typeof(expr),
      value = expr
    ))
  }
  pv_binding_abort(
    "PV_BIND_E081",
    "The formula contains a language node unsupported by the binding AST.",
    "formula"
  )
}

pv_binding_formula_ast <- function(formula, data = NULL) {
  info <- pv_binding_formula_rhs_info(formula, data = data)
  list(
    node_type = "formula",
    outcome_placeholder = "OUTCOME",
    data_variables = info$data_variables,
    rhs = pv_binding_language_ast(info$rhs)
  )
}

pv_binding_predictor_type_id <- function(x) {
  if (is.ordered(x)) {
    return("ordered_factor")
  }
  if (is.factor(x)) {
    return("factor")
  }
  if (inherits(x, "POSIXct")) {
    return("POSIXct")
  }
  if (inherits(x, "Date")) {
    return("Date")
  }
  if (is.character(x) && is.null(attr(x, "class", exact = TRUE))) {
    return("character")
  }
  if (is.logical(x) && is.null(attr(x, "class", exact = TRUE))) {
    return("logical")
  }
  if (is.integer(x) && is.null(attr(x, "class", exact = TRUE))) {
    return("integer")
  }
  if (is.double(x) && is.null(attr(x, "class", exact = TRUE))) {
    return("double")
  }
  pv_binding_abort(
    "PV_BIND_E081",
    "A formula predictor has an unsupported type or class.",
    "predictor"
  )
}

pv_binding_validate_predictor <- function(
  x,
  allow_character = FALSE,
  component = "predictor"
) {
  type_id <- pv_binding_predictor_type_id(x)
  if (type_id == "character" && !allow_character) {
    pv_binding_abort(
      "PV_BIND_E081",
      "Character formula predictors are forbidden; convert them to an explicit factor with declared level order.",
      component
    )
  }
  if (anyNA(x)) {
    pv_binding_abort(
      "PV_BIND_E081",
      "Formula predictors must not contain missing values.",
      component
    )
  }
  if (type_id %in% c("double", "Date", "POSIXct") && any(!is.finite(as.double(x)))) {
    pv_binding_abort(
      "PV_BIND_E081",
      "Numeric, Date, and POSIXct predictors must contain only finite values.",
      component
    )
  }
  if (type_id == "character") {
    values <- enc2utf8(x)
    if (any(!validUTF8(values))) {
      pv_binding_abort(
        "PV_BIND_E081",
        "Character predictors must have valid UTF-8 encoding.",
        component
      )
    }
  }
  type_id
}

pv_binding_ordered_labels_hash <- function(labels, domain, role) {
  pv_binding_hash_payload(
    list(role = role, count = length(labels), ordered_labels = as.character(labels)),
    domain
  )
}

# Survey files carry cosmetic metadata from their SAS, SPSS, or Stata origin:
# nearly every PISA column arrives with a `label`, and many with a `format.*`
# string. Those attributes are not part of the analysed content, and the
# projections record the structural properties they do depend on -- type, factor
# levels, orderedness, timezone -- as explicit manifest fields. Reducing a plain
# atomic column to its values and names keeps the binding hash stable across
# files that differ only in that metadata. Anything carrying a class is returned
# untouched, so factors, Dates, and POSIXct still reach the canonicaliser's own
# typed branches and any unaccounted attribute there is still refused.
pv_binding_column_values <- function(x) {
  if (!is.atomic(x) || !is.null(attr(x, "class", exact = TRUE))) {
    return(x)
  }
  value_names <- names(x)
  attributes(x) <- NULL
  if (!is.null(value_names)) {
    names(x) <- value_names
  }
  x
}

pv_binding_predictor_projection <- function(data, formula) {
  info <- pv_binding_formula_rhs_info(formula, data = data)
  predictor_names <- info$data_variables
  columns <- lapply(predictor_names, function(name) {
    values <- pv_binding_column_values(data[[name]])
    type_id <- pv_binding_validate_predictor(values)
    factor_levels <- if (is.factor(values)) levels(values) else character()
    timezone <- if (inherits(values, "POSIXct")) attr(values, "tzone", exact = TRUE) else NULL
    if (is.null(timezone)) {
      timezone <- character()
    }
    list(
      name = name,
      type_id = type_id,
      ordered = is.ordered(values),
      level_count = length(factor_levels),
      ordered_level_labels_hash = pv_binding_ordered_labels_hash(
        factor_levels,
        "predictor",
        "raw_factor_levels"
      ),
      timezone_hash = pv_binding_ordered_labels_hash(timezone, "predictor", "timezone"),
      values_hash = pv_binding_hash_payload(
        list(role = "raw_predictor_values", name = name, values = values),
        "predictor"
      )
    )
  })
  names(columns) <- NULL
  schema_payload <- lapply(columns, function(column) {
    column[setdiff(names(column), "values_hash")]
  })
  ordered_value_payload <- lapply(predictor_names, function(name) {
    list(name = name, values = pv_binding_column_values(data[[name]]))
  })

  list(
    source_rule = "rhs_ast_first_occurrence_v1",
    row_count = nrow(data),
    column_count = length(predictor_names),
    names = predictor_names,
    columns = columns,
    schema_hash = pv_binding_hash_payload(schema_payload, "predictor"),
    ordered_values_hash = pv_binding_hash_payload(ordered_value_payload, "predictor")
  )
}

pv_binding_safe_model_eval <- function(expr, stage) {
  tryCatch(
    withCallingHandlers(
      expr,
      warning = function(warning) {
        pv_binding_abort(
          "PV_BIND_E081",
          sprintf("Response-free %s evaluation emitted a warning.", stage),
          "model_matrix"
        )
      }
    ),
    pvstackr_binding_error = function(error) stop(error),
    error = function(error) {
      pv_binding_abort(
        "PV_BIND_E081",
        sprintf("Response-free %s construction failed.", stage),
        "model_matrix"
      )
    }
  )
}

pv_binding_numeric_matrix_hash <- function(x, domain, role) {
  if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x))) {
    pv_binding_abort(
      "PV_BIND_E081",
      "A resolved model or contrast matrix must be finite and numeric.",
      "model_matrix"
    )
  }
  pv_binding_hash_payload(
    list(
      role = role,
      dimensions = as.integer(dim(x)),
      values_column_major = as.double(x)
    ),
    domain
  )
}

pv_binding_keyed_matrix_values_hash <- function(x, row_digests, role) {
  if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x)) ||
      !is.character(row_digests) || length(row_digests) != nrow(x) ||
      anyNA(row_digests) || any(!grepl("^sha256:[0-9a-f]{64}$", row_digests))) {
    pv_binding_abort(
      "PV_BIND_E081",
      "Keyed model-matrix values require finite matrix rows and canonical row digests.",
      "model_matrix"
    )
  }
  order_index <- order(row_digests, method = "radix")
  keyed_rows <- lapply(order_index, function(index) {
    list(
      row_digest = row_digests[[index]],
      values = as.double(x[index, , drop = TRUE])
    )
  })
  pv_binding_hash_payload(
    list(role = role, column_count = ncol(x), keyed_rows = keyed_rows),
    "model_matrix"
  )
}

pv_binding_keyed_offset_hash <- function(offset, row_digests) {
  if (is.null(offset)) {
    return(pv_binding_hash_payload(
      list(role = "offset_membership_values", has_offset = FALSE),
      "model_matrix"
    ))
  }
  if (!is.numeric(offset) || any(!is.finite(offset)) || length(offset) != length(row_digests)) {
    pv_binding_abort("PV_BIND_E081", "Keyed offset values must be finite and row aligned.", "model_matrix")
  }
  order_index <- order(row_digests, method = "radix")
  keyed_rows <- lapply(order_index, function(index) {
    list(row_digest = row_digests[[index]], value = as.double(offset[[index]]))
  })
  pv_binding_hash_payload(
    list(role = "offset_membership_values", has_offset = TRUE, keyed_rows = keyed_rows),
    "model_matrix"
  )
}

pv_binding_xlevel_projection <- function(frame) {
  factor_names <- names(frame)[vapply(frame, is.factor, logical(1))]
  specs <- lapply(factor_names, function(name) {
    values <- frame[[name]]
    list(
      variable = name,
      ordered = is.ordered(values),
      level_count = length(levels(values)),
      ordered_level_labels_hash = pv_binding_ordered_labels_hash(
        levels(values),
        "model_matrix",
        "resolved_xlevels"
      )
    )
  })
  names(specs) <- NULL
  specs
}

pv_binding_resolved_contrast_matrices <- function(frame, model_matrix) {
  contrast_map <- attr(model_matrix, "contrasts", exact = TRUE)
  if (is.null(contrast_map)) {
    return(list())
  }
  contrast_names <- names(contrast_map)
  resolved <- lapply(seq_along(contrast_map), function(index) {
    variable <- contrast_names[[index]]
    frame_value <- frame[[variable]]
    if (is.factor(frame_value) || is.logical(frame_value)) {
      pv_binding_safe_model_eval(
        stats::contrasts(frame_value, contrasts = TRUE),
        "factor contrast"
      )
    } else {
      matrix(numeric(), nrow = 0L, ncol = 0L)
    }
  })
  names(resolved) <- NULL
  resolved
}

pv_binding_contrast_projection <- function(
  frame,
  model_matrix,
  resolved_contrasts = NULL
) {
  contrast_map <- attr(model_matrix, "contrasts", exact = TRUE)
  if (is.null(contrast_map)) {
    if (!is.null(resolved_contrasts) && length(resolved_contrasts) != 0L) {
      pv_binding_abort(
        "PV_BIND_E081",
        "Resolved contrast matrices do not align with the model matrix.",
        "model_matrix"
      )
    }
    return(list())
  }
  if (is.null(resolved_contrasts)) {
    resolved_contrasts <- pv_binding_resolved_contrast_matrices(frame, model_matrix)
  }
  if (!is.list(resolved_contrasts) ||
      length(resolved_contrasts) != length(contrast_map) ||
      !is.null(names(resolved_contrasts))) {
    pv_binding_abort(
      "PV_BIND_E081",
      "Resolved contrast matrices must be an unnamed list aligned with the model matrix.",
      "model_matrix"
    )
  }
  contrast_names <- names(contrast_map)
  specs <- lapply(seq_along(contrast_map), function(index) {
    variable <- contrast_names[[index]]
    generator <- as.character(contrast_map[[index]])
    frame_value <- frame[[variable]]
    resolved <- resolved_contrasts[[index]]
    if (!(is.factor(frame_value) || is.logical(frame_value)) || !is.matrix(resolved) ||
        !is.numeric(resolved) || any(!is.finite(resolved))) {
      pv_binding_abort(
        "PV_BIND_E081",
        "Resolved factor contrasts must be finite numeric matrices aligned with factor variables.",
        "model_matrix"
      )
    }
    resolved_dimnames <- dimnames(resolved)
    row_labels <- if (is.null(resolved_dimnames) || is.null(resolved_dimnames[[1L]])) {
      character()
    } else {
      resolved_dimnames[[1L]]
    }
    col_labels <- if (is.null(resolved_dimnames) || is.null(resolved_dimnames[[2L]])) {
      character()
    } else {
      resolved_dimnames[[2L]]
    }
    list(
      variable = variable,
      generator_label_hash = pv_binding_ordered_labels_hash(
        generator,
        "factor_contrast",
        "generator"
      ),
      dimensions = as.integer(dim(resolved)),
      row_label_count = length(row_labels),
      row_label_hash = pv_binding_ordered_labels_hash(
        row_labels,
        "factor_contrast",
        "row_labels"
      ),
      column_label_count = length(col_labels),
      column_label_hash = pv_binding_ordered_labels_hash(
        col_labels,
        "factor_contrast",
        "column_labels"
      ),
      values_hash = pv_binding_numeric_matrix_hash(
        resolved,
        "factor_contrast",
        "resolved_contrast_values"
      )
    )
  })
  names(specs) <- NULL
  specs
}

pv_binding_model_bundle_fields <- function() {
  c(
    "bundle_schema_version", "formula_ast", "formula_environment",
    "predictor_projection",
    "row_count", "row_order_hash", "term_labels", "intercept", "frame",
    "model_matrix", "offset", "xlevels", "contrast_option_labels",
    "resolved_contrasts", "contrasts", "frame_values_hash",
    "model_matrix_values_hash",
    "offset_values_hash", "bundle_hash"
  )
}

pv_binding_model_bundle_frame_hash <- function(frame) {
  if (!is.data.frame(frame) || is.null(names(frame)) || anyNA(names(frame)) ||
      any(!nzchar(names(frame))) || anyDuplicated(names(frame))) {
    pv_binding_abort(
      "PV_BIND_E081",
      "A resolved model bundle frame must be a uniquely named data frame.",
      "model_matrix"
    )
  }
  columns <- lapply(seq_along(frame), function(index) {
    value <- frame[[index]]
    if (is.factor(value)) {
      return(list(
        kind = "factor",
        storage_type = typeof(value),
        class_labels = class(value),
        ordered = is.ordered(value),
        levels = levels(value),
        values = unname(as.integer(value))
      ))
    }
    if (is.matrix(value)) {
      flattened <- as.vector(value)
      attributes(flattened) <- NULL
      if (!typeof(flattened) %in% c("logical", "integer", "double", "character", "raw") ||
          (is.numeric(flattened) && any(!is.finite(flattened)))) {
        pv_binding_abort(
          "PV_BIND_E081",
          "Resolved model-frame matrix columns use an unsupported value type.",
          "model_matrix"
        )
      }
      return(list(
        kind = "matrix",
        storage_type = typeof(value),
        class_labels = class(value),
        dimensions = as.integer(dim(value)),
        dimnames = dimnames(value),
        values = flattened
      ))
    }
    if (!is.atomic(value) ||
        !typeof(value) %in% c("logical", "integer", "double", "character", "raw")) {
      pv_binding_abort(
        "PV_BIND_E081",
        "Resolved model-frame columns use an unsupported value type.",
        "model_matrix"
      )
    }
    flattened <- value
    attributes(flattened) <- NULL
    if (is.numeric(flattened) && any(!is.finite(flattened))) {
      pv_binding_abort(
        "PV_BIND_E081",
        "Resolved numeric model-frame columns must be finite.",
        "model_matrix"
      )
    }
    list(
      kind = "atomic",
      storage_type = typeof(value),
      class_labels = class(value),
      values = flattened
    )
  })
  names(columns) <- NULL
  pv_binding_hash_payload(
    list(
      role = "resolved_model_frame_values",
      row_count = as.integer(nrow(frame)),
      column_names = names(frame),
      columns = columns
    ),
    "model_matrix"
  )
}

pv_binding_model_bundle_row_order_hash <- function(predictor_projection, row_count) {
  pv_binding_hash_payload(
    list(
      role = "resolved_model_bundle_raw_row_order",
      row_count = as.integer(row_count),
      predictor_names = predictor_projection$names,
      predictor_schema_hash = predictor_projection$schema_hash,
      predictor_values_hash = predictor_projection$ordered_values_hash
    ),
    "model_matrix"
  )
}

pv_binding_model_bundle_hash_payload <- function(bundle) {
  contrast_map <- attr(bundle$model_matrix, "contrasts", exact = TRUE)
  if (is.null(contrast_map)) {
    contrast_map <- character()
  }
  list(
    bundle_schema_version = bundle$bundle_schema_version,
    formula_ast_hash = pv_binding_hash_payload(bundle$formula_ast, "formula"),
    predictor_projection = bundle$predictor_projection,
    row_count = bundle$row_count,
    row_order_hash = bundle$row_order_hash,
    term_labels = bundle$term_labels,
    intercept = bundle$intercept,
    model_matrix_schema = list(
      dimensions = as.integer(dim(bundle$model_matrix)),
      column_names = colnames(bundle$model_matrix),
      assign = as.integer(attr(bundle$model_matrix, "assign", exact = TRUE)),
      contrasts = contrast_map
    ),
    model_matrix_values_hash = bundle$model_matrix_values_hash,
    offset_values_hash = bundle$offset_values_hash,
    xlevels = bundle$xlevels,
    contrast_option_labels = bundle$contrast_option_labels,
    contrasts = bundle$contrasts,
    frame_values_hash = bundle$frame_values_hash
  )
}

pv_binding_resolve_model_bundle <- function(data, formula) {
  predictor_projection <- pv_binding_predictor_projection(data, formula)
  info <- pv_binding_formula_rhs_info(formula, data = data)
  rhs_formula <- stats::as.formula(call("~", info$rhs), env = environment(formula))
  terms_object <- pv_binding_safe_model_eval(
    stats::terms(rhs_formula, data = data),
    "terms"
  )
  frame <- pv_binding_safe_model_eval(
    stats::model.frame(
      terms_object,
      data = data,
      na.action = stats::na.fail,
      drop.unused.levels = FALSE
    ),
    "model frame"
  )
  model_matrix <- pv_binding_safe_model_eval(
    stats::model.matrix(terms_object, data = frame),
    "model matrix"
  )
  if (!is.numeric(model_matrix) || any(!is.finite(model_matrix))) {
    pv_binding_abort(
      "PV_BIND_E081",
      "The resolved model matrix must contain only finite numeric values.",
      "model_matrix"
    )
  }
  offset <- stats::model.offset(frame)
  if (!is.null(offset) && (!is.numeric(offset) || any(!is.finite(offset)))) {
    pv_binding_abort(
      "PV_BIND_E081",
      "The resolved offset must contain only finite numeric values.",
      "model_matrix"
    )
  }
  if (!is.null(offset)) {
    offset <- unname(as.double(offset))
  }
  term_labels <- unname(attr(terms_object, "term.labels", exact = TRUE))
  xlevels <- pv_binding_xlevel_projection(frame)
  resolved_contrasts <- pv_binding_resolved_contrast_matrices(frame, model_matrix)
  contrasts <- pv_binding_contrast_projection(
    frame,
    model_matrix,
    resolved_contrasts = resolved_contrasts
  )
  option_labels <- getOption("contrasts")
  if (is.null(option_labels)) {
    option_labels <- character()
  }
  option_labels <- unname(as.character(option_labels))
  bundle <- list(
    bundle_schema_version = "pvstackr_resolved_model_bundle_v1",
    formula_ast = pv_binding_formula_ast(formula, data = data),
    formula_environment = environment(formula),
    predictor_projection = predictor_projection,
    row_count = as.integer(nrow(data)),
    row_order_hash = pv_binding_model_bundle_row_order_hash(
      predictor_projection,
      nrow(data)
    ),
    term_labels = term_labels,
    intercept = as.integer(attr(terms_object, "intercept", exact = TRUE)),
    frame = frame,
    model_matrix = model_matrix,
    offset = offset,
    xlevels = xlevels,
    contrast_option_labels = option_labels,
    resolved_contrasts = resolved_contrasts,
    contrasts = contrasts,
    frame_values_hash = pv_binding_model_bundle_frame_hash(frame),
    model_matrix_values_hash = pv_binding_numeric_matrix_hash(
      model_matrix,
      "model_matrix",
      "resolved_model_matrix_values"
    ),
    offset_values_hash = pv_binding_hash_payload(
      list(role = "offset_values", values = offset),
      "model_matrix"
    )
  )
  bundle$bundle_hash <- pv_binding_hash_payload(
    pv_binding_model_bundle_hash_payload(bundle),
    "model_matrix"
  )
  pv_binding_validate_model_bundle(bundle, data, formula)
  bundle
}

pv_binding_validate_model_bundle <- function(bundle, data, formula) {
  pv_binding_assert_exact_fields(
    bundle,
    pv_binding_model_bundle_fields(),
    "resolved model bundle",
    "PV_BIND_E081"
  )
  if (!identical(bundle$bundle_schema_version, "pvstackr_resolved_model_bundle_v1")) {
    pv_binding_abort("PV_BIND_E081", "The resolved model bundle schema is unsupported.", "model_matrix")
  }
  if (!is.environment(bundle$formula_environment) ||
      !identical(bundle$formula_environment, environment(formula))) {
    pv_binding_abort(
      "PV_BIND_E081",
      "The resolved model bundle cannot be reused across formula evaluation environments.",
      "model_matrix"
    )
  }
  current_ast <- pv_binding_formula_ast(formula, data = data)
  current_predictor <- pv_binding_predictor_projection(data, formula)
  current_row_count <- as.integer(nrow(data))
  current_row_hash <- pv_binding_model_bundle_row_order_hash(
    current_predictor,
    current_row_count
  )
  if (!identical(bundle$formula_ast, current_ast) ||
      !identical(bundle$predictor_projection, current_predictor) ||
      !identical(bundle$row_count, current_row_count) ||
      !identical(bundle$row_order_hash, current_row_hash)) {
    pv_binding_abort(
      "PV_BIND_E081",
      "The resolved model bundle does not match the supplied formula and raw predictor rows.",
      "model_matrix"
    )
  }
  frame <- bundle$frame
  model_matrix <- bundle$model_matrix
  offset <- bundle$offset
  frame_attribute_names <- names(attributes(frame))
  required_frame_attributes <- c("names", "terms", "row.names", "class")
  model_matrix_attribute_names <- names(attributes(model_matrix))
  required_matrix_attributes <- c("dim", "dimnames", "assign")
  allowed_matrix_attributes <- c(required_matrix_attributes, "contrasts")
  if (!is.data.frame(frame) || nrow(frame) != current_row_count ||
      !setequal(frame_attribute_names, required_frame_attributes) ||
      !is.matrix(model_matrix) || !is.numeric(model_matrix) ||
      any(!is.finite(model_matrix)) || nrow(model_matrix) != current_row_count ||
      ncol(model_matrix) < 1L || is.null(colnames(model_matrix)) ||
      anyNA(colnames(model_matrix)) || any(!nzchar(colnames(model_matrix))) ||
      anyDuplicated(colnames(model_matrix)) ||
      !all(required_matrix_attributes %in% model_matrix_attribute_names) ||
      any(!model_matrix_attribute_names %in% allowed_matrix_attributes) ||
      anyDuplicated(model_matrix_attribute_names)) {
    pv_binding_abort(
      "PV_BIND_E081",
      "The resolved model bundle frame and model matrix are malformed or row-misaligned.",
      "model_matrix"
    )
  }
  assign <- attr(model_matrix, "assign", exact = TRUE)
  if (!is.integer(assign) || length(assign) != ncol(model_matrix) ||
      anyNA(assign) || any(assign < 0L)) {
    pv_binding_abort(
      "PV_BIND_E081",
      "The resolved model bundle model-matrix assignment is malformed.",
      "model_matrix"
    )
  }
  if (!(is.null(offset) ||
      (is.numeric(offset) && is.null(attributes(offset)) &&
        length(offset) == current_row_count && all(is.finite(offset))))) {
    pv_binding_abort(
      "PV_BIND_E081",
      "The resolved model bundle offset is malformed or row-misaligned.",
      "model_matrix"
    )
  }
  frame_offset <- pv_binding_safe_model_eval(
    stats::model.offset(frame),
    "bundle offset validation"
  )
  if (!is.null(frame_offset)) {
    frame_offset <- unname(as.double(frame_offset))
  }
  if (!identical(frame_offset, offset)) {
    pv_binding_abort(
      "PV_BIND_E081",
      "The resolved model bundle frame and stored offset disagree.",
      "model_matrix"
    )
  }

  info <- pv_binding_formula_rhs_info(formula, data = data)
  rhs_formula <- stats::as.formula(call("~", info$rhs), env = environment(formula))
  terms_object <- pv_binding_safe_model_eval(
    stats::terms(rhs_formula, data = data),
    "bundle terms validation"
  )
  expected_term_labels <- unname(attr(terms_object, "term.labels", exact = TRUE))
  expected_intercept <- as.integer(attr(terms_object, "intercept", exact = TRUE))
  if (!is.character(bundle$term_labels) || !is.null(attributes(bundle$term_labels)) ||
      !identical(bundle$term_labels, expected_term_labels) ||
      !identical(bundle$intercept, expected_intercept)) {
    pv_binding_abort(
      "PV_BIND_E081",
      "The resolved model bundle terms do not match the supplied formula.",
      "model_matrix"
    )
  }
  expected_xlevels <- pv_binding_xlevel_projection(frame)
  if (!is.list(bundle$resolved_contrasts) ||
      !is.null(attributes(bundle$resolved_contrasts)) ||
      any(!vapply(bundle$resolved_contrasts, function(value) {
        is.matrix(value) && setequal(
          names(attributes(value)),
          c("dim", "dimnames")
        )
      }, logical(1)))) {
    pv_binding_abort(
      "PV_BIND_E081",
      "The resolved contrast matrices carry a malformed transient representation.",
      "model_matrix"
    )
  }
  expected_contrasts <- pv_binding_contrast_projection(
    frame,
    model_matrix,
    resolved_contrasts = bundle$resolved_contrasts
  )
  if (!identical(bundle$xlevels, expected_xlevels) ||
      !identical(bundle$contrasts, expected_contrasts) ||
      !is.character(bundle$contrast_option_labels) ||
      !is.null(attributes(bundle$contrast_option_labels)) ||
      length(bundle$contrast_option_labels) != 2L) {
    pv_binding_abort(
      "PV_BIND_E081",
      "The resolved factor and contrast state is not canonical for the bundle.",
      "model_matrix"
    )
  }
  expected_matrix_hash <- pv_binding_numeric_matrix_hash(
    model_matrix,
    "model_matrix",
    "resolved_model_matrix_values"
  )
  expected_offset_hash <- pv_binding_hash_payload(
    list(role = "offset_values", values = offset),
    "model_matrix"
  )
  expected_bundle_hash <- pv_binding_hash_payload(
    pv_binding_model_bundle_hash_payload(bundle),
    "model_matrix"
  )
  expected_frame_hash <- pv_binding_model_bundle_frame_hash(frame)
  if (!identical(bundle$frame_values_hash, expected_frame_hash) ||
      !identical(bundle$model_matrix_values_hash, expected_matrix_hash) ||
      !identical(bundle$offset_values_hash, expected_offset_hash) ||
      !identical(bundle$bundle_hash, expected_bundle_hash)) {
    pv_binding_abort(
      "PV_BIND_E081",
      "The resolved model bundle matrix, offset, or integrity hash is stale.",
      "model_matrix"
    )
  }
  invisible(bundle)
}

pv_binding_model_matrix_projection <- function(
  data,
  formula,
  row_digests = NULL,
  model_bundle = NULL
) {
  if (is.null(model_bundle)) {
    model_bundle <- pv_binding_resolve_model_bundle(data, formula)
  }
  pv_binding_validate_model_bundle(model_bundle, data, formula)
  predictor_projection <- model_bundle$predictor_projection
  model_matrix <- model_bundle$model_matrix
  offset <- model_bundle$offset
  if (is.null(row_digests)) {
    predictor_names <- predictor_projection$names
    if (length(predictor_names) > 0L) {
      row_digests <- pv_binding_row_identity_digests(
        data,
        bound_cols = predictor_names
      )$digests
    } else {
      row_digests <- vapply(seq_len(nrow(data)), function(index) {
        pv_binding_hash_payload(
          list(role = "model_row_position", index = as.integer(index)),
          "row"
        )
      }, character(1))
    }
  }
  term_labels <- model_bundle$term_labels
  xlevels <- model_bundle$xlevels
  contrasts <- model_bundle$contrasts
  reportable_fe_names <- paste0("b_", colnames(model_matrix))
  reportable_fe_names[reportable_fe_names == "b_(Intercept)"] <- "b_Intercept"

  list(
    builder_id = "stats_model_frame_matrix_v1",
    na_action = "na.fail",
    row_count = nrow(model_matrix),
    column_count = ncol(model_matrix),
    fe_colname_count = ncol(model_matrix),
    ordered_fe_colnames_hash = pv_binding_ordered_labels_hash(
      colnames(model_matrix),
      "model_matrix",
      "fe_colnames"
    ),
    ordered_reportable_fe_names_hash = pv_binding_ordered_labels_hash(
      reportable_fe_names,
      "model_matrix",
      "reportable_fe_names"
    ),
    assign = as.integer(attr(model_matrix, "assign", exact = TRUE)),
    intercept = model_bundle$intercept,
    term_count = length(term_labels),
    ordered_term_labels_hash = pv_binding_ordered_labels_hash(
      term_labels,
      "model_matrix",
      "term_labels"
    ),
    has_offset = !is.null(offset),
    offset_values_hash = model_bundle$offset_values_hash,
    xlevel_count = length(xlevels),
    xlevels = xlevels,
    contrast_count = length(contrasts),
    contrasts = contrasts,
    values_hash = model_bundle$model_matrix_values_hash,
    membership_values_hash = pv_binding_keyed_matrix_values_hash(
      model_matrix,
      row_digests,
      "resolved_model_matrix_membership_values"
    ),
    offset_membership_hash = pv_binding_keyed_offset_hash(offset, row_digests),
    predictor_schema_hash = predictor_projection$schema_hash,
    predictor_values_hash = predictor_projection$ordered_values_hash
  )
}

pv_binding_validate_column_spec <- function(data, columns, label, min_length = 1L) {
  if (!is.data.frame(data) || is.null(names(data)) || anyNA(names(data)) ||
      any(!nzchar(names(data))) || anyDuplicated(names(data))) {
    pv_binding_abort(
      "PV_BIND_E081",
      "Binding inputs require a data frame with unique non-empty column names.",
      label
    )
  }
  if (!is.character(columns) || length(columns) < min_length || anyNA(columns) ||
      any(!nzchar(columns)) || anyDuplicated(columns)) {
    pv_binding_abort(
      "PV_BIND_E081",
      sprintf("%s columns must be unique non-empty names.", label),
      label
    )
  }
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0L) {
    pv_binding_abort(
      "PV_BIND_E081",
      sprintf("%s columns are missing from data: %s.", label, paste(missing, collapse = ", ")),
      label
    )
  }
  columns
}

pv_binding_column_schema <- function(data, columns, component) {
  specs <- lapply(columns, function(name) {
    values <- data[[name]]
    type_id <- pv_binding_validate_predictor(
      values,
      allow_character = TRUE,
      component = component
    )
    factor_levels <- if (is.factor(values)) levels(values) else character()
    list(
      name = name,
      type_id = type_id,
      ordered = is.ordered(values),
      level_count = length(factor_levels),
      ordered_level_labels_hash = pv_binding_ordered_labels_hash(
        factor_levels,
        component,
        "column_levels"
      )
    )
  })
  names(specs) <- NULL
  specs
}

pv_binding_row_identity_digests <- function(data, id_cols = NULL, bound_cols = NULL) {
  if (is.null(id_cols)) {
    id_cols <- character()
  }
  if (!is.character(id_cols) || anyNA(id_cols) || any(!nzchar(id_cols)) || anyDuplicated(id_cols)) {
    pv_binding_abort("PV_BIND_E081", "id_cols must be unique non-empty names when supplied.", "row")
  }
  if (length(id_cols) > 0L) {
    identity_columns <- pv_binding_validate_column_spec(data, id_cols, "row", min_length = 1L)
    identity_mode <- "declared_id"
  } else {
    if (is.null(bound_cols)) {
      pv_binding_abort(
        "PV_BIND_E081",
        "bound_cols are required when declared row IDs are absent.",
        "row"
      )
    }
    identity_columns <- pv_binding_validate_column_spec(data, bound_cols, "row", min_length = 1L)
    identity_mode <- "bound_row_digest"
  }
  schema <- pv_binding_column_schema(data, identity_columns, "row")
  digests <- vapply(seq_len(nrow(data)), function(index) {
    row_values <- lapply(identity_columns, function(name) data[[name]][index])
    names(row_values) <- identity_columns
    pv_binding_hash_payload(
      list(role = "row_identity", columns = identity_columns, values = row_values),
      "row"
    )
  }, character(1))
  if (identity_mode == "declared_id" && anyDuplicated(digests)) {
    pv_binding_abort(
      "PV_BIND_E081",
      "Declared ID columns must jointly identify unique rows.",
      "row"
    )
  }
  list(
    identity_mode = identity_mode,
    identity_columns = identity_columns,
    schema = schema,
    digests = unname(digests)
  )
}

pv_binding_row_projection <- function(data, id_cols = NULL, bound_cols = NULL) {
  identity <- pv_binding_row_identity_digests(data, id_cols = id_cols, bound_cols = bound_cols)
  declared_ids <- if (identity$identity_mode == "declared_id") identity$identity_columns else character()
  list(
    n = nrow(data),
    identity_mode = identity$identity_mode,
    id_cols = declared_ids,
    id_schema = list(
      identity_column_count = length(identity$identity_columns),
      identity_columns = identity$identity_columns,
      columns = identity$schema
    ),
    row_order_hash = pv_binding_hash_payload(
      list(role = "row_order", row_digests = identity$digests),
      "row"
    ),
    row_membership_hash = pv_binding_hash_payload(
      list(role = "row_membership", row_digests = sort(identity$digests, method = "radix")),
      "row"
    ),
    analysis_row_policy = "na.fail"
  )
}

pv_binding_keyed_values_hash <- function(data, columns, row_digests, domain, role) {
  order_index <- order(row_digests, method = "radix")
  rows <- lapply(order_index, function(index) {
    values <- lapply(columns, function(name) data[[name]][index])
    names(values) <- columns
    list(row_digest = row_digests[[index]], values = values)
  })
  pv_binding_hash_payload(list(role = role, keyed_rows = rows), domain)
}

pv_binding_pv_projection <- function(
  data,
  pv_cols,
  id_cols = NULL,
  bound_cols = NULL,
  min_length = 2L
) {
  pv_cols <- pv_binding_validate_column_spec(data, pv_cols, "pv", min_length = min_length)
  for (name in pv_cols) {
    values <- data[[name]]
    if (!is.numeric(values) || any(!is.finite(values))) {
      pv_binding_abort("PV_BIND_E081", "Plausible values must be finite numeric vectors.", "pv")
    }
  }
  identity <- pv_binding_row_identity_digests(data, id_cols = id_cols, bound_cols = bound_cols)
  per_column_hashes <- vapply(pv_cols, function(name) {
    pv_binding_hash_payload(
      list(role = "pv_column_values", name = name,
           values = pv_binding_column_values(data[[name]])),
      "pv"
    )
  }, character(1))
  ordered_payload <- lapply(pv_cols, function(name) {
    list(name = name, values = pv_binding_column_values(data[[name]]))
  })
  list(
    M = length(pv_cols),
    names = pv_cols,
    types = vapply(data[pv_cols], typeof, character(1)),
    per_column_hashes = per_column_hashes,
    ordered_values_hash = pv_binding_hash_payload(
      list(role = "pv_ordered_values", columns = ordered_payload),
      "pv"
    ),
    membership_values_hash = pv_binding_keyed_values_hash(
      data,
      pv_cols,
      identity$digests,
      "pv",
      "pv_membership_values"
    )
  )
}

pv_binding_weight_projection <- function(
  data,
  weight_col,
  rep_weight_cols,
  fay_k = 0.5,
  id_cols = NULL,
  bound_cols = NULL,
  validation_policy_id = "finite_strictly_positive_v1",
  target_transform_id = "raw_wls_weights_v1",
  stack_transform_id = "mean_one_then_divide_by_M_v1"
) {
  weight_col <- pv_binding_validate_column_spec(data, weight_col, "weight", min_length = 1L)
  if (length(weight_col) != 1L) {
    pv_binding_abort("PV_BIND_E081", "Exactly one base weight column is required.", "weight")
  }
  rep_weight_cols <- pv_binding_validate_column_spec(data, rep_weight_cols, "weight", min_length = 2L)
  if (any(c(weight_col, rep_weight_cols) %in% c(weight_col, rep_weight_cols)[duplicated(c(weight_col, rep_weight_cols))])) {
    pv_binding_abort("PV_BIND_E081", "Base and replicate weight names must be distinct.", "weight")
  }
  weight_names <- c(weight_col, rep_weight_cols)
  for (name in weight_names) {
    values <- data[[name]]
    if (!is.numeric(values) || any(!is.finite(values)) || any(values <= 0)) {
      pv_binding_abort("PV_BIND_E081", "Weights must be finite and strictly positive.", "weight")
    }
  }
  if (!is.numeric(fay_k) || length(fay_k) != 1L || !is.finite(fay_k) || fay_k < 0 || fay_k >= 1) {
    pv_binding_abort("PV_BIND_E081", "fay_k must be finite in [0, 1).", "weight")
  }
  expected_ids <- c(
    validation_policy_id = "finite_strictly_positive_v1",
    target_transform_id = "raw_wls_weights_v1",
    stack_transform_id = "mean_one_then_divide_by_M_v1"
  )
  supplied_ids <- c(
    validation_policy_id = validation_policy_id,
    target_transform_id = target_transform_id,
    stack_transform_id = stack_transform_id
  )
  if (!identical(supplied_ids, expected_ids)) {
    pv_binding_abort("PV_BIND_E081", "Weight policy and transform IDs are unsupported.", "weight")
  }
  identity <- pv_binding_row_identity_digests(data, id_cols = id_cols, bound_cols = bound_cols)
  replicate_payload <- lapply(rep_weight_cols, function(name) {
    list(name = name, values = pv_binding_column_values(data[[name]]))
  })
  list(
    base_name = weight_col[[1L]],
    replicate_names = rep_weight_cols,
    replicate_count = length(rep_weight_cols),
    base_values_hash = pv_binding_hash_payload(
      list(role = "base_weight_values", name = weight_col[[1L]],
           values = pv_binding_column_values(data[[weight_col]])),
      "weight"
    ),
    replicate_values_hash = pv_binding_hash_payload(
      list(role = "replicate_weight_values", columns = replicate_payload),
      "weight"
    ),
    membership_values_hash = pv_binding_keyed_values_hash(
      data,
      weight_names,
      identity$digests,
      "weight",
      "weight_membership_values"
    ),
    fay_k = as.double(fay_k),
    validation_policy_id = validation_policy_id,
    target_transform_id = target_transform_id,
    stack_transform_id = stack_transform_id
  )
}

pv_binding_factor_contrast_component <- function(
  data,
  formula,
  model_bundle = NULL,
  model_projection = NULL
) {
  if (is.null(model_bundle)) {
    model_bundle <- pv_binding_resolve_model_bundle(data, formula)
  }
  pv_binding_validate_model_bundle(model_bundle, data, formula)
  if (is.null(model_projection)) {
    model_projection <- pv_binding_model_matrix_projection(
      data,
      formula,
      model_bundle = model_bundle
    )
  }
  if (!identical(model_projection$contrasts, model_bundle$contrasts)) {
    pv_binding_abort(
      "PV_BIND_E081",
      "Factor-contrast projection does not match the resolved model bundle.",
      "factor_contrast"
    )
  }
  variables <- vapply(model_projection$contrasts, `[[`, character(1), "variable")
  option_labels <- model_bundle$contrast_option_labels
  option_snapshot <- list(
    count = length(option_labels),
    ordered_label_hash = pv_binding_ordered_labels_hash(
      option_labels,
      "factor_contrast",
      "global_contrast_options"
    )
  )
  resolved_payload <- list(
    variables = variables,
    option_snapshot = option_snapshot,
    resolved_specs = model_projection$contrasts
  )
  list(
    variables = variables,
    option_snapshot = option_snapshot,
    resolved_specs = model_projection$contrasts,
    resolved_hash = pv_binding_hash_payload(resolved_payload, "factor_contrast")
  )
}

pv_binding_validate_string_vector <- function(x, label, min_length = 1L) {
  if (!is.character(x) || length(x) < min_length || anyNA(x) || any(!nzchar(x)) ||
      any(!validUTF8(enc2utf8(x))) || any(trimws(x) != x) || anyDuplicated(x)) {
    pv_binding_abort(
      "PV_BIND_E081",
      sprintf("%s must contain unique non-empty names.", label),
      "estimand"
    )
  }
  x
}

pv_binding_validate_metadata_string <- function(x, label, component) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x) ||
      !validUTF8(enc2utf8(x)) || trimws(x) != x) {
    pv_binding_abort(
      "PV_BIND_E081",
      sprintf("%s must be a trimmed non-empty UTF-8 string.", label),
      component
    )
  }
  x
}

pv_binding_validate_estimand_contrast_component <- function(component) {
  required <- c(
    "contrast_id", "transformation_class", "input_fe_count",
    "ordered_input_fe_names_hash", "output_count", "ordered_output_names_hash",
    "matrix_dimensions", "matrix_hash", "coverage_inheritance"
  )
  if (!is.list(component) || !identical(names(component), required)) {
    pv_binding_abort("PV_BIND_E081", "Estimand contrast fields must exactly match the binding schema.", "estimand_contrast")
  }
  pv_binding_validate_metadata_string(component$contrast_id, "contrast_id", "estimand_contrast")
  if (!identical(component$transformation_class, "affine")) {
    pv_binding_abort("PV_BIND_E081", "Estimand contrast transformation_class must be affine.", "estimand_contrast")
  }
  if (!is.integer(component$input_fe_count) || length(component$input_fe_count) != 1L ||
      is.na(component$input_fe_count) || component$input_fe_count < 1L ||
      !is.integer(component$output_count) || length(component$output_count) != 1L ||
      is.na(component$output_count) || component$output_count < 1L) {
    pv_binding_abort("PV_BIND_E081", "Estimand contrast input/output counts must be positive integers.", "estimand_contrast")
  }
  for (field in c("ordered_input_fe_names_hash", "ordered_output_names_hash")) {
    if (!is.character(component[[field]]) || length(component[[field]]) != 1L ||
        is.na(component[[field]]) || !grepl("^sha256:[0-9a-f]{64}$", component[[field]])) {
      pv_binding_abort("PV_BIND_E081", "Estimand contrast ordered-name hashes must be canonical SHA-256 digests.", "estimand_contrast")
    }
  }
  dimensions <- component$matrix_dimensions
  if (!is.integer(dimensions) || length(dimensions) != 2L || anyNA(dimensions) || any(dimensions < 1L) ||
      !identical(dimensions, c(component$output_count, component$input_fe_count))) {
    pv_binding_abort("PV_BIND_E081", "Estimand contrast dimensions must align with output and input names.", "estimand_contrast")
  }
  if (!is.character(component$matrix_hash) || length(component$matrix_hash) != 1L ||
      is.na(component$matrix_hash) || !grepl("^sha256:[0-9a-f]{64}$", component$matrix_hash)) {
    pv_binding_abort("PV_BIND_E081", "Estimand contrast matrix_hash must be a canonical SHA-256 digest.", "estimand_contrast")
  }
  if (!is.logical(component$coverage_inheritance) || length(component$coverage_inheritance) != 1L ||
      is.na(component$coverage_inheritance)) {
    pv_binding_abort("PV_BIND_E081", "Estimand contrast coverage_inheritance must be a logical scalar.", "estimand_contrast")
  }
  invisible(component)
}

pv_binding_estimand_contrast_projection <- function(
  fe_names,
  L = NULL,
  contrast_id = NULL,
  output_names = NULL,
  transformation_class = "affine",
  coverage_inheritance = TRUE
) {
  fe_names <- pv_binding_validate_string_vector(fe_names, "fe_names")
  if (is.null(L)) {
    L <- diag(length(fe_names))
    dimnames(L) <- list(fe_names, fe_names)
    if (is.null(contrast_id)) {
      contrast_id <- "identity_fixed_effect_v1"
    }
  }
  if (!is.matrix(L) || !is.numeric(L) || any(!is.finite(L)) || ncol(L) != length(fe_names)) {
    pv_binding_abort("PV_BIND_E081", "Estimand contrast L must be a finite numeric matrix aligned to fixed effects.", "estimand_contrast")
  }
  if (is.null(colnames(L)) || !identical(colnames(L), fe_names)) {
    pv_binding_abort("PV_BIND_E081", "Estimand contrast columns must exactly match fe_names in order.", "estimand_contrast")
  }
  if (is.null(output_names)) {
    output_names <- rownames(L)
  }
  output_names <- pv_binding_validate_string_vector(output_names, "output_names")
  if (length(output_names) != nrow(L) || (!is.null(rownames(L)) && !identical(rownames(L), output_names))) {
    pv_binding_abort("PV_BIND_E081", "Estimand contrast rows must exactly match output_names.", "estimand_contrast")
  }
  if (!is.character(contrast_id) || length(contrast_id) != 1L || is.na(contrast_id) || !nzchar(contrast_id)) {
    pv_binding_abort("PV_BIND_E081", "contrast_id must be a non-empty string.", "estimand_contrast")
  }
  if (!identical(transformation_class, "affine")) {
    pv_binding_abort("PV_BIND_E081", "Only affine estimand contrasts are supported by this binding component.", "estimand_contrast")
  }
  if (!is.logical(coverage_inheritance) || length(coverage_inheritance) != 1L || is.na(coverage_inheritance)) {
    pv_binding_abort("PV_BIND_E081", "coverage_inheritance must be a non-missing logical scalar.", "estimand_contrast")
  }
  list(
    contrast_id = contrast_id,
    transformation_class = transformation_class,
    input_fe_count = as.integer(length(fe_names)),
    ordered_input_fe_names_hash = pv_binding_ordered_labels_hash(
      fe_names,
      "estimand_contrast",
      "input_fe_names"
    ),
    output_count = as.integer(length(output_names)),
    ordered_output_names_hash = pv_binding_ordered_labels_hash(
      output_names,
      "estimand_contrast",
      "output_names"
    ),
    matrix_dimensions = as.integer(dim(L)),
    matrix_hash = pv_binding_numeric_matrix_hash(L, "estimand_contrast", "affine_L"),
    coverage_inheritance = coverage_inheritance
  ) |> pv_binding_validate_estimand_contrast_component()
}

pv_binding_family_link_projection <- function(
  family_id,
  link_id,
  response_support_id = "real",
  dispersion_role = "estimated"
) {
  fields <- list(
    family_id = pv_binding_validate_metadata_string(family_id, "family_id", "family_link"),
    link_id = pv_binding_validate_metadata_string(link_id, "link_id", "family_link"),
    response_support_id = pv_binding_validate_metadata_string(response_support_id, "response_support_id", "family_link"),
    dispersion_role = pv_binding_validate_metadata_string(dispersion_role, "dispersion_role", "family_link")
  )
  fields
}

pv_binding_estimand_projection <- function(
  fe_names,
  estimand_contrast,
  interval_role,
  estimand_id = "brr_fay_fixed_effects",
  target_source = "external_brr_fay_rubin",
  target_engine_id = "lm_wls_brr_fay_v1",
  parameter_scope = "fixed_effect",
  coverage_claim_allowed = FALSE
) {
  fe_names <- pv_binding_validate_string_vector(fe_names, "fe_names")
  pv_binding_validate_estimand_contrast_component(estimand_contrast)
  expected_input_hash <- pv_binding_ordered_labels_hash(
    fe_names,
    "estimand_contrast",
    "input_fe_names"
  )
  if (!identical(estimand_contrast$input_fe_count, as.integer(length(fe_names))) ||
      !identical(estimand_contrast$ordered_input_fe_names_hash, expected_input_hash)) {
    pv_binding_abort("PV_BIND_E081", "estimand_contrast must be a canonical affine contrast aligned to fe_names.", "estimand")
  }
  string_fields <- list(
    estimand_id = pv_binding_validate_metadata_string(estimand_id, "estimand_id", "estimand"),
    target_source = pv_binding_validate_metadata_string(target_source, "target_source", "estimand"),
    target_engine_id = pv_binding_validate_metadata_string(target_engine_id, "target_engine_id", "estimand"),
    parameter_scope = pv_binding_validate_metadata_string(parameter_scope, "parameter_scope", "estimand"),
    interval_role = pv_binding_validate_metadata_string(interval_role, "interval_role", "estimand")
  )
  if (!is.logical(coverage_claim_allowed) || length(coverage_claim_allowed) != 1L || is.na(coverage_claim_allowed)) {
    pv_binding_abort("PV_BIND_E081", "coverage_claim_allowed must be a non-missing logical scalar.", "estimand")
  }
  if (coverage_claim_allowed && !isTRUE(estimand_contrast$coverage_inheritance)) {
    pv_binding_abort("PV_BIND_E081", "Coverage cannot be inherited through a non-inheriting contrast.", "estimand")
  }
  c(
    string_fields["estimand_id"],
    string_fields["target_source"],
    string_fields["target_engine_id"],
    string_fields["parameter_scope"],
    list(fe_names = fe_names),
    string_fields["interval_role"],
    list(
      coverage_claim_allowed = coverage_claim_allowed,
      estimand_contrast_hash = pv_binding_hash_payload(estimand_contrast, "estimand_contrast")
    )
  )
}

pv_binding_component_hash <- function(component, component_name) {
  if (!is.list(component) || !is.character(component_name) || length(component_name) != 1L ||
      is.na(component_name) || !component_name %in% pv_binding_hash_domains()) {
    pv_binding_abort("PV_BIND_E081", "A component hash requires a named binding domain and list payload.", "manifest")
  }
  pv_binding_hash_payload(component, component_name)
}

pv_binding_formula_component <- function(data, formula, model_projection = NULL) {
  if (is.null(model_projection)) {
    model_projection <- pv_binding_model_matrix_projection(data, formula)
  }
  ast <- pv_binding_formula_ast(formula, data = data)
  specials <- if (isTRUE(model_projection$has_offset)) "offset" else character()
  list(
    outcome_placeholder = "OUTCOME",
    data_variable_count = length(ast$data_variables),
    data_variables = ast$data_variables,
    rhs_ast_hash = pv_binding_hash_payload(ast$rhs, "formula"),
    term_count = model_projection$term_count,
    ordered_term_labels_hash = model_projection$ordered_term_labels_hash,
    intercept = model_projection$intercept,
    specials_count = length(specials),
    ordered_specials_hash = pv_binding_ordered_labels_hash(
      specials,
      "formula",
      "specials"
    ),
    offset_hash = model_projection$offset_values_hash
  )
}

pv_binding_component_names <- function() {
  c(
    "row", "pv", "predictor", "formula", "model_matrix", "weight",
    "factor_contrast", "estimand_contrast", "family_link", "estimand"
  )
}

pv_binding_component_hash_names <- function() {
  paste0(pv_binding_component_names(), "_hash")
}

pv_binding_component_field_registry <- function() {
  list(
    row = c(
      "n", "identity_mode", "id_cols", "id_schema", "row_order_hash",
      "row_membership_hash", "analysis_row_policy"
    ),
    pv = c(
      "M", "names", "types", "per_column_hashes", "ordered_values_hash",
      "membership_values_hash"
    ),
    predictor = c(
      "source_rule", "row_count", "column_count", "names", "columns",
      "schema_hash", "ordered_values_hash", "membership_values_hash"
    ),
    formula = c(
      "outcome_placeholder", "data_variable_count", "data_variables",
      "rhs_ast_hash", "term_count", "ordered_term_labels_hash", "intercept",
      "specials_count", "ordered_specials_hash", "offset_hash"
    ),
    model_matrix = c(
      "builder_id", "na_action", "row_count", "column_count", "fe_colname_count",
      "ordered_fe_colnames_hash", "ordered_reportable_fe_names_hash", "assign", "intercept",
      "term_count", "ordered_term_labels_hash", "has_offset", "offset_values_hash",
      "xlevel_count", "xlevels", "contrast_count", "contrasts", "values_hash",
      "membership_values_hash", "offset_membership_hash", "predictor_schema_hash",
      "predictor_values_hash"
    ),
    weight = c(
      "base_name", "replicate_names", "replicate_count", "base_values_hash",
      "replicate_values_hash", "membership_values_hash", "fay_k",
      "validation_policy_id", "target_transform_id", "stack_transform_id"
    ),
    factor_contrast = c("variables", "option_snapshot", "resolved_specs", "resolved_hash"),
    estimand_contrast = c(
      "contrast_id", "transformation_class", "input_fe_count",
      "ordered_input_fe_names_hash", "output_count", "ordered_output_names_hash",
      "matrix_dimensions", "matrix_hash", "coverage_inheritance"
    ),
    family_link = c("family_id", "link_id", "response_support_id", "dispersion_role"),
    estimand = c(
      "estimand_id", "target_source", "target_engine_id", "parameter_scope",
      "fe_names", "interval_role", "coverage_claim_allowed", "estimand_contrast_hash"
    )
  )
}

pv_binding_manifest_root_fields <- function() {
  c(
    "contract_id", "manifest_schema_version", "canonicalizer_id",
    "hash_algorithm_id", "binding_scope", "components", "component_hashes",
    "manifest_hash"
  )
}

pv_binding_manifest_constants <- function() {
  list(
    contract_id = "pvstackr_data_binding_v1",
    manifest_schema_version = "0.2.0",
    canonicalizer_id = "pvstackr_c14n_tlv_v1",
    hash_algorithm_id = "sha256",
    binding_scope = "target_fit_compatibility"
  )
}

pv_binding_assert_exact_fields <- function(x, fields, label, code = "PV_BIND_E002") {
  if (!is.list(x) || !identical(names(x), fields) ||
      !identical(attributes(x), list(names = fields))) {
    pv_binding_abort(
      code,
      sprintf("%s fields and order do not exactly match the binding schema.", label),
      label
    )
  }
  invisible(x)
}

pv_binding_assert_manifest_attributes_recursive <- function(x, label = "manifest") {
  if (is.null(x)) {
    return(invisible(x))
  }
  value_names <- names(x)
  expected_attributes <- if (is.null(value_names)) NULL else list(names = value_names)
  if (!identical(attributes(x), expected_attributes)) {
    pv_binding_abort(
      "PV_BIND_E002",
      sprintf("%s contains attributes outside its exact binding schema.", label),
      "manifest"
    )
  }
  if (is.list(x)) {
    for (index in seq_along(x)) {
      child_label <- if (is.null(value_names)) {
        sprintf("%s item %d", label, index)
      } else {
        sprintf("%s.%s", label, value_names[[index]])
      }
      pv_binding_assert_manifest_attributes_recursive(x[[index]], child_label)
    }
  }
  invisible(x)
}

pv_binding_assert_manifest_atomic_names <- function(x, label = "manifest") {
  if (is.list(x)) {
    for (index in seq_along(x)) {
      pv_binding_assert_manifest_atomic_names(x[[index]], label)
    }
  } else if (!is.null(names(x))) {
    pv_binding_abort(
      "PV_BIND_E002",
      sprintf("%s contains names on an atomic field that is not keyed by schema.", label),
      "manifest"
    )
  }
  invisible(x)
}

pv_binding_assert_unnamed_manifest_list <- function(x, label) {
  if (!is.list(x) || !is.null(attributes(x))) {
    pv_binding_abort(
      "PV_BIND_E002",
      sprintf("%s must be an unnamed attribute-free list.", label),
      "manifest"
    )
  }
  invisible(x)
}

pv_binding_validate_manifest_attribute_schema <- function(manifest) {
  pv_binding_assert_manifest_attributes_recursive(manifest)

  atomic_view <- manifest
  atomic_view$components$pv$types <- unname(atomic_view$components$pv$types)
  atomic_view$components$pv$per_column_hashes <- unname(
    atomic_view$components$pv$per_column_hashes
  )
  pv_binding_assert_manifest_atomic_names(atomic_view)

  containers <- list(
    "row id-schema columns" = manifest$components$row$id_schema$columns,
    "predictor columns" = manifest$components$predictor$columns,
    "model-matrix xlevels" = manifest$components$model_matrix$xlevels,
    "model-matrix contrasts" = manifest$components$model_matrix$contrasts,
    "factor-contrast resolved specs" = manifest$components$factor_contrast$resolved_specs
  )
  for (label in names(containers)) {
    pv_binding_assert_unnamed_manifest_list(containers[[label]], label)
  }
  invisible(manifest)
}

pv_binding_manifest_hash_payload <- function(manifest) {
  pv_binding_assert_manifest_attributes_recursive(manifest)
  fields <- pv_binding_manifest_root_fields()
  payload <- manifest[setdiff(fields, "manifest_hash")]
  if (!is.null(manifest$model_bundle_hash)) {
    payload$model_bundle_hash <- manifest$model_bundle_hash
  }
  if (!is.null(manifest$migration)) {
    migration <- manifest$migration
    if (!is.list(migration)) {
      pv_binding_abort("PV_BIND_E002", "Manifest migration metadata must be a list.", "manifest")
    }
    migration$warnings <- NULL
    payload$migration <- migration
  }
  payload
}

pv_binding_legacy_hash_algorithm_id <- function() {
  "adler32_r_serialize_v2"
}

pv_binding_validate_optional_legacy_hashes <- function(x) {
  allowed <- c("algorithm_id", "design_hash", "target_hash")
  if (!is.list(x) || !identical(attributes(x), list(names = names(x))) ||
      is.null(names(x)) || anyNA(names(x)) || any(!nzchar(names(x))) ||
      anyDuplicated(names(x)) || !all(names(x) %in% allowed) ||
      !identical(names(x), allowed[allowed %in% names(x)]) ||
      !"algorithm_id" %in% names(x) ||
      !any(c("design_hash", "target_hash") %in% names(x))) {
    pv_binding_abort(
      "PV_BIND_E002",
      "legacy_hashes must contain algorithm_id followed by applicable design_hash and target_hash fields.",
      "manifest"
    )
  }
  if (!is.character(x$algorithm_id) || length(x$algorithm_id) != 1L ||
      !is.null(attributes(x$algorithm_id)) ||
      !identical(x$algorithm_id, pv_binding_legacy_hash_algorithm_id())) {
    pv_binding_abort(
      "PV_BIND_E002",
      "legacy_hashes algorithm_id is unsupported.",
      "manifest"
    )
  }
  hash_fields <- intersect(c("design_hash", "target_hash"), names(x))
  valid <- vapply(x[hash_fields], function(value) {
    is.character(value) && length(value) == 1L && is.null(attributes(value)) &&
      !is.na(value) && grepl("^[0-9a-f]{8}$", value)
  }, logical(1))
  if (!all(valid)) {
    pv_binding_abort(
      "PV_BIND_E002",
      "legacy_hashes design_hash and target_hash values must be attribute-free scalar 8-hex strings.",
      "manifest"
    )
  }
  invisible(x)
}

pv_binding_validate_optional_created_at <- function(x) {
  if (!is.character(x) || length(x) != 1L || !is.null(attributes(x)) || is.na(x) ||
      !grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}(\\.[0-9]{1,6})?Z$", x)) {
    pv_binding_abort("PV_BIND_E002", "created_at must be one UTC ISO-8601 timestamp.", "manifest")
  }
  whole_second <- sub("\\.[0-9]{1,6}Z$", "Z", x)
  parsed <- suppressWarnings(as.POSIXct(
    whole_second,
    format = "%Y-%m-%dT%H:%M:%SZ",
    tz = "UTC"
  ))
  roundtrip <- if (length(parsed) == 1L && !is.na(parsed)) {
    format(parsed, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC", usetz = FALSE)
  } else {
    NA_character_
  }
  if (is.na(roundtrip) || !identical(roundtrip, whole_second)) {
    pv_binding_abort(
      "PV_BIND_E002",
      "created_at must name a real UTC calendar date and time.",
      "manifest"
    )
  }
  invisible(x)
}

pv_binding_migration_from_schema_ids <- function() {
  "0.1.0"
}

pv_binding_migration_function_ids <- function() {
  "pv_binding_revalidate_brr_target_v1"
}

pv_binding_migration_warning_codes <- function() {
  c(
    "PV_BIND_MIGRATION_LEGACY_HASH_ONLY",
    "PV_BIND_MIGRATION_SOURCE_METADATA_DROPPED"
  )
}

pv_binding_validate_optional_migration <- function(x) {
  fields <- c(
    "migration_from_schema", "migration_function", "binding_revalidated",
    "inspection_only", "warnings"
  )
  pv_binding_assert_exact_fields(x, fields, "migration")
  for (field in c("migration_from_schema", "migration_function")) {
    value <- x[[field]]
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !is.null(attributes(value))) {
      pv_binding_abort(
        "PV_BIND_E002",
        "Migration schema and function must be attribute-free scalar registered identifiers.",
        "manifest"
      )
    }
  }
  pv_binding_validate_metadata_string(x$migration_from_schema, "migration_from_schema", "manifest")
  pv_binding_validate_metadata_string(x$migration_function, "migration_function", "manifest")
  if (!x$migration_from_schema %in% pv_binding_migration_from_schema_ids() ||
      !x$migration_function %in% pv_binding_migration_function_ids()) {
    pv_binding_abort(
      "PV_BIND_E002",
      "Migration schema and function must be exact registered identifiers.",
      "manifest"
    )
  }
  if (!is.logical(x$binding_revalidated) || length(x$binding_revalidated) != 1L ||
      is.na(x$binding_revalidated) || !is.logical(x$inspection_only) ||
      length(x$inspection_only) != 1L || is.na(x$inspection_only) ||
      identical(x$binding_revalidated, x$inspection_only)) {
    pv_binding_abort(
      "PV_BIND_E002",
      "Migration must be exactly one of revalidated-reportable or inspection-only.",
      "manifest"
    )
  }
  warning_codes <- pv_binding_migration_warning_codes()
  if (!is.character(x$warnings) || !is.null(attributes(x$warnings)) ||
      anyNA(x$warnings) || anyDuplicated(x$warnings) ||
      !all(x$warnings %in% warning_codes) ||
      !identical(x$warnings, warning_codes[warning_codes %in% x$warnings])) {
    pv_binding_abort(
      "PV_BIND_E002",
      "Migration warnings must be empty or an ordered unique subset of registered warning codes.",
      "manifest"
    )
  }
  invisible(x)
}

pv_binding_validate_optional_manifest_fields <- function(manifest) {
  if ("legacy_hashes" %in% names(manifest)) {
    pv_binding_validate_optional_legacy_hashes(manifest$legacy_hashes)
  }
  if ("created_at" %in% names(manifest)) {
    pv_binding_validate_optional_created_at(manifest$created_at)
  }
  if ("model_bundle_hash" %in% names(manifest)) {
    pv_binding_schema_hash(
      manifest$model_bundle_hash,
      "model_bundle_hash",
      "manifest"
    )
  }
  if ("migration" %in% names(manifest)) {
    pv_binding_validate_optional_migration(manifest$migration)
  }
  invisible(manifest)
}

pv_binding_manifest_assert_reportable <- function(manifest) {
  pv_binding_manifest_validate(manifest)
  migration <- manifest$migration
  if (!is.null(migration) &&
      (!identical(migration$binding_revalidated, TRUE) ||
        !identical(migration$inspection_only, FALSE))) {
    pv_binding_abort(
      "PV_BIND_E080",
      "Inspection-only migrated targets are not reportable until raw-input revalidation succeeds.",
      "manifest"
    )
  }
  invisible(manifest)
}

pv_binding_stack_long_data_hash <- function(manifest) {
  pv_binding_manifest_validate(manifest)
  model_matrix <- manifest$components$model_matrix
  internal_columns <- sprintf(
    "pvstackrMM%03d",
    seq_len(model_matrix$column_count)
  )
  hash_columns <- c(
    ".pvstackr_y", ".pvstackr_pv", ".pvstackr_row", ".pvstackr_weight",
    internal_columns,
    if (isTRUE(model_matrix$has_offset)) "pvstackrOffset" else NULL
  )
  relevant_components <- c(
    "row_hash", "pv_hash", "model_matrix_hash", "weight_hash"
  )
  payload <- list(
    contract_id = manifest$contract_id,
    topology = "single_long_fit",
    materialization_id = "authenticated_component_projection_v1",
    component_hashes = manifest$component_hashes[relevant_components],
    hash_columns = hash_columns
  )
  substr(
    sub(
      "^sha256:",
      "",
      pv_binding_hash_payload(payload, "manifest")
    ),
    1L,
    8L
  )
}

pv_binding_schema_abort <- function(detail, component) {
  pv_binding_abort("PV_BIND_E002", detail, component)
}

pv_binding_schema_count <- function(x, label, component, minimum = 0L) {
  if (!is.integer(x) || length(x) != 1L || is.na(x) || x < minimum) {
    pv_binding_schema_abort(sprintf("%s must be an integer count.", label), component)
  }
  x
}

pv_binding_schema_hash <- function(x, label, component) {
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      !grepl("^sha256:[0-9a-f]{64}$", x)) {
    pv_binding_schema_abort(sprintf("%s must be a canonical SHA-256 digest.", label), component)
  }
  x
}

pv_binding_schema_names <- function(x, label, component, minimum = 0L) {
  if (!is.character(x) || length(x) < minimum || anyNA(x) || any(!nzchar(x)) ||
      any(trimws(x) != x) || any(!validUTF8(enc2utf8(x))) || anyDuplicated(x)) {
    pv_binding_schema_abort(sprintf("%s must contain unique trimmed UTF-8 names.", label), component)
  }
  x
}

pv_binding_validate_level_schema <- function(spec, fields, component, allow_character) {
  pv_binding_assert_exact_fields(spec, fields, paste0(component, " column"))
  pv_binding_validate_metadata_string(spec$name, "column name", component)
  allowed <- c("logical", "integer", "double", "factor", "ordered_factor", "Date", "POSIXct")
  if (allow_character) {
    allowed <- c(allowed, "character")
  }
  if (!is.character(spec$type_id) || length(spec$type_id) != 1L ||
      !spec$type_id %in% allowed || !is.logical(spec$ordered) ||
      length(spec$ordered) != 1L || is.na(spec$ordered) ||
      !identical(spec$ordered, identical(spec$type_id, "ordered_factor"))) {
    pv_binding_schema_abort("Column type and ordered metadata are inconsistent.", component)
  }
  pv_binding_schema_count(spec$level_count, "level_count", component)
  is_factor <- spec$type_id %in% c("factor", "ordered_factor")
  if ((is_factor && spec$level_count < 1L) || (!is_factor && spec$level_count != 0L)) {
    pv_binding_schema_abort("Column level count is inconsistent with its type.", component)
  }
  pv_binding_schema_hash(spec$ordered_level_labels_hash, "ordered_level_labels_hash", component)
  if ("timezone_hash" %in% fields) {
    pv_binding_schema_hash(spec$timezone_hash, "timezone_hash", component)
  }
  if ("values_hash" %in% fields) {
    pv_binding_schema_hash(spec$values_hash, "values_hash", component)
  }
  invisible(spec)
}

pv_binding_validate_contrast_spec_schema <- function(spec, component = "factor_contrast") {
  fields <- c(
    "variable", "generator_label_hash", "dimensions", "row_label_count",
    "row_label_hash", "column_label_count", "column_label_hash", "values_hash"
  )
  pv_binding_assert_exact_fields(spec, fields, "resolved contrast")
  pv_binding_validate_metadata_string(spec$variable, "contrast variable", component)
  for (field in c("generator_label_hash", "row_label_hash", "column_label_hash", "values_hash")) {
    pv_binding_schema_hash(spec[[field]], field, component)
  }
  if (!is.integer(spec$dimensions) || length(spec$dimensions) != 2L ||
      anyNA(spec$dimensions) || any(spec$dimensions < 0L)) {
    pv_binding_schema_abort("Contrast dimensions must be two non-negative integers.", component)
  }
  pv_binding_schema_count(spec$row_label_count, "row_label_count", component)
  pv_binding_schema_count(spec$column_label_count, "column_label_count", component)
  if (spec$row_label_count > spec$dimensions[[1L]] ||
      spec$column_label_count > spec$dimensions[[2L]]) {
    pv_binding_schema_abort("Contrast label counts cannot exceed dimensions.", component)
  }
  invisible(spec)
}

pv_binding_validate_component_deep <- function(component, component_name) {
  if (identical(component_name, "row")) {
    pv_binding_schema_count(component$n, "n", "row", minimum = 1L)
    if (!is.character(component$identity_mode) || length(component$identity_mode) != 1L ||
        !component$identity_mode %in% c("declared_id", "bound_row_digest")) {
      pv_binding_schema_abort("row identity_mode is unsupported.", "row")
    }
    pv_binding_schema_names(component$id_cols, "id_cols", "row")
    pv_binding_assert_exact_fields(
      component$id_schema,
      c("identity_column_count", "identity_columns", "columns"),
      "row id_schema"
    )
    count <- pv_binding_schema_count(
      component$id_schema$identity_column_count,
      "identity_column_count",
      "row",
      minimum = 1L
    )
    columns <- pv_binding_schema_names(
      component$id_schema$identity_columns,
      "identity_columns",
      "row",
      minimum = 1L
    )
    if (length(columns) != count || !is.list(component$id_schema$columns) ||
        length(component$id_schema$columns) != count) {
      pv_binding_schema_abort("row identity schema counts do not align.", "row")
    }
    fields <- c("name", "type_id", "ordered", "level_count", "ordered_level_labels_hash")
    for (spec in component$id_schema$columns) {
      pv_binding_validate_level_schema(spec, fields, "row", allow_character = TRUE)
    }
    spec_names <- vapply(component$id_schema$columns, `[[`, character(1), "name")
    if (!identical(spec_names, columns) ||
        (component$identity_mode == "declared_id" && !identical(component$id_cols, columns)) ||
        (component$identity_mode == "bound_row_digest" && length(component$id_cols) != 0L)) {
      pv_binding_schema_abort("row identity fields do not agree in order.", "row")
    }
    pv_binding_schema_hash(component$row_order_hash, "row_order_hash", "row")
    pv_binding_schema_hash(component$row_membership_hash, "row_membership_hash", "row")
    if (!identical(component$analysis_row_policy, "na.fail")) {
      pv_binding_schema_abort("row analysis policy must be na.fail.", "row")
    }
  } else if (identical(component_name, "pv")) {
    M <- pv_binding_schema_count(component$M, "M", "pv", minimum = 1L)
    declared_names <- pv_binding_schema_names(component$names, "PV names", "pv", minimum = 1L)
    if (length(declared_names) != M || !is.character(component$types) ||
        !identical(names(component$types), declared_names) || length(component$types) != M ||
        any(!component$types %in% c("integer", "double")) ||
        !is.character(component$per_column_hashes) ||
        !identical(names(component$per_column_hashes), declared_names) ||
        length(component$per_column_hashes) != M) {
      pv_binding_schema_abort("PV names, types, and per-column hashes must align to M.", "pv")
    }
    for (hash in unname(component$per_column_hashes)) {
      pv_binding_schema_hash(hash, "per_column_hash", "pv")
    }
    pv_binding_schema_hash(component$ordered_values_hash, "ordered_values_hash", "pv")
    pv_binding_schema_hash(component$membership_values_hash, "membership_values_hash", "pv")
  } else if (identical(component_name, "predictor")) {
    if (!identical(component$source_rule, "rhs_ast_first_occurrence_v1")) {
      pv_binding_schema_abort("predictor source rule is unsupported.", "predictor")
    }
    pv_binding_schema_count(component$row_count, "row_count", "predictor", minimum = 1L)
    count <- pv_binding_schema_count(component$column_count, "column_count", "predictor")
    declared_names <- pv_binding_schema_names(component$names, "predictor names", "predictor")
    if (length(declared_names) != count || !is.list(component$columns) || length(component$columns) != count) {
      pv_binding_schema_abort("predictor columns must align to column_count.", "predictor")
    }
    fields <- c(
      "name", "type_id", "ordered", "level_count", "ordered_level_labels_hash",
      "timezone_hash", "values_hash"
    )
    for (spec in component$columns) {
      pv_binding_validate_level_schema(spec, fields, "predictor", allow_character = FALSE)
    }
    if (count > 0L && !identical(vapply(component$columns, `[[`, character(1), "name"), declared_names)) {
      pv_binding_schema_abort("predictor column names are not in declared order.", "predictor")
    }
    for (field in c("schema_hash", "ordered_values_hash", "membership_values_hash")) {
      pv_binding_schema_hash(component[[field]], field, "predictor")
    }
    schema_payload <- lapply(component$columns, function(column) {
      column[setdiff(names(column), "values_hash")]
    })
    if (!identical(
      component$schema_hash,
      pv_binding_hash_payload(schema_payload, "predictor")
    )) {
      pv_binding_schema_abort("predictor schema_hash is not canonical for its column schema.", "predictor")
    }
  } else if (identical(component_name, "formula")) {
    if (!identical(component$outcome_placeholder, "OUTCOME")) {
      pv_binding_schema_abort("formula outcome placeholder must be OUTCOME.", "formula")
    }
    variable_count <- pv_binding_schema_count(component$data_variable_count, "data_variable_count", "formula")
    variables <- pv_binding_schema_names(component$data_variables, "data_variables", "formula")
    if (length(variables) != variable_count) {
      pv_binding_schema_abort("formula data variables do not align to count.", "formula")
    }
    pv_binding_schema_count(component$term_count, "term_count", "formula")
    if (!is.integer(component$intercept) || length(component$intercept) != 1L ||
        is.na(component$intercept) || !component$intercept %in% c(0L, 1L)) {
      pv_binding_schema_abort("formula intercept must be integer 0 or 1.", "formula")
    }
    pv_binding_schema_count(component$specials_count, "specials_count", "formula")
    for (field in c(
      "rhs_ast_hash", "ordered_term_labels_hash", "ordered_specials_hash", "offset_hash"
    )) {
      pv_binding_schema_hash(component[[field]], field, "formula")
    }
  } else if (identical(component_name, "model_matrix")) {
    if (!identical(component$builder_id, "stats_model_frame_matrix_v1") ||
        !identical(component$na_action, "na.fail")) {
      pv_binding_schema_abort("model-matrix builder and NA policy are unsupported.", "model_matrix")
    }
    row_count <- pv_binding_schema_count(component$row_count, "row_count", "model_matrix", minimum = 1L)
    column_count <- pv_binding_schema_count(component$column_count, "column_count", "model_matrix", minimum = 1L)
    fe_count <- pv_binding_schema_count(component$fe_colname_count, "fe_colname_count", "model_matrix", minimum = 1L)
    if (column_count != fe_count || !is.integer(component$assign) ||
        length(component$assign) != column_count || anyNA(component$assign) || any(component$assign < 0L)) {
      pv_binding_schema_abort("model-matrix column counts and assign vector do not align.", "model_matrix")
    }
    if (!is.integer(component$intercept) || length(component$intercept) != 1L ||
        !component$intercept %in% c(0L, 1L) || !is.logical(component$has_offset) ||
        length(component$has_offset) != 1L || is.na(component$has_offset)) {
      pv_binding_schema_abort("model-matrix intercept/offset metadata is malformed.", "model_matrix")
    }
    pv_binding_schema_count(component$term_count, "term_count", "model_matrix")
    xlevel_count <- pv_binding_schema_count(component$xlevel_count, "xlevel_count", "model_matrix")
    contrast_count <- pv_binding_schema_count(component$contrast_count, "contrast_count", "model_matrix")
    if (!is.list(component$xlevels) || length(component$xlevels) != xlevel_count ||
        !is.list(component$contrasts) || length(component$contrasts) != contrast_count) {
      pv_binding_schema_abort("model-matrix xlevel/contrast counts do not align.", "model_matrix")
    }
    xlevel_fields <- c("variable", "ordered", "level_count", "ordered_level_labels_hash")
    for (spec in component$xlevels) {
      pv_binding_assert_exact_fields(spec, xlevel_fields, "xlevel spec")
      pv_binding_validate_metadata_string(spec$variable, "xlevel variable", "model_matrix")
      if (!is.logical(spec$ordered) || length(spec$ordered) != 1L || is.na(spec$ordered)) {
        pv_binding_schema_abort("xlevel ordered flag must be logical.", "model_matrix")
      }
      pv_binding_schema_count(spec$level_count, "xlevel level_count", "model_matrix", minimum = 1L)
      pv_binding_schema_hash(spec$ordered_level_labels_hash, "xlevel labels hash", "model_matrix")
    }
    for (spec in component$contrasts) {
      pv_binding_validate_contrast_spec_schema(spec, "model_matrix")
    }
    for (field in c(
      "ordered_fe_colnames_hash", "ordered_reportable_fe_names_hash",
      "ordered_term_labels_hash", "offset_values_hash", "values_hash",
      "membership_values_hash", "offset_membership_hash", "predictor_schema_hash",
      "predictor_values_hash"
    )) {
      pv_binding_schema_hash(component[[field]], field, "model_matrix")
    }
    invisible(row_count)
  } else if (identical(component_name, "weight")) {
    pv_binding_validate_metadata_string(component$base_name, "base weight name", "weight")
    reps <- pv_binding_schema_names(component$replicate_names, "replicate names", "weight", minimum = 2L)
    count <- pv_binding_schema_count(component$replicate_count, "replicate_count", "weight", minimum = 2L)
    if (length(reps) != count || component$base_name %in% reps ||
        !is.numeric(component$fay_k) || length(component$fay_k) != 1L ||
        is.na(component$fay_k) || !is.finite(component$fay_k) ||
        component$fay_k < 0 || component$fay_k >= 1) {
      pv_binding_schema_abort("weight names/count/Fay metadata is malformed.", "weight")
    }
    for (field in c("base_values_hash", "replicate_values_hash", "membership_values_hash")) {
      pv_binding_schema_hash(component[[field]], field, "weight")
    }
    if (!identical(component$validation_policy_id, "finite_strictly_positive_v1") ||
        !identical(component$target_transform_id, "raw_wls_weights_v1") ||
        !identical(component$stack_transform_id, "mean_one_then_divide_by_M_v1")) {
      pv_binding_schema_abort("weight validation/transform IDs are unsupported.", "weight")
    }
  } else if (identical(component_name, "factor_contrast")) {
    variables <- pv_binding_schema_names(component$variables, "contrast variables", "factor_contrast")
    pv_binding_assert_exact_fields(
      component$option_snapshot,
      c("count", "ordered_label_hash"),
      "factor contrast option snapshot"
    )
    option_count <- pv_binding_schema_count(
      component$option_snapshot$count,
      "option count",
      "factor_contrast",
      minimum = 2L
    )
    if (option_count != 2L) {
      pv_binding_schema_abort("factor contrast option snapshot must contain two ordered option labels.", "factor_contrast")
    }
    pv_binding_schema_hash(component$option_snapshot$ordered_label_hash, "option labels hash", "factor_contrast")
    if (!is.list(component$resolved_specs) || length(component$resolved_specs) != length(variables)) {
      pv_binding_schema_abort("factor contrast specs must align to variables.", "factor_contrast")
    }
    for (spec in component$resolved_specs) {
      pv_binding_validate_contrast_spec_schema(spec, "factor_contrast")
    }
    if (length(variables) > 0L &&
        !identical(vapply(component$resolved_specs, `[[`, character(1), "variable"), variables)) {
      pv_binding_schema_abort("factor contrast variables are not in resolved order.", "factor_contrast")
    }
    pv_binding_schema_hash(component$resolved_hash, "resolved_hash", "factor_contrast")
    expected_resolved_hash <- pv_binding_hash_payload(
      list(
        variables = variables,
        option_snapshot = component$option_snapshot,
        resolved_specs = component$resolved_specs
      ),
      "factor_contrast"
    )
    if (!identical(component$resolved_hash, expected_resolved_hash)) {
      pv_binding_schema_abort("factor contrast resolved_hash is not canonical.", "factor_contrast")
    }
  } else if (identical(component_name, "estimand_contrast")) {
    pv_binding_validate_estimand_contrast_component(component)
  } else if (identical(component_name, "family_link")) {
    rebuilt <- do.call(pv_binding_family_link_projection, component)
    if (!identical(rebuilt, component)) {
      pv_binding_schema_abort("family/link component is not canonical.", "family_link")
    }
  } else if (identical(component_name, "estimand")) {
    pv_binding_validate_string_vector(component$fe_names, "fe_names")
    for (field in c("estimand_id", "target_source", "target_engine_id", "parameter_scope", "interval_role")) {
      pv_binding_validate_metadata_string(component[[field]], field, "estimand")
    }
    if (!is.logical(component$coverage_claim_allowed) ||
        length(component$coverage_claim_allowed) != 1L || is.na(component$coverage_claim_allowed)) {
      pv_binding_schema_abort("estimand coverage flag must be logical.", "estimand")
    }
    pv_binding_schema_hash(component$estimand_contrast_hash, "estimand_contrast_hash", "estimand")
  } else {
    pv_binding_schema_abort("Unknown binding component.", "manifest")
  }
  invisible(component)
}

pv_binding_manifest_build <- function(
  data,
  formula,
  pv_cols,
  weight_col,
  rep_weight_cols,
  fay_k,
  id_cols = NULL,
  family_link,
  estimand_contrast = NULL,
  estimand_metadata,
  model_bundle = NULL
) {
  info <- pv_binding_formula_rhs_info(formula, data = data)
  if (is.null(model_bundle)) {
    model_bundle <- pv_binding_resolve_model_bundle(data, formula)
  }
  pv_binding_validate_model_bundle(model_bundle, data, formula)
  pv_cols <- pv_binding_validate_column_spec(data, pv_cols, "pv", min_length = 1L)
  weight_col <- pv_binding_validate_column_spec(data, weight_col, "weight", min_length = 1L)
  rep_weight_cols <- pv_binding_validate_column_spec(data, rep_weight_cols, "weight", min_length = 2L)
  bound_cols <- unique(c(pv_cols, info$data_variables, weight_col, rep_weight_cols))
  row_identity <- pv_binding_row_identity_digests(
    data,
    id_cols = id_cols,
    bound_cols = bound_cols
  )
  row <- pv_binding_row_projection(data, id_cols = id_cols, bound_cols = bound_cols)
  pv <- pv_binding_pv_projection(
    data,
    pv_cols,
    id_cols = id_cols,
    bound_cols = bound_cols,
    min_length = 1L
  )
  predictor <- model_bundle$predictor_projection
  predictor$membership_values_hash <- pv_binding_keyed_values_hash(
    data,
    predictor$names,
    row_identity$digests,
    "predictor",
    "predictor_membership_values"
  )
  model_matrix <- pv_binding_model_matrix_projection(
    data,
    formula,
    row_identity$digests,
    model_bundle = model_bundle
  )
  formula_component <- pv_binding_formula_component(data, formula, model_matrix)
  weight <- pv_binding_weight_projection(
    data,
    weight_col,
    rep_weight_cols,
    fay_k = fay_k,
    id_cols = id_cols,
    bound_cols = bound_cols
  )
  factor_contrast <- pv_binding_factor_contrast_component(
    data,
    formula,
    model_bundle = model_bundle,
    model_projection = model_matrix
  )

  family_fields <- c("family_id", "link_id", "response_support_id", "dispersion_role")
  pv_binding_assert_exact_fields(family_link, family_fields, "family_link", "PV_BIND_E081")
  family_link <- do.call(pv_binding_family_link_projection, family_link)

  estimand_fields <- c(
    "estimand_id", "target_source", "target_engine_id", "parameter_scope",
    "fe_names", "interval_role", "coverage_claim_allowed"
  )
  pv_binding_assert_exact_fields(estimand_metadata, estimand_fields, "estimand", "PV_BIND_E081")
  estimand_metadata$fe_names <- pv_binding_validate_string_vector(
    estimand_metadata$fe_names,
    "fe_names"
  )
  expected_fe_hash <- pv_binding_ordered_labels_hash(
    estimand_metadata$fe_names,
    "model_matrix",
    "reportable_fe_names"
  )
  if (!identical(expected_fe_hash, model_matrix$ordered_reportable_fe_names_hash)) {
    pv_binding_abort(
      "PV_BIND_E081",
      "Estimand fixed-effect names must exactly match the resolved model-matrix names in order.",
      "estimand"
    )
  }
  if (is.null(estimand_contrast)) {
    estimand_contrast <- pv_binding_estimand_contrast_projection(estimand_metadata$fe_names)
  } else {
    pv_binding_validate_estimand_contrast_component(estimand_contrast)
  }
  estimand <- do.call(
    pv_binding_estimand_projection,
    c(
      list(
        fe_names = estimand_metadata$fe_names,
        estimand_contrast = estimand_contrast,
        interval_role = estimand_metadata$interval_role
      ),
      estimand_metadata[c(
        "estimand_id", "target_source", "target_engine_id", "parameter_scope",
        "coverage_claim_allowed"
      )]
    )
  )

  components <- list(
    row = row,
    pv = pv,
    predictor = predictor,
    formula = formula_component,
    model_matrix = model_matrix,
    weight = weight,
    factor_contrast = factor_contrast,
    estimand_contrast = estimand_contrast,
    family_link = family_link,
    estimand = estimand
  )
  component_hashes <- stats::setNames(
    Map(pv_binding_component_hash, components, names(components)),
    pv_binding_component_hash_names()
  )
  constants <- pv_binding_manifest_constants()
  manifest <- c(
    constants,
    list(
      components = components,
      component_hashes = component_hashes
    )
  )
  manifest$manifest_hash <- pv_binding_hash_payload(
    pv_binding_manifest_hash_payload(manifest),
    "manifest"
  )
  pv_binding_manifest_validate(manifest)
  manifest
}

pv_binding_manifest_validate <- function(manifest) {
  if (is.null(manifest)) {
    pv_binding_abort("PV_BIND_E001", "A binding manifest is required.", "manifest")
  }
  required <- pv_binding_manifest_root_fields()
  optional <- c(
    "legacy_hashes", "created_at", "model_bundle_hash", "migration"
  )
  actual_optional <- if (is.list(manifest) && length(names(manifest)) > length(required)) {
    names(manifest)[-(seq_along(required))]
  } else {
    character()
  }
  if (!is.list(manifest) || anyDuplicated(names(manifest)) ||
      !identical(names(manifest)[seq_along(required)], required) ||
      any(!names(manifest) %in% c(required, optional)) ||
      !identical(actual_optional, optional[optional %in% actual_optional])) {
    pv_binding_abort("PV_BIND_E002", "Manifest root fields and order do not match schema 0.2.0.", "manifest")
  }
  pv_binding_assert_manifest_attributes_recursive(manifest)
  pv_binding_validate_optional_manifest_fields(manifest)
  constants <- pv_binding_manifest_constants()
  if (!identical(manifest$contract_id, constants$contract_id) ||
      !identical(manifest$manifest_schema_version, constants$manifest_schema_version) ||
      !identical(manifest$binding_scope, constants$binding_scope)) {
    pv_binding_abort("PV_BIND_E002", "The binding contract, schema version, or scope is unsupported.", "manifest")
  }
  if (!identical(manifest$canonicalizer_id, constants$canonicalizer_id)) {
    pv_binding_abort("PV_BIND_E003", "The binding canonicalizer is unsupported.", "manifest")
  }
  if (!identical(manifest$hash_algorithm_id, constants$hash_algorithm_id)) {
    pv_binding_abort("PV_BIND_E004", "The binding hash algorithm is unsupported.", "manifest")
  }
  pv_binding_assert_exact_fields(
    manifest$components,
    pv_binding_component_names(),
    "components"
  )
  field_registry <- pv_binding_component_field_registry()
  for (name in pv_binding_component_names()) {
    pv_binding_assert_exact_fields(
      manifest$components[[name]],
      field_registry[[name]],
      paste0(name, " component")
    )
  }
  pv_binding_assert_exact_fields(
    manifest$component_hashes,
    pv_binding_component_hash_names(),
    "component_hashes"
  )
  hashes <- unlist(manifest$component_hashes, use.names = FALSE)
  if (!is.character(hashes) || any(!grepl("^sha256:[0-9a-f]{64}$", hashes))) {
    pv_binding_abort("PV_BIND_E005", "Stored component hashes are malformed.", "manifest")
  }
  expected_manifest_hash <- pv_binding_hash_payload(
    pv_binding_manifest_hash_payload(manifest),
    "manifest"
  )
  if (!is.character(manifest$manifest_hash) || length(manifest$manifest_hash) != 1L ||
      is.na(manifest$manifest_hash) || !identical(manifest$manifest_hash, expected_manifest_hash)) {
    pv_binding_abort(
      "PV_BIND_E005",
      "The manifest self-hash does not match its stored payload.",
      "manifest",
      expected_hash = expected_manifest_hash,
      observed_hash = if (is.character(manifest$manifest_hash) && length(manifest$manifest_hash) == 1L && !is.na(manifest$manifest_hash)) manifest$manifest_hash else "malformed"
    )
  }
  recomputed <- stats::setNames(
    Map(
      pv_binding_component_hash,
      manifest$components,
      names(manifest$components)
    ),
    pv_binding_component_hash_names()
  )
  if (!identical(manifest$component_hashes, recomputed)) {
    pv_binding_abort("PV_BIND_E005", "Stored component hashes do not match current component projections.", "manifest")
  }
  for (name in pv_binding_component_names()) {
    pv_binding_validate_component_deep(manifest$components[[name]], name)
  }
  pv_binding_validate_manifest_attribute_schema(manifest)
  estimand_component <- manifest$components$estimand
  rebuilt_estimand <- pv_binding_estimand_projection(
    fe_names = estimand_component$fe_names,
    estimand_contrast = manifest$components$estimand_contrast,
    interval_role = estimand_component$interval_role,
    estimand_id = estimand_component$estimand_id,
    target_source = estimand_component$target_source,
    target_engine_id = estimand_component$target_engine_id,
    parameter_scope = estimand_component$parameter_scope,
    coverage_claim_allowed = estimand_component$coverage_claim_allowed
  )
  if (!identical(rebuilt_estimand, estimand_component)) {
    pv_binding_abort("PV_BIND_E002", "Estimand component metadata is not canonical.", "estimand")
  }
  if (!identical(manifest$components$row$n, manifest$components$predictor$row_count) ||
      !identical(manifest$components$row$n, manifest$components$model_matrix$row_count) ||
      !identical(manifest$components$pv$names, names(manifest$components$pv$per_column_hashes))) {
    pv_binding_abort("PV_BIND_E002", "Cross-component row/PV dimensions do not align.", "manifest")
  }
  expected_fe_hash <- pv_binding_ordered_labels_hash(
    manifest$components$estimand$fe_names,
    "model_matrix",
    "reportable_fe_names"
  )
  if (!identical(
      expected_fe_hash,
      manifest$components$model_matrix$ordered_reportable_fe_names_hash
    ) || !identical(
      manifest$components$model_matrix$predictor_schema_hash,
      manifest$components$predictor$schema_hash
    ) || !identical(
      manifest$components$model_matrix$predictor_values_hash,
      manifest$components$predictor$ordered_values_hash
    ) || !identical(
      manifest$components$formula$intercept,
      manifest$components$model_matrix$intercept
    ) || !identical(
      manifest$components$formula$ordered_term_labels_hash,
      manifest$components$model_matrix$ordered_term_labels_hash
    ) || !identical(
      manifest$components$formula$offset_hash,
      manifest$components$model_matrix$offset_values_hash
    ) || !identical(
      manifest$components$factor_contrast$resolved_specs,
      manifest$components$model_matrix$contrasts
    )) {
    pv_binding_abort("PV_BIND_E002", "Cross-component formula/design/estimand projections do not align.", "manifest")
  }
  invisible(manifest)
}

pv_binding_manifest_compare <- function(expected, observed) {
  pv_binding_manifest_validate(expected)
  pv_binding_manifest_validate(observed)
  old <- expected$components
  new <- observed$components
  codes <- character()
  add <- function(code, when) {
    if (isTRUE(when)) {
      codes <<- c(codes, code)
    }
  }

  add("PV_BIND_E010", !identical(old$row$n, new$row$n))
  membership_equal <- identical(old$row$row_membership_hash, new$row$row_membership_hash)
  order_equal <- identical(old$row$row_order_hash, new$row$row_order_hash)
  add("PV_BIND_E011", !membership_equal)
  add("PV_BIND_E012", membership_equal && !order_equal)
  pure_reorder <- membership_equal && !order_equal

  add("PV_BIND_E020", !identical(old$pv[c("M", "names", "types")], new$pv[c("M", "names", "types")]))
  add("PV_BIND_E021", !identical(old$pv$membership_values_hash, new$pv$membership_values_hash))
  canonical_predictor_columns <- function(component, include_values) {
    columns <- component$columns
    names(columns) <- vapply(columns, `[[`, character(1), "name")
    columns <- columns[sort(names(columns), method = "radix")]
    if (!include_values) {
      columns <- lapply(columns, function(column) column[setdiff(names(column), "values_hash")])
    } else {
      columns <- lapply(columns, `[[`, "values_hash")
    }
    columns
  }
  predictor_schema_changed <- !identical(
    canonical_predictor_columns(old$predictor, FALSE),
    canonical_predictor_columns(new$predictor, FALSE)
  )
  predictor_values_changed <- if (pure_reorder) {
    !identical(old$predictor$membership_values_hash, new$predictor$membership_values_hash)
  } else {
    !identical(
      canonical_predictor_columns(old$predictor, TRUE),
      canonical_predictor_columns(new$predictor, TRUE)
    )
  }
  add("PV_BIND_E030", predictor_schema_changed)
  add("PV_BIND_E031", predictor_values_changed)

  formula_semantic_fields <- c(
    "outcome_placeholder", "data_variable_count", "data_variables", "rhs_ast_hash",
    "term_count", "ordered_term_labels_hash", "intercept", "specials_count",
    "ordered_specials_hash"
  )
  add("PV_BIND_E040", !identical(old$formula[formula_semantic_fields], new$formula[formula_semantic_fields]))
  model_schema_fields <- c(
    "builder_id", "na_action", "row_count", "column_count", "fe_colname_count",
    "ordered_fe_colnames_hash", "ordered_reportable_fe_names_hash", "assign", "intercept",
    "term_count", "ordered_term_labels_hash", "has_offset", "xlevel_count", "xlevels",
    "contrast_count", "predictor_schema_hash"
  )
  add("PV_BIND_E041", !identical(old$model_matrix[model_schema_fields], new$model_matrix[model_schema_fields]))
  model_ordered_value_changed <- !identical(
    old$model_matrix[c("values_hash", "offset_values_hash")],
    new$model_matrix[c("values_hash", "offset_values_hash")]
  )
  model_membership_changed <- !identical(
    old$model_matrix[c("membership_values_hash", "offset_membership_hash")],
    new$model_matrix[c("membership_values_hash", "offset_membership_hash")]
  )
  add(
    "PV_BIND_E042",
    model_membership_changed || (model_ordered_value_changed && !pure_reorder)
  )
  add("PV_BIND_E043", !identical(old$factor_contrast, new$factor_contrast))

  add("PV_BIND_E050", !identical(
    old$weight[c("base_name", "replicate_names")],
    new$weight[c("base_name", "replicate_names")]
  ))
  add("PV_BIND_E051", !identical(old$weight$membership_values_hash, new$weight$membership_values_hash))
  add("PV_BIND_E052", !identical(
    old$weight[c("validation_policy_id", "target_transform_id", "stack_transform_id")],
    new$weight[c("validation_policy_id", "target_transform_id", "stack_transform_id")]
  ))
  add("PV_BIND_E053", !identical(
    old$weight[c("replicate_count", "fay_k")],
    new$weight[c("replicate_count", "fay_k")]
  ))
  add("PV_BIND_E060", !identical(
    old$family_link[c("family_id", "response_support_id", "dispersion_role")],
    new$family_link[c("family_id", "response_support_id", "dispersion_role")]
  ))
  add("PV_BIND_E061", !identical(old$family_link$link_id, new$family_link$link_id))
  add("PV_BIND_E070", !identical(
    old$estimand[c(
      "estimand_id", "target_source", "target_engine_id", "parameter_scope",
      "interval_role", "coverage_claim_allowed"
    )],
    new$estimand[c(
      "estimand_id", "target_source", "target_engine_id", "parameter_scope",
      "interval_role", "coverage_claim_allowed"
    )]
  ))
  add("PV_BIND_E071", !identical(old$estimand$fe_names, new$estimand$fe_names))
  add("PV_BIND_E072", !identical(old$estimand_contrast, new$estimand_contrast) ||
      !identical(old$estimand$estimand_contrast_hash, new$estimand$estimand_contrast_hash))

  mismatched_components <- pv_binding_component_names()[
    !vapply(
      pv_binding_component_hash_names(),
      function(name) identical(expected$component_hashes[[name]], observed$component_hashes[[name]]),
      logical(1)
    )
  ]
  component_code_sets <- list(
    row = c("PV_BIND_E010", "PV_BIND_E011", "PV_BIND_E012"),
    pv = c("PV_BIND_E020", "PV_BIND_E021"),
    predictor = c("PV_BIND_E030", "PV_BIND_E031"),
    formula = "PV_BIND_E040",
    model_matrix = c("PV_BIND_E041", "PV_BIND_E042"),
    weight = c("PV_BIND_E050", "PV_BIND_E051", "PV_BIND_E052", "PV_BIND_E053"),
    factor_contrast = "PV_BIND_E043",
    estimand_contrast = "PV_BIND_E072",
    family_link = c("PV_BIND_E060", "PV_BIND_E061"),
    estimand = c("PV_BIND_E070", "PV_BIND_E071", "PV_BIND_E072")
  )
  fallback_codes <- c(
    row = "PV_BIND_E011",
    pv = "PV_BIND_E021",
    predictor = "PV_BIND_E031",
    formula = "PV_BIND_E040",
    model_matrix = "PV_BIND_E042",
    weight = "PV_BIND_E051",
    factor_contrast = "PV_BIND_E043",
    estimand_contrast = "PV_BIND_E072",
    family_link = "PV_BIND_E060",
    estimand = "PV_BIND_E070"
  )
  for (component in mismatched_components) {
    represented_by_row_order <- pure_reorder && component %in% c(
      "pv", "predictor", "model_matrix", "weight"
    ) && if (component == "pv") {
      identical(old$pv$membership_values_hash, new$pv$membership_values_hash)
    } else if (component == "predictor") {
      identical(old$predictor$membership_values_hash, new$predictor$membership_values_hash)
    } else if (component == "model_matrix") {
      identical(
        old$model_matrix[c("membership_values_hash", "offset_membership_hash")],
        new$model_matrix[c("membership_values_hash", "offset_membership_hash")]
      )
    } else {
      identical(old$weight$membership_values_hash, new$weight$membership_values_hash)
    }
    represented_by_formula <- component == "predictor" &&
      "PV_BIND_E040" %in% codes && !predictor_schema_changed && !predictor_values_changed
    represented_by_model_value <- component == "formula" &&
      "PV_BIND_E042" %in% codes &&
      identical(old$formula[formula_semantic_fields], new$formula[formula_semantic_fields])
    if (length(intersect(codes, component_code_sets[[component]])) == 0L) {
      if (!represented_by_row_order && !represented_by_formula && !represented_by_model_value) {
        codes <- c(codes, unname(fallback_codes[[component]]))
      }
    }
  }
  registry_order <- names(pv_binding_error_registry())
  all_codes <- registry_order[registry_order %in% unique(codes)]
  list(
    ok = length(mismatched_components) == 0L && length(all_codes) == 0L,
    all_codes = all_codes,
    primary_code = if (length(all_codes) == 0L) NA_character_ else all_codes[[1L]],
    mismatched_components = mismatched_components
  )
}

pv_binding_proof_fields <- function() {
  c(
    "contract_id", "manifest_schema_version", "target_manifest_hash",
    "current_manifest_hash", "verification_policy"
  )
}

pv_binding_proof_policy <- function() {
  "recomputed_current_inputs_exact_component_match_v1"
}

pv_binding_manifest_assert_match <- function(target_manifest, current_manifest) {
  comparison <- pv_binding_manifest_compare(target_manifest, current_manifest)
  if (!isTRUE(comparison$ok)) {
    if (!is.character(comparison$primary_code) || length(comparison$primary_code) != 1L ||
        is.na(comparison$primary_code) || !comparison$primary_code %in% names(pv_binding_error_registry())) {
      pv_binding_abort(
        "PV_BIND_E005",
        "Binding comparison failed without a canonical primary mismatch code.",
        "manifest"
      )
    }
    safe_components <- paste(comparison$mismatched_components, collapse = ", ")
    pv_binding_abort(
      comparison$primary_code,
      sprintf(
        "Target and current binding manifests differ in canonical component(s): %s.",
        safe_components
      ),
      component = "manifest",
      expected_hash = target_manifest$manifest_hash,
      observed_hash = current_manifest$manifest_hash,
      all_codes = comparison$all_codes
    )
  }
  invisible(comparison)
}

pv_binding_proof_validate <- function(proof, target_manifest = NULL) {
  pv_binding_assert_exact_fields(
    proof,
    pv_binding_proof_fields(),
    "binding proof",
    "PV_BIND_E002"
  )
  if (any(!vapply(proof, function(value) is.null(attributes(value)), logical(1)))) {
    pv_binding_abort(
      "PV_BIND_E002",
      "Binding proof scalar fields must be attribute-free.",
      "binding_proof"
    )
  }
  constants <- pv_binding_manifest_constants()
  if (!identical(proof$contract_id, constants$contract_id) ||
      !identical(proof$manifest_schema_version, constants$manifest_schema_version) ||
      !identical(proof$verification_policy, pv_binding_proof_policy())) {
    pv_binding_abort(
      "PV_BIND_E002",
      "Binding proof contract, schema, or verification policy is unsupported.",
      "binding_proof"
    )
  }
  hashes <- c(proof$target_manifest_hash, proof$current_manifest_hash)
  if (!is.character(hashes) || length(hashes) != 2L || anyNA(hashes) ||
      any(!grepl("^sha256:[0-9a-f]{64}$", hashes)) ||
      !identical(proof$target_manifest_hash, proof$current_manifest_hash)) {
    pv_binding_abort(
      "PV_BIND_E005",
      "Binding proof hashes must be identical canonical SHA-256 manifest hashes.",
      "binding_proof",
      expected_hash = proof$target_manifest_hash,
      observed_hash = proof$current_manifest_hash
    )
  }
  if (!is.null(target_manifest)) {
    pv_binding_manifest_assert_reportable(target_manifest)
    if (!identical(proof$target_manifest_hash, target_manifest$manifest_hash)) {
      pv_binding_abort(
        "PV_BIND_E005",
        "Binding proof does not link to the validated target manifest.",
        "binding_proof",
        expected_hash = target_manifest$manifest_hash,
        observed_hash = proof$target_manifest_hash
      )
    }
  }
  invisible(proof)
}

pv_binding_proof_build <- function(target_manifest, current_manifest) {
  pv_binding_manifest_assert_reportable(target_manifest)
  pv_binding_manifest_assert_match(target_manifest, current_manifest)
  constants <- pv_binding_manifest_constants()
  proof <- list(
    contract_id = constants$contract_id,
    manifest_schema_version = constants$manifest_schema_version,
    target_manifest_hash = target_manifest$manifest_hash,
    current_manifest_hash = current_manifest$manifest_hash,
    verification_policy = pv_binding_proof_policy()
  )
  pv_binding_proof_validate(proof, target_manifest = target_manifest)
  proof
}

pv_binding_retained_raw_inputs_validate <- function(data, formula, target) {
  if (!is.data.frame(data) || !inherits(formula, "formula")) {
    pv_binding_abort(
      "PV_BIND_E081",
      "Retained binding attestation requires a data frame and formula.",
      "retained_design"
    )
  }
  validate_pvstackr_brr_target(target)
  manifest <- target$binding_manifest
  pv_binding_manifest_assert_reportable(manifest)
  expected <- manifest$components
  info <- pv_binding_formula_rhs_info(formula, data = data)
  bound_cols <- unique(c(
    target$pv_cols,
    info$data_variables,
    target$weight_col,
    target$rep_weight_cols
  ))
  row_identity <- pv_binding_row_identity_digests(
    data,
    id_cols = target$id_cols,
    bound_cols = bound_cols
  )
  observed <- list(
    row = pv_binding_row_projection(
      data,
      id_cols = target$id_cols,
      bound_cols = bound_cols
    ),
    pv = pv_binding_pv_projection(
      data,
      target$pv_cols,
      id_cols = target$id_cols,
      bound_cols = bound_cols,
      min_length = 1L
    ),
    predictor = pv_binding_predictor_projection(data, formula),
    weight = pv_binding_weight_projection(
      data,
      target$weight_col,
      target$rep_weight_cols,
      fay_k = target$fay_k,
      id_cols = target$id_cols,
      bound_cols = bound_cols
    )
  )
  observed$predictor$membership_values_hash <- pv_binding_keyed_values_hash(
    data,
    observed$predictor$names,
    row_identity$digests,
    "predictor",
    "predictor_membership_values"
  )

  codes <- character()
  add <- function(code, condition) {
    if (isTRUE(condition)) {
      codes <<- c(codes, code)
    }
  }
  add("PV_BIND_E010", !identical(expected$row$n, observed$row$n))
  membership_equal <- identical(
    expected$row$row_membership_hash,
    observed$row$row_membership_hash
  )
  order_equal <- identical(
    expected$row$row_order_hash,
    observed$row$row_order_hash
  )
  add("PV_BIND_E011", !membership_equal)
  add("PV_BIND_E012", membership_equal && !order_equal)
  pure_reorder <- membership_equal && !order_equal

  add("PV_BIND_E020", !identical(
    expected$pv[c("M", "names", "types")],
    observed$pv[c("M", "names", "types")]
  ))
  pv_values_equal <- if (pure_reorder) {
    identical(expected$pv$membership_values_hash, observed$pv$membership_values_hash)
  } else {
    identical(expected$pv$ordered_values_hash, observed$pv$ordered_values_hash)
  }
  add("PV_BIND_E021", !pv_values_equal)

  predictor_schema_fields <- c(
    "source_rule", "row_count", "column_count", "names", "schema_hash"
  )
  add("PV_BIND_E030", !identical(
    expected$predictor[predictor_schema_fields],
    observed$predictor[predictor_schema_fields]
  ))
  predictor_values_equal <- if (pure_reorder) {
    identical(
      expected$predictor$membership_values_hash,
      observed$predictor$membership_values_hash
    )
  } else {
    identical(
      expected$predictor$ordered_values_hash,
      observed$predictor$ordered_values_hash
    )
  }
  add("PV_BIND_E031", !predictor_values_equal)

  observed_ast <- pv_binding_formula_ast(formula, data = data)
  add("PV_BIND_E040", !identical(
    list(
      data_variable_count = expected$formula$data_variable_count,
      data_variables = expected$formula$data_variables,
      rhs_ast_hash = expected$formula$rhs_ast_hash
    ),
    list(
      data_variable_count = length(observed_ast$data_variables),
      data_variables = observed_ast$data_variables,
      rhs_ast_hash = pv_binding_hash_payload(observed_ast$rhs, "formula")
    )
  ))

  add("PV_BIND_E050", !identical(
    expected$weight[c("base_name", "replicate_names")],
    observed$weight[c("base_name", "replicate_names")]
  ))
  add("PV_BIND_E051", !identical(
    expected$weight$membership_values_hash,
    observed$weight$membership_values_hash
  ))
  add("PV_BIND_E052", !identical(
    expected$weight[c(
      "validation_policy_id", "target_transform_id", "stack_transform_id"
    )],
    observed$weight[c(
      "validation_policy_id", "target_transform_id", "stack_transform_id"
    )]
  ))
  add("PV_BIND_E053", !identical(
    expected$weight[c("replicate_count", "fay_k")],
    observed$weight[c("replicate_count", "fay_k")]
  ))

  registry <- names(pv_binding_error_registry())
  codes <- registry[registry %in% unique(codes)]
  if (length(codes) > 0L) {
    pv_binding_abort(
      codes[[1L]],
      "Retained raw design inputs no longer match the authenticated target binding.",
      "retained_design",
      all_codes = codes
    )
  }
  invisible(TRUE)
}

pv_binding_target_content_fields <- function() {
  c(
    "content_schema_version", "manifest_hash", "target_source", "target_engine_id",
    "M", "R", "fay_k", "fe_names", "per_pv", "interval_policy",
    "target_policy", "numeric_comparison_policy", "derived", "target_content_hash"
  )
}

pv_binding_target_per_pv_fields <- function() {
  c("pv_index", "pv_col", "beta", "U")
}

pv_binding_target_interval_fields <- function() {
  c("df_method", "df_complete", "conf_level", "interval_role", "coverage_claim_allowed")
}

pv_binding_target_policy_fields <- function() {
  c(
    "replicate_weight_role", "target_repair", "fixed_effects_only", "df_method",
    "interval_role", "coverage_claim_allowed"
  )
}

pv_binding_target_derived_fields <- function() {
  c("beta", "U_bar", "B", "T_MI", "se", "df", "df_classic")
}

pv_binding_target_numeric_policy <- function() {
  list(
    policy_id = "absolute_relative_v1",
    absolute_tolerance = 1e-12,
    relative_tolerance = 1e-10,
    infinity_policy = "positive_only_for_df"
  )
}

pv_binding_target_matrix_tolerance <- function() {
  1e-10
}

pv_binding_target_abort <- function(detail, component = "target_content") {
  pv_binding_abort("PV_BIND_E090", detail, component)
}

pv_binding_target_validate_string <- function(x, label) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x) ||
      !validUTF8(enc2utf8(x)) || trimws(x) != x) {
    pv_binding_target_abort(sprintf("%s must be a trimmed non-empty UTF-8 string.", label))
  }
  x
}

pv_binding_target_validate_scalar_integer <- function(x, label, lower) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != floor(x) || x < lower || x > .Machine$integer.max) {
    pv_binding_target_abort(sprintf("%s must be a finite integer at least %d.", label, lower))
  }
  as.integer(x)
}

pv_binding_target_validate_fay <- function(fay_k) {
  if (!is.numeric(fay_k) || length(fay_k) != 1L || is.na(fay_k) ||
      !is.finite(fay_k) || fay_k < 0 || fay_k >= 1) {
    pv_binding_target_abort("fay_k must be finite in [0, 1).")
  }
  as.double(fay_k)
}

pv_binding_target_validate_raw_vector <- function(x, fe_names, label) {
  if (!is.numeric(x) || !identical(names(x), fe_names) || any(!is.finite(x))) {
    pv_binding_target_abort(sprintf("%s must be a finite named vector aligned to fixed effects.", label))
  }
  as.double(x) |> stats::setNames(fe_names)
}

pv_binding_target_validate_raw_matrix <- function(x, fe_names, label) {
  p <- length(fe_names)
  if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x)) ||
      !identical(dim(x), c(p, p)) || !identical(rownames(x), fe_names) ||
      !identical(colnames(x), fe_names)) {
    pv_binding_target_abort(sprintf("%s must be a finite named square matrix aligned to fixed effects.", label))
  }
  if (!isTRUE(all.equal(x, t(x), tolerance = 1e-12, check.attributes = FALSE))) {
    pv_binding_target_abort(sprintf("%s must be symmetric.", label))
  }
  eigenvalues <- tryCatch(
    eigen(0.5 * (x + t(x)), symmetric = TRUE, only.values = TRUE)$values,
    error = function(error) pv_binding_target_abort(sprintf("%s eigendecomposition failed.", label))
  )
  if (min(eigenvalues) < -1e-10) {
    pv_binding_target_abort(sprintf("%s must be positive semidefinite.", label))
  }
  out <- 0.5 * (x + t(x))
  dimnames(out) <- list(fe_names, fe_names)
  out
}

pv_binding_target_normalize_df_complete <- function(df_complete, fe_names, df_method) {
  if (!is.character(df_method) || length(df_method) != 1L || is.na(df_method) ||
      !df_method %in% c("classic", "barnard_rubin")) {
    pv_binding_target_abort("df_method must be classic or barnard_rubin.")
  }
  if (identical(df_method, "classic")) {
    if (!is.null(df_complete)) {
      pv_binding_target_abort("Classic Rubin intervals require df_complete = NULL.")
    }
    return(NULL)
  }
  if (!is.numeric(df_complete) || length(df_complete) < 1L || any(is.na(df_complete)) ||
      any(is.nan(df_complete)) || any(df_complete <= 0) || any(df_complete == -Inf)) {
    pv_binding_target_abort("Barnard-Rubin df_complete must contain positive finite values or +Inf.")
  }
  if (length(df_complete) == 1L) {
    return(stats::setNames(rep(as.double(df_complete), length(fe_names)), fe_names))
  }
  if (!identical(length(df_complete), length(fe_names)) || is.null(names(df_complete)) ||
      anyDuplicated(names(df_complete)) || !setequal(names(df_complete), fe_names)) {
    pv_binding_target_abort("Named df_complete must align exactly to fixed effects.")
  }
  as.double(df_complete[fe_names]) |> stats::setNames(fe_names)
}

pv_binding_target_interval_policy <- function(
  df_method,
  df_complete,
  conf_level,
  interval_role,
  coverage_claim_allowed,
  fe_names
) {
  df_complete <- pv_binding_target_normalize_df_complete(df_complete, fe_names, df_method)
  if (!is.numeric(conf_level) || length(conf_level) != 1L || is.na(conf_level) ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    pv_binding_target_abort("conf_level must be finite in (0, 1).")
  }
  interval_role <- pv_binding_target_validate_string(interval_role, "interval_role")
  if (!is.logical(coverage_claim_allowed) || length(coverage_claim_allowed) != 1L ||
      is.na(coverage_claim_allowed)) {
    pv_binding_target_abort("coverage_claim_allowed must be a logical scalar.")
  }
  expected <- if (identical(df_method, "barnard_rubin")) {
    list(interval_role = "coverage_barnard_rubin", coverage_claim_allowed = TRUE)
  } else {
    list(interval_role = "descriptive_classic_rubin", coverage_claim_allowed = FALSE)
  }
  if (!identical(interval_role, expected$interval_role) ||
      !identical(coverage_claim_allowed, expected$coverage_claim_allowed)) {
    pv_binding_target_abort("Interval role and coverage policy are inconsistent with df_method.")
  }
  list(
    df_method = df_method,
    df_complete = df_complete,
    conf_level = as.double(conf_level),
    interval_role = interval_role,
    coverage_claim_allowed = coverage_claim_allowed
  )
}

pv_binding_target_policy_projection <- function(
  df_method,
  interval_role,
  coverage_claim_allowed,
  replicate_weight_role = "external_design_variance_only",
  target_repair = "forbidden",
  fixed_effects_only = TRUE
) {
  replicate_weight_role <- pv_binding_target_validate_string(
    replicate_weight_role,
    "replicate_weight_role"
  )
  target_repair <- pv_binding_target_validate_string(target_repair, "target_repair")
  if (!identical(replicate_weight_role, "external_design_variance_only") ||
      !identical(target_repair, "forbidden") || !identical(fixed_effects_only, TRUE)) {
    pv_binding_target_abort("Target policy must preserve external-design, no-repair, fixed-effect semantics.")
  }
  list(
    replicate_weight_role = replicate_weight_role,
    target_repair = target_repair,
    fixed_effects_only = fixed_effects_only,
    df_method = df_method,
    interval_role = interval_role,
    coverage_claim_allowed = coverage_claim_allowed
  )
}

pv_binding_target_numeric_equal <- function(expected, observed, allow_positive_infinity = FALSE) {
  if (!is.numeric(expected) || !is.numeric(observed) || typeof(expected) != typeof(observed) ||
      !identical(attributes(expected), attributes(observed)) || !identical(dim(expected), dim(observed))) {
    return(FALSE)
  }
  if (length(expected) != length(observed) || any(is.na(expected)) || any(is.na(observed)) ||
      any(is.nan(expected)) || any(is.nan(observed))) {
    return(FALSE)
  }
  expected_infinite <- is.infinite(expected)
  observed_infinite <- is.infinite(observed)
  if (!identical(expected_infinite, observed_infinite)) {
    return(FALSE)
  }
  if (any(expected_infinite)) {
    if (!allow_positive_infinity || any(expected[expected_infinite] < 0) ||
        any(observed[observed_infinite] < 0) ||
        !identical(expected[expected_infinite], observed[observed_infinite])) {
      return(FALSE)
    }
  }
  finite <- !expected_infinite
  policy <- pv_binding_target_numeric_policy()
  difference <- abs(expected[finite] - observed[finite])
  scale <- pmax(abs(expected[finite]), abs(observed[finite]))
  all(difference <= policy$absolute_tolerance + policy$relative_tolerance * scale)
}

pv_binding_target_derived_validate <- function(derived) {
  pv_binding_assert_exact_fields(
    derived,
    pv_binding_target_derived_fields(),
    "target derived",
    "PV_BIND_E090"
  )
  fe_names <- names(derived$beta)
  if (!is.character(fe_names) || length(fe_names) < 1L || anyNA(fe_names) ||
      any(!nzchar(fe_names)) || anyDuplicated(fe_names) || any(!is.finite(derived$beta))) {
    pv_binding_target_abort("Derived beta must be finite, named, and unique.")
  }
  for (field in c("U_bar", "B")) {
    pv_binding_target_validate_raw_matrix(derived[[field]], fe_names, paste0("derived ", field))
  }
  T_MI <- pv_binding_target_validate_raw_matrix(derived$T_MI, fe_names, "derived T_MI")
  eigenvalues <- eigen(T_MI, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigenvalues <= 0)) {
    pv_binding_target_abort("Derived T_MI must be strictly positive definite.")
  }
  if (!is.numeric(derived$se) || !identical(names(derived$se), fe_names) ||
      any(!is.finite(derived$se)) || any(derived$se <= 0)) {
    pv_binding_target_abort("Derived standard errors must be finite, strictly positive, and aligned.")
  }
  for (field in c("df", "df_classic")) {
    value <- derived[[field]]
    if (!is.numeric(value) || !identical(names(value), fe_names) || any(is.na(value)) ||
        any(is.nan(value)) || any(value <= 0) || any(value == -Inf)) {
      pv_binding_target_abort("Derived degrees of freedom must be positive finite values or +Inf.")
    }
  }
  invisible(derived)
}

pv_binding_target_compare_derived <- function(expected, observed) {
  pv_binding_target_derived_validate(expected)
  pv_binding_target_derived_validate(observed)
  fields <- pv_binding_target_derived_fields()
  equal <- vapply(fields, function(field) {
    pv_binding_target_numeric_equal(
      expected[[field]],
      observed[[field]],
      allow_positive_infinity = field %in% c("df", "df_classic")
    )
  }, logical(1))
  if (!all(equal)) {
    pv_binding_target_abort("Stored derived target fields do not match canonical Rubin recomputation within policy tolerance.")
  }
  invisible(TRUE)
}

pv_binding_target_content_compute <- function(
  per_pv,
  M,
  R,
  fay_k,
  df_method,
  df_complete,
  conf_level,
  target_source,
  target_engine_id,
  target_policy,
  manifest_hash
) {
  M <- pv_binding_target_validate_scalar_integer(M, "M", 1L)
  R <- pv_binding_target_validate_scalar_integer(R, "R", 2L)
  fay_k <- pv_binding_target_validate_fay(fay_k)
  if (!is.list(per_pv) || length(per_pv) != M) {
    pv_binding_target_abort("per_pv length must exactly equal M.")
  }
  target_source <- pv_binding_target_validate_string(target_source, "target_source")
  target_engine_id <- pv_binding_target_validate_string(target_engine_id, "target_engine_id")
  if (!identical(target_source, "external_brr_fay_rubin") ||
      !identical(target_engine_id, "lm_wls_brr_fay_v1")) {
    pv_binding_target_abort("Target source and engine IDs are unsupported for BRR-Fay content.")
  }
  if (!is.character(manifest_hash) || length(manifest_hash) != 1L || is.na(manifest_hash) ||
      !grepl("^sha256:[0-9a-f]{64}$", manifest_hash)) {
    pv_binding_target_abort("manifest_hash must be a canonical SHA-256 digest.")
  }

  first_beta <- per_pv[[1L]]$beta
  fe_names <- names(first_beta)
  if (!is.character(fe_names) || length(fe_names) < 1L || anyNA(fe_names) ||
      any(!nzchar(fe_names)) || anyDuplicated(fe_names)) {
    pv_binding_target_abort("The first per-PV beta must define unique fixed-effect names.")
  }
  canonical_per_pv <- lapply(seq_along(per_pv), function(index) {
    item <- per_pv[[index]]
    if (!is.list(item) || !all(c("beta", "U") %in% names(item))) {
      pv_binding_target_abort("Each per-PV primitive must contain beta and U.")
    }
    if (!is.null(item$fe_names) && !identical(item$fe_names, fe_names)) {
      pv_binding_target_abort("Per-PV fixed-effect names must match in exact order.")
    }
    if (!is.null(item$R)) {
      item_R <- pv_binding_target_validate_scalar_integer(item$R, "per-PV R", 2L)
      if (!identical(item_R, R)) {
        pv_binding_target_abort("Per-PV replicate count must match R.")
      }
    }
    if (!is.null(item$fay_k)) {
      item_fay <- pv_binding_target_validate_fay(item$fay_k)
      if (!identical(item_fay, fay_k)) {
        pv_binding_target_abort("Per-PV Fay coefficient must match fay_k.")
      }
    }
    pv_col <- item$pv_col
    if (is.null(pv_col)) {
      pv_col <- paste0("PV", index)
    }
    pv_col <- pv_binding_target_validate_string(pv_col, "pv_col")
    list(
      pv_index = as.integer(index),
      pv_col = pv_col,
      beta = pv_binding_target_validate_raw_vector(item$beta, fe_names, "per-PV beta"),
      U = pv_binding_target_validate_raw_matrix(item$U, fe_names, "per-PV U")
    )
  })
  pv_labels <- vapply(canonical_per_pv, `[[`, character(1), "pv_col")
  if (anyDuplicated(pv_labels)) {
    pv_binding_target_abort("Per-PV labels must be unique and ordered.")
  }
  beta_rows <- do.call(rbind, lapply(canonical_per_pv, `[[`, "beta"))
  colnames(beta_rows) <- fe_names
  pv_binding_assert_exact_fields(
    target_policy,
    pv_binding_target_policy_fields(),
    "target policy",
    "PV_BIND_E090"
  )
  interval_role <- target_policy$interval_role
  coverage_claim_allowed <- target_policy$coverage_claim_allowed
  interval_policy <- pv_binding_target_interval_policy(
    df_method = df_method,
    df_complete = df_complete,
    conf_level = conf_level,
    interval_role = interval_role,
    coverage_claim_allowed = coverage_claim_allowed,
    fe_names = fe_names
  )
  canonical_policy <- pv_binding_target_policy_projection(
    df_method = interval_policy$df_method,
    interval_role = interval_policy$interval_role,
    coverage_claim_allowed = interval_policy$coverage_claim_allowed,
    replicate_weight_role = target_policy$replicate_weight_role,
    target_repair = target_policy$target_repair,
    fixed_effects_only = target_policy$fixed_effects_only
  )
  if (!identical(target_policy, canonical_policy)) {
    pv_binding_target_abort("Target policy fields are inconsistent with interval primitives.")
  }

  pooled <- tryCatch(
    rubin_pool_matrix(
      beta = beta_rows,
      U = lapply(canonical_per_pv, `[[`, "U"),
      orientation = "rows_pv",
      conf_level = interval_policy$conf_level,
      allow_m1 = M == 1L,
      df_method = interval_policy$df_method,
      df_complete = interval_policy$df_complete
    ),
    error = function(error) {
      if (inherits(error, "pvstackr_binding_error")) {
        stop(error)
      }
      pv_binding_target_abort("Canonical Rubin recomputation failed.")
    }
  )
  derived <- list(
    beta = pooled$beta,
    U_bar = pooled$U_bar,
    B = pooled$B,
    T_MI = pooled$T_MI,
    se = pooled$se,
    df = pooled$df,
    df_classic = pooled$df_classic
  )
  pv_binding_target_derived_validate(derived)
  list(
    content_schema_version = "pvstackr_target_content_v1",
    manifest_hash = manifest_hash,
    target_source = target_source,
    target_engine_id = target_engine_id,
    M = M,
    R = R,
    fay_k = fay_k,
    fe_names = fe_names,
    per_pv = canonical_per_pv,
    interval_policy = interval_policy,
    target_policy = canonical_policy,
    numeric_comparison_policy = pv_binding_target_numeric_policy(),
    derived = derived
  )
}

pv_binding_target_content_hash_payload <- function(target_content) {
  target_content[setdiff(pv_binding_target_content_fields(), "target_content_hash")]
}

pv_binding_target_content_build <- function(
  per_pv,
  M,
  R,
  fay_k,
  df_method,
  df_complete = NULL,
  conf_level,
  target_source,
  target_engine_id,
  target_policy,
  manifest_hash,
  stored_derived = NULL
) {
  content <- pv_binding_target_content_compute(
    per_pv = per_pv,
    M = M,
    R = R,
    fay_k = fay_k,
    df_method = df_method,
    df_complete = df_complete,
    conf_level = conf_level,
    target_source = target_source,
    target_engine_id = target_engine_id,
    target_policy = target_policy,
    manifest_hash = manifest_hash
  )
  if (!is.null(stored_derived)) {
    pv_binding_target_compare_derived(content$derived, stored_derived)
  }
  content$target_content_hash <- pv_binding_hash_payload(
    pv_binding_target_content_hash_payload(content),
    "target_content"
  )
  pv_binding_target_content_validate(content)
  content
}

pv_binding_target_content_validate <- function(target_content) {
  pv_binding_assert_exact_fields(
    target_content,
    pv_binding_target_content_fields(),
    "target_content",
    "PV_BIND_E090"
  )
  pv_binding_assert_exact_fields(
    target_content$interval_policy,
    pv_binding_target_interval_fields(),
    "target interval policy",
    "PV_BIND_E090"
  )
  pv_binding_assert_exact_fields(
    target_content$target_policy,
    pv_binding_target_policy_fields(),
    "target policy",
    "PV_BIND_E090"
  )
  pv_binding_assert_exact_fields(
    target_content$numeric_comparison_policy,
    names(pv_binding_target_numeric_policy()),
    "target numeric policy",
    "PV_BIND_E090"
  )
  pv_binding_target_derived_validate(target_content$derived)
  if (!is.list(target_content$per_pv) || length(target_content$per_pv) != target_content$M) {
    pv_binding_target_abort("Stored per_pv primitives do not align with M.")
  }
  for (item in target_content$per_pv) {
    pv_binding_assert_exact_fields(
      item,
      pv_binding_target_per_pv_fields(),
      "target per_pv primitive",
      "PV_BIND_E090"
    )
  }
  expected_hash <- pv_binding_hash_payload(
    pv_binding_target_content_hash_payload(target_content),
    "target_content"
  )
  if (!is.character(target_content$target_content_hash) ||
      length(target_content$target_content_hash) != 1L ||
      is.na(target_content$target_content_hash) ||
      !identical(target_content$target_content_hash, expected_hash)) {
    pv_binding_target_abort("The target content self-hash is stale or malformed.")
  }
  recomputed <- pv_binding_target_content_compute(
    per_pv = target_content$per_pv,
    M = target_content$M,
    R = target_content$R,
    fay_k = target_content$fay_k,
    df_method = target_content$interval_policy$df_method,
    df_complete = target_content$interval_policy$df_complete,
    conf_level = target_content$interval_policy$conf_level,
    target_source = target_content$target_source,
    target_engine_id = target_content$target_engine_id,
    target_policy = target_content$target_policy,
    manifest_hash = target_content$manifest_hash
  )
  if (!identical(target_content$numeric_comparison_policy, recomputed$numeric_comparison_policy)) {
    pv_binding_target_abort("The target numeric comparison policy is not canonical.")
  }
  pv_binding_target_compare_derived(recomputed$derived, target_content$derived)
  canonical_fields <- setdiff(names(recomputed), "derived")
  if (!identical(target_content[canonical_fields], recomputed[canonical_fields])) {
    pv_binding_target_abort("Stored target primitives or policy are not canonical.")
  }
  invisible(target_content)
}

pv_binding_target_manifest_validate <- function(target_content, manifest) {
  pv_binding_target_content_validate(target_content)
  pv_binding_manifest_validate(manifest)
  components <- manifest$components

  if (!identical(target_content$manifest_hash, manifest$manifest_hash)) {
    pv_binding_abort(
      "PV_BIND_E090",
      "Target content is not linked to this validated binding manifest.",
      "target_content"
    )
  }
  target_pv_names <- vapply(target_content$per_pv, `[[`, character(1), "pv_col")
  if (!identical(target_content$M, components$pv$M) ||
      !identical(target_pv_names, components$pv$names)) {
    pv_binding_abort(
      "PV_BIND_E020",
      "Target plausible-value count, labels, or order do not match the binding manifest.",
      "target_content"
    )
  }
  if (!identical(target_content$R, components$weight$replicate_count) ||
      !identical(target_content$fay_k, components$weight$fay_k)) {
    pv_binding_abort(
      "PV_BIND_E053",
      "Target replicate count or Fay coefficient does not match the binding manifest.",
      "target_content"
    )
  }
  if (!identical(target_content$fe_names, components$estimand$fe_names)) {
    pv_binding_abort(
      "PV_BIND_E071",
      "Target fixed-effect names do not match the manifest estimand in exact order.",
      "target_content"
    )
  }
  if (!identical(target_content$target_source, components$estimand$target_source) ||
      !identical(target_content$target_engine_id, components$estimand$target_engine_id) ||
      !identical(target_content$interval_policy$interval_role, components$estimand$interval_role) ||
      !identical(
        target_content$interval_policy$coverage_claim_allowed,
        components$estimand$coverage_claim_allowed
      ) || !identical(components$estimand$estimand_id, "brr_fay_fixed_effects") ||
      !identical(components$estimand$parameter_scope, "fixed_effect")) {
    pv_binding_abort(
      "PV_BIND_E070",
      "Target source, engine, interval, coverage, or estimand semantics do not match the manifest.",
      "target_content"
    )
  }
  expected_family_link <- pv_binding_family_link_projection(
    family_id = "gaussian",
    link_id = "identity",
    response_support_id = "real",
    dispersion_role = "estimated"
  )
  if (!identical(
    components$family_link[c(
      "family_id", "response_support_id", "dispersion_role"
    )],
    expected_family_link[c(
      "family_id", "response_support_id", "dispersion_role"
    )]
  )) {
    pv_binding_abort("PV_BIND_E060", "BRR-Fay target content requires canonical Gaussian real-support estimated-dispersion semantics.", "target_content")
  }
  if (!identical(components$family_link$link_id, expected_family_link$link_id)) {
    pv_binding_abort("PV_BIND_E061", "BRR-Fay target content requires identity-link manifest semantics.", "target_content")
  }
  if (!identical(target_content$target_policy$replicate_weight_role, "external_design_variance_only") ||
      !identical(target_content$target_policy$target_repair, "forbidden") ||
      !identical(target_content$target_policy$fixed_effects_only, TRUE) ||
      !identical(components$weight$target_transform_id, "raw_wls_weights_v1") ||
      !identical(components$estimand_contrast$contrast_id, "identity_fixed_effect_v1") ||
      !identical(
        components$estimand_contrast$input_fe_count,
        as.integer(length(target_content$fe_names))
      ) || !identical(
        components$estimand_contrast$output_count,
        as.integer(length(target_content$fe_names))
      )) {
    pv_binding_abort(
      "PV_BIND_E070",
      "Target policy or fixed-effect contrast semantics do not match the manifest.",
      "target_content"
    )
  }
  expected_fe_hash <- pv_binding_ordered_labels_hash(
    target_content$fe_names,
    "estimand_contrast",
    "input_fe_names"
  )
  if (!identical(
    expected_fe_hash,
    components$estimand_contrast$ordered_input_fe_names_hash
  )) {
    pv_binding_abort(
      "PV_BIND_E071",
      "Target fixed-effect hash does not match the manifest contrast basis.",
      "target_content"
    )
  }
  expected_contrast <- pv_binding_estimand_contrast_projection(
    target_content$fe_names
  )
  expected_contrast_hash <- pv_binding_hash_payload(
    expected_contrast,
    "estimand_contrast"
  )
  if (!identical(components$estimand_contrast, expected_contrast) ||
      !identical(
        components$estimand$estimand_contrast_hash,
        expected_contrast_hash
      )) {
    pv_binding_abort(
      "PV_BIND_E072",
      "BRR-Fay target content requires the exact identity fixed-effect contrast and its canonical hash link.",
      "target_content"
    )
  }
  invisible(target_content)
}

pv_binding_legacy_revalidation_abort <- function(detail) {
  pv_binding_abort("PV_BIND_E080", detail, "target_content")
}

pv_binding_target_content_from_brr_target <- function(
  target,
  manifest,
  conf_level = NULL,
  data = NULL,
  formula = NULL
) {
  tryCatch(
    validate_pvstackr_brr_target(target),
    error = function(error) {
      if (inherits(error, "pvstackr_binding_error")) {
        stop(error)
      }
      pv_binding_target_abort("Current BRR-Fay target validation failed.")
    }
  )
  if (identical(target$schema_version, "0.2.0")) {
    pv_binding_target_manifest_validate(target$target_content, manifest)
    if (!is.null(conf_level) && !identical(as.double(conf_level), target$conf_level)) {
      pv_binding_target_abort(
        "Schema-0.2 target content confidence level cannot be overridden."
      )
    }
    return(target$target_content)
  }
  if (is.null(data) || is.null(formula)) {
    pv_binding_legacy_revalidation_abort(
      "Schema-0.1 BRR-Fay targets require complete original data and formula revalidation."
    )
  }
  if (!is.data.frame(data) || !inherits(formula, "formula")) {
    pv_binding_legacy_revalidation_abort(
      "Legacy target raw revalidation requires a data frame and formula."
    )
  }
  if (is.null(conf_level)) {
    conf_level <- target$conf_level
  }
  if (is.null(conf_level)) {
    pv_binding_legacy_revalidation_abort(
      "Legacy target revalidation requires an explicit or stored conf_level."
    )
  }
  retained_conf_level <- target$conf_level
  if (!is.null(retained_conf_level)) {
    supplied_conf_level <- tryCatch(
      pv_assert_probability(conf_level, "conf_level"),
      error = function(error) NULL
    )
    retained_conf_level <- tryCatch(
      pv_assert_probability(retained_conf_level, "target$conf_level"),
      error = function(error) NULL
    )
    if (is.null(supplied_conf_level) || is.null(retained_conf_level) ||
        !identical(as.double(supplied_conf_level), as.double(retained_conf_level))) {
      pv_binding_legacy_revalidation_abort(
        "The supplied confidence level does not match the value retained on the legacy target."
      )
    }
    conf_level <- retained_conf_level
  }
  role_fields <- c("pv_cols", "weight_col", "rep_weight_cols", "id_cols")
  if (!identical(target$schema_version, "0.1.0") ||
      !all(role_fields %in% names(target)) ||
      !is.character(target$pv_cols) || length(target$pv_cols) < 1L ||
      length(target$pv_cols) != target$M ||
      !is.null(attributes(target$pv_cols)) ||
      anyNA(target$pv_cols) || any(!nzchar(target$pv_cols)) ||
      anyDuplicated(target$pv_cols) ||
      !is.character(target$weight_col) || length(target$weight_col) != 1L ||
      !is.null(attributes(target$weight_col)) || is.na(target$weight_col) ||
      !nzchar(target$weight_col) ||
      !is.character(target$rep_weight_cols) || length(target$rep_weight_cols) < 2L ||
      !is.null(attributes(target$rep_weight_cols)) ||
      anyNA(target$rep_weight_cols) || any(!nzchar(target$rep_weight_cols)) ||
      anyDuplicated(target$rep_weight_cols) ||
      !(is.null(target$id_cols) ||
        (is.character(target$id_cols) && is.null(attributes(target$id_cols)) &&
          !anyNA(target$id_cols) &&
          all(nzchar(target$id_cols)) && !anyDuplicated(target$id_cols))) ||
      !identical(
        unname(vapply(target$per_pv, `[[`, character(1), "pv_col")),
        unname(target$pv_cols)
      )) {
    pv_binding_legacy_revalidation_abort(
      "Legacy target role specifications are missing, malformed, or inconsistent."
    )
  }
  if (!identical(target$engine, "lm")) {
    pv_binding_legacy_revalidation_abort(
      "The schema-0.1 target engine cannot be independently revalidated."
    )
  }

  pv_binding_manifest_validate(manifest)
  raw_revalidation <- tryCatch(
    {
      target_formula_ast <- pv_binding_formula_ast(target$formula, data = data)
      supplied_formula_ast <- pv_binding_formula_ast(formula, data = data)
      if (!identical(target_formula_ast, supplied_formula_ast)) {
        pv_binding_legacy_revalidation_abort(
          "The supplied formula does not exactly match the formula stored on the legacy target."
        )
      }
      fresh_bundle <- pv_binding_resolve_model_bundle(
        data = data,
        formula = formula
      )
      fe_names <- pv_normalize_fe_names(colnames(fresh_bundle$model_matrix))
      estimand_metadata <- list(
        estimand_id = "brr_fay_fixed_effects",
        target_source = target$target_source,
        target_engine_id = "lm_wls_brr_fay_v1",
        parameter_scope = "fixed_effect",
        fe_names = fe_names,
        interval_role = target$interval_role,
        coverage_claim_allowed = target$coverage_claim_allowed
      )
      fresh_manifest <- pv_binding_manifest_build(
        data = data,
        formula = formula,
        pv_cols = target$pv_cols,
        weight_col = target$weight_col,
        rep_weight_cols = target$rep_weight_cols,
        fay_k = target$fay_k,
        id_cols = target$id_cols,
        family_link = pv_binding_family_link_projection("gaussian", "identity"),
        estimand_contrast = NULL,
        estimand_metadata = estimand_metadata,
        model_bundle = fresh_bundle
      )
      fresh_manifest <- c(
        fresh_manifest,
        list(model_bundle_hash = fresh_bundle$bundle_hash)
      )
      fresh_manifest$manifest_hash <- pv_binding_hash_payload(
        pv_binding_manifest_hash_payload(fresh_manifest),
        "manifest"
      )
      pv_binding_manifest_validate(fresh_manifest)
      manifest_comparison <- pv_binding_manifest_compare(fresh_manifest, manifest)
      if (!isTRUE(manifest_comparison$ok) ||
          !identical(fresh_manifest$component_hashes, manifest$component_hashes)) {
        pv_binding_legacy_revalidation_abort(
          "The supplied manifest does not exactly match raw-input recomputation."
        )
      }

      df_complete_input <- if (identical(target$df_method, "classic")) {
        NULL
      } else {
        target$df_complete
      }
      rebuilt_target <- pv_brr_target(
        data = data,
        formula = formula,
        pv_cols = target$pv_cols,
        weight_col = target$weight_col,
        rep_weight_cols = target$rep_weight_cols,
        fay_k = target$fay_k,
        id_cols = target$id_cols,
        conf_level = conf_level,
        allow_m1 = identical(target$M, 1L),
        df_method = target$df_method,
        df_complete = df_complete_input,
        engine = "lm",
        verbose = FALSE
      )
      exact_fields <- c(
        "per_pv", pv_binding_target_derived_fields(), "M", "R", "fay_k",
        "df_method", "interval_role", "coverage_claim_allowed",
        "fe_names", "target_source", "engine", "policy"
      )
      rebuilt_legacy_hashes <- rebuilt_target$binding_manifest$legacy_hashes
      rebuilt_legacy_df_complete <- if (identical(rebuilt_target$df_method, "classic")) {
        stats::setNames(rep(NA_real_, length(rebuilt_target$fe_names)), rebuilt_target$fe_names)
      } else {
        rebuilt_target$df_complete
      }
      if (!identical(target[exact_fields], rebuilt_target[exact_fields]) ||
          !identical(target$df_complete, rebuilt_legacy_df_complete) ||
          !identical(target$design_hash, rebuilt_legacy_hashes$design_hash) ||
          !identical(target$target_hash, rebuilt_legacy_hashes$target_hash)) {
        pv_binding_legacy_revalidation_abort(
          "Legacy target primitives do not exactly match independent raw-input reconstruction."
        )
      }
      list(target = rebuilt_target, manifest = fresh_manifest)
    },
    pvstackr_binding_error = function(error) {
      if (identical(error$code, "PV_BIND_E080")) {
        stop(error)
      }
      pv_binding_legacy_revalidation_abort(
        "Legacy target raw-input recomputation did not complete successfully."
      )
    },
    error = function(error) {
      pv_binding_legacy_revalidation_abort(
        "Legacy target raw-input recomputation did not complete successfully."
      )
    }
  )

  rebuilt_target <- raw_revalidation$target
  df_complete <- if (identical(rebuilt_target$df_method, "classic")) {
    NULL
  } else {
    rebuilt_target$df_complete
  }
  content <- pv_binding_target_content_build(
    per_pv = rebuilt_target$per_pv,
    M = rebuilt_target$M,
    R = rebuilt_target$R,
    fay_k = rebuilt_target$fay_k,
    df_method = rebuilt_target$df_method,
    df_complete = df_complete,
    conf_level = conf_level,
    target_source = rebuilt_target$target_source,
    target_engine_id = "lm_wls_brr_fay_v1",
    target_policy = rebuilt_target$policy,
    manifest_hash = manifest$manifest_hash,
    stored_derived = rebuilt_target[pv_binding_target_derived_fields()]
  )
  pv_binding_target_manifest_validate(content, manifest)
  content
}

pv_binding_raw_join <- function(parts) {
  if (length(parts) == 0L) {
    return(raw())
  }
  if (!all(vapply(parts, is.raw, logical(1)))) {
    pv_binding_abort("PV_BIND_E081", "Canonical byte parts must all be raw vectors.", "canonicalizer")
  }
  do.call(c, unname(parts))
}

pv_binding_u64be <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 0 || x != floor(x) || x > .Machine$integer.max) {
    pv_binding_abort(
      "PV_BIND_E081",
      "A TLV length must be an exact non-negative integer within the supported range.",
      "canonicalizer"
    )
  }
  value <- as.double(x)
  out <- raw(8L)
  for (index in 8:1) {
    out[[index]] <- as.raw(value %% 256)
    value <- floor(value / 256)
  }
  out
}

pv_binding_utf8_raw <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    pv_binding_abort("PV_BIND_E081", "A canonical string must be a non-missing character scalar.", "canonicalizer")
  }
  value <- enc2utf8(x)
  if (!validUTF8(value)) {
    pv_binding_abort("PV_BIND_E081", "A canonical string must have valid UTF-8 encoding.", "canonicalizer")
  }
  charToRaw(value)
}

pv_binding_tlv <- function(tag, payload = raw()) {
  tag_raw <- pv_binding_utf8_raw(tag)
  if (!is.raw(payload)) {
    pv_binding_abort("PV_BIND_E081", "A TLV payload must be a raw vector.", "canonicalizer")
  }
  pv_binding_raw_join(list(
    pv_binding_u64be(length(tag_raw)),
    tag_raw,
    pv_binding_u64be(length(payload)),
    payload
  ))
}

pv_binding_names_bytes <- function(x) {
  value_names <- names(x)
  if (is.null(value_names)) {
    return(pv_binding_tlv("names_absent"))
  }
  if (length(value_names) != length(x) || anyNA(value_names)) {
    pv_binding_abort("PV_BIND_E081", "Vector names must be complete and non-missing when present.", "canonicalizer")
  }
  encoded <- lapply(value_names, function(value) {
    pv_binding_tlv("name", pv_binding_utf8_raw(value))
  })
  pv_binding_tlv(
    "names_present",
    pv_binding_raw_join(c(list(pv_binding_u64be(length(value_names))), encoded))
  )
}

pv_binding_vector_bytes <- function(tag, x, encode_one) {
  elements <- lapply(seq_along(x), function(index) encode_one(x[[index]]))
  payload <- pv_binding_raw_join(c(
    list(pv_binding_u64be(length(x)), pv_binding_names_bytes(x)),
    elements
  ))
  pv_binding_tlv(tag, payload)
}

pv_binding_character_bytes <- function(x) {
  pv_binding_vector_bytes("character", x, function(value) {
    if (is.na(value)) {
      pv_binding_tlv("na")
    } else {
      pv_binding_tlv("value", pv_binding_utf8_raw(value))
    }
  })
}

pv_binding_logical_bytes <- function(x) {
  pv_binding_vector_bytes("logical", x, function(value) {
    if (is.na(value)) {
      pv_binding_tlv("na")
    } else if (isTRUE(value)) {
      pv_binding_tlv("true")
    } else {
      pv_binding_tlv("false")
    }
  })
}

pv_binding_integer_bytes <- function(x) {
  pv_binding_vector_bytes("integer", x, function(value) {
    if (is.na(value)) {
      pv_binding_tlv("na")
    } else {
      pv_binding_tlv("value", writeBin(as.integer(value), raw(), size = 4L, endian = "big"))
    }
  })
}

pv_binding_double_bytes <- function(x) {
  pv_binding_vector_bytes("double", x, function(value) {
    if (is.nan(value)) {
      pv_binding_tlv("nan")
    } else if (is.na(value)) {
      pv_binding_tlv("na")
    } else if (is.infinite(value)) {
      pv_binding_tlv(if (value > 0) "positive_infinity" else "negative_infinity")
    } else {
      if (value == 0) {
        value <- 0
      }
      pv_binding_tlv("value", writeBin(as.double(value), raw(), size = 8L, endian = "big"))
    }
  })
}

pv_binding_factor_bytes <- function(x) {
  payload <- list(
    ordered = is.ordered(x),
    levels = levels(x),
    codes = as.integer(x)
  )
  pv_binding_tlv("factor", pv_binding_canonical_bytes(payload))
}

pv_binding_date_bytes <- function(x) {
  payload <- list(class = class(x), values = as.double(x))
  pv_binding_tlv("date", pv_binding_canonical_bytes(payload))
}

pv_binding_posixct_bytes <- function(x) {
  timezone <- attr(x, "tzone", exact = TRUE)
  if (is.null(timezone)) {
    timezone <- character()
  }
  payload <- list(class = class(x), timezone = timezone, values = as.double(x))
  pv_binding_tlv("posixct", pv_binding_canonical_bytes(payload))
}

pv_binding_matrix_bytes <- function(x) {
  extra_attributes <- setdiff(names(attributes(x)), c("dim", "dimnames"))
  if (length(extra_attributes) > 0L) {
    pv_binding_abort(
      "PV_BIND_E081",
      "Matrix attributes must be extracted into explicit schema fields before canonicalization.",
      "canonicalizer"
    )
  }
  if (length(dim(x)) != 2L) {
    pv_binding_abort("PV_BIND_E081", "Only two-dimensional matrices are supported.", "canonicalizer")
  }
  values <- as.vector(x)
  names(values) <- NULL
  payload <- list(
    storage_type = typeof(x),
    dimensions = as.integer(dim(x)),
    dimnames = dimnames(x),
    values = values
  )
  pv_binding_tlv("matrix", pv_binding_canonical_bytes(payload))
}

pv_binding_list_bytes <- function(x) {
  value_names <- names(x)
  expected_attributes <- if (is.null(value_names)) NULL else list(names = value_names)
  if (!identical(attributes(x), expected_attributes)) {
    pv_binding_abort(
      "PV_BIND_E081",
      "Canonical lists may carry only their exact names attribute.",
      "canonicalizer"
    )
  }
  if (!is.null(value_names)) {
    if (length(value_names) != length(x) || anyNA(value_names) || any(!nzchar(value_names)) || anyDuplicated(value_names)) {
      pv_binding_abort(
        "PV_BIND_E081",
        "A named canonical list must have complete, non-empty, unique field names in schema order.",
        "canonicalizer"
      )
    }
  }
  elements <- lapply(x, pv_binding_canonical_bytes)
  payload <- pv_binding_raw_join(c(
    list(pv_binding_u64be(length(x)), pv_binding_names_bytes(x)),
    elements
  ))
  pv_binding_tlv("list", payload)
}

pv_binding_assert_atomic_attributes <- function(x) {
  value_names <- names(x)
  expected_attributes <- if (is.null(value_names)) NULL else list(names = value_names)
  if (!identical(attributes(x), expected_attributes)) {
    pv_binding_abort(
      "PV_BIND_E081",
      "Atomic values may carry only their exact names attribute before canonicalization.",
      "canonicalizer"
    )
  }
  invisible(x)
}

pv_binding_canonical_bytes <- function(x) {
  if (is.null(x)) {
    return(pv_binding_tlv("null"))
  }
  if (is.data.frame(x)) {
    pv_binding_abort("PV_BIND_E081", "Data frames require an explicit schema projection before canonicalization.", "canonicalizer")
  }
  if (is.matrix(x)) {
    return(pv_binding_matrix_bytes(x))
  }
  if (is.ordered(x) || is.factor(x)) {
    return(pv_binding_factor_bytes(x))
  }
  if (inherits(x, "POSIXct")) {
    return(pv_binding_posixct_bytes(x))
  }
  if (inherits(x, "Date")) {
    return(pv_binding_date_bytes(x))
  }
  if (is.raw(x)) {
    pv_binding_assert_atomic_attributes(x)
    return(pv_binding_tlv("raw", pv_binding_raw_join(list(
      pv_binding_u64be(length(x)),
      pv_binding_names_bytes(x),
      x
    ))))
  }
  if (is.character(x)) {
    pv_binding_assert_atomic_attributes(x)
    return(pv_binding_character_bytes(x))
  }
  if (is.logical(x)) {
    pv_binding_assert_atomic_attributes(x)
    return(pv_binding_logical_bytes(x))
  }
  if (is.integer(x)) {
    pv_binding_assert_atomic_attributes(x)
    return(pv_binding_integer_bytes(x))
  }
  if (is.double(x)) {
    pv_binding_assert_atomic_attributes(x)
    return(pv_binding_double_bytes(x))
  }
  if (is.list(x) && is.null(attr(x, "class", exact = TRUE))) {
    return(pv_binding_list_bytes(x))
  }
  pv_binding_abort(
    "PV_BIND_E081",
    "The value has a type or class unsupported by the binding canonicalizer.",
    "canonicalizer"
  )
}

pv_binding_sha256_raw <- function(bytes) {
  if (!is.raw(bytes)) {
    pv_binding_abort("PV_BIND_E081", "SHA-256 input must be one raw byte vector.", "hash")
  }
  if (!requireNamespace("digest", quietly = TRUE)) {
    pv_binding_abort("PV_BIND_E004", "The SHA-256 provider is unavailable.", "hash")
  }
  paste0("sha256:", digest::digest(bytes, algo = "sha256", serialize = FALSE))
}

pv_binding_hash_payload <- function(payload, domain) {
  allowed <- pv_binding_hash_domains()
  if (!is.character(domain) || length(domain) != 1L || is.na(domain) || !domain %in% allowed) {
    pv_binding_abort("PV_BIND_E004", "The binding hash domain is unsupported.", "hash")
  }
  prefix <- charToRaw(paste0("pvstackr-binding-v1/", domain))
  bytes <- pv_binding_raw_join(list(
    prefix,
    as.raw(0L),
    pv_binding_canonical_bytes(payload)
  ))
  pv_binding_sha256_raw(bytes)
}
