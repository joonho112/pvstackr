pv_compare_collect_fits <- function(..., fits = NULL) {
  dots <- list(...)
  if (!is.null(fits)) {
    if (length(dots) > 0L) {
      pv_abort("Supply fits either through `...` or `fits`, not both.")
    }
    dots <- fits
  }
  if (length(dots) == 1L && is.list(dots[[1L]]) && !inherits(dots[[1L]], "pvstackr_fit")) {
    dots <- dots[[1L]]
  }
  if (!is.list(dots) || length(dots) < 2L) {
    pv_abort("`pv_compare_methods()` requires at least two `pvstackr_fit` objects.")
  }
  for (i in seq_along(dots)) {
    if (!inherits(dots[[i]], "pvstackr_fit")) {
      pv_abort("Every method comparison input must be a `pvstackr_fit` object.")
    }
    validate_pvstackr_fit(dots[[i]])
  }
  labels <- names(dots)
  if (is.null(labels)) {
    labels <- rep("", length(dots))
  }
  missing <- !nzchar(labels)
  labels[missing] <- vapply(dots[missing], `[[`, character(1), "method")
  if (anyDuplicated(labels)) {
    labels <- make.unique(labels, sep = "_")
  }
  names(dots) <- labels
  dots
}

pv_compare_terms <- function(fits) {
  terms <- unique(unlist(lapply(fits, function(fit) {
    if (is.data.frame(fit$estimates) && "term" %in% names(fit$estimates)) {
      as.character(fit$estimates$term)
    } else {
      character()
    }
  }), use.names = FALSE))
  if (length(terms) == 0L) {
    pv_abort("At least one compared fit must contain reportable estimate terms.")
  }
  terms
}

pv_compare_reference_label <- function(fits, reference_method = NULL) {
  labels <- names(fits)
  methods <- vapply(fits, `[[`, character(1), "method")
  statuses <- vapply(fits, `[[`, character(1), "status")
  if (!is.null(reference_method)) {
    reference_method <- pv_assert_scalar_string(reference_method, "reference_method")
    label_hit <- labels == reference_method
    method_hit <- methods == reference_method
    hits <- which(label_hit | method_hit)
    if (length(hits) != 1L) {
      pv_abort("`reference_method` must identify exactly one compared fit by label or method.")
    }
    if (identical(statuses[[hits]], "blocked")) {
      pv_abort("`reference_method` must not identify a blocked fit.")
    }
    return(labels[[hits]])
  }
  preferred <- which(methods == "per_pv" & statuses != "blocked")
  if (length(preferred) > 0L) {
    return(labels[[preferred[[1L]]]])
  }
  available <- which(statuses != "blocked")
  if (length(available) == 0L) {
    pv_abort("At least one compared fit must be non-blocked to serve as a reference.")
  }
  labels[[available[[1L]]]]
}

pv_compare_fit_n <- function(fit) {
  if (identical(fit$method, "per_pv") && !is.null(fit$target$M)) {
    return(as.integer(fit$target$M))
  }
  if (!is.null(fit$stack_fit$meta$n_fits)) {
    return(as.integer(fit$stack_fit$meta$n_fits))
  }
  if (!is.null(fit$diagnostics$reference$M)) {
    return(as.integer(fit$diagnostics$reference$M))
  }
  if (!is.null(fit$diagnostics$psis)) {
    return(1L)
  }
  NA_integer_
}

pv_compare_normalize_timings <- function(timings) {
  if (is.null(timings)) {
    return(NULL)
  }
  if (is.list(timings) && !is.data.frame(timings)) {
    timings <- unlist(timings, use.names = TRUE)
  }
  if (!is.numeric(timings)) {
    pv_abort("`timings` must be a named numeric vector or list.")
  }
  if (is.null(names(timings)) || any(!nzchar(names(timings)))) {
    pv_abort("`timings` must be named by method label or method id.")
  }
  if (any(!is.finite(timings) | timings < 0)) {
    pv_abort("`timings` values must be finite non-negative elapsed seconds.")
  }
  timings
}

pv_compare_elapsed <- function(label, method, timings) {
  if (is.null(timings)) {
    return(NA_real_)
  }
  if (label %in% names(timings)) {
    return(unname(timings[[label]]))
  }
  if (method %in% names(timings)) {
    return(unname(timings[[method]]))
  }
  NA_real_
}

pv_compare_fit_target_hash <- function(fit) {
  if (!is.null(fit$target$target_hash)) {
    return(fit$target$target_hash)
  }
  if (!is.null(fit$ccc$target_hash)) {
    return(fit$ccc$target_hash)
  }
  if (!is.null(fit$diagnostics$preflight$target_hash)) {
    return(fit$diagnostics$preflight$target_hash)
  }
  if (!is.null(fit$diagnostics$pooling$target_hash)) {
    return(fit$diagnostics$pooling$target_hash)
  }
  if (!is.null(fit$provenance$target_hash)) {
    return(fit$provenance$target_hash)
  }
  if (is.data.frame(fit$estimates) && "target_hash" %in% names(fit$estimates) && nrow(fit$estimates) > 0L) {
    hashes <- unique(as.character(fit$estimates$target_hash))
    hashes <- hashes[nzchar(hashes) & !is.na(hashes)]
    if (length(hashes) > 0L) {
      return(hashes[[1L]])
    }
  }
  NA_character_
}

pv_compare_fit_target_source <- function(fit) {
  if (!is.null(fit$target$target_source)) {
    return(fit$target$target_source)
  }
  if (!is.null(fit$ccc$target_source)) {
    return(fit$ccc$target_source)
  }
  if (!is.null(fit$diagnostics$preflight$target_source)) {
    return(fit$diagnostics$preflight$target_source)
  }
  if (!is.null(fit$diagnostics$pooling$target_source)) {
    return(fit$diagnostics$pooling$target_source)
  }
  if (is.data.frame(fit$estimates) && "target_source" %in% names(fit$estimates) && nrow(fit$estimates) > 0L) {
    sources <- unique(as.character(fit$estimates$target_source))
    sources <- sources[nzchar(sources) & !is.na(sources)]
    if (length(sources) > 0L) {
      return(paste(sources, collapse = ","))
    }
  }
  NA_character_
}

pv_compare_fit_pooling_hash <- function(fit) {
  if (!is.null(fit$target$pooling_hash)) {
    return(fit$target$pooling_hash)
  }
  if (!is.null(fit$diagnostics$pooling$pooling_hash)) {
    return(fit$diagnostics$pooling$pooling_hash)
  }
  if (!is.null(fit$target$provenance$pooling_hash)) {
    return(fit$target$provenance$pooling_hash)
  }
  if (!is.null(fit$provenance$pooling_hash)) {
    return(fit$provenance$pooling_hash)
  }
  if (!is.null(fit$diagnostics$pooling$target_hash)) {
    return(fit$diagnostics$pooling$target_hash)
  }
  if (is.data.frame(fit$estimates) && "pooling_hash" %in% names(fit$estimates) && nrow(fit$estimates) > 0L) {
    hashes <- unique(as.character(fit$estimates$pooling_hash))
    hashes <- hashes[nzchar(hashes) & !is.na(hashes)]
    if (length(hashes) > 0L) {
      return(hashes[[1L]])
    }
  }
  NA_character_
}

pv_compare_fit_pooling_source <- function(fit) {
  if (!is.null(fit$diagnostics$pooling$pooling_source)) {
    return(fit$diagnostics$pooling$pooling_source)
  }
  if (is.data.frame(fit$estimates) && "pooling_source" %in% names(fit$estimates) && nrow(fit$estimates) > 0L) {
    sources <- unique(as.character(fit$estimates$pooling_source))
    sources <- sources[nzchar(sources) & !is.na(sources)]
    if (length(sources) > 0L) {
      return(paste(sources, collapse = ","))
    }
  }
  if (!is.null(fit$diagnostics$pooling$target_source)) {
    return(fit$diagnostics$pooling$target_source)
  }
  if (identical(fit$method, "per_pv") && !is.null(fit$target$target_source)) {
    return(fit$target$target_source)
  }
  NA_character_
}

pv_compare_shared_hash <- function(hashes) {
  hashes <- as.character(hashes)
  out <- rep(FALSE, length(hashes))
  valid <- !is.na(hashes) & nzchar(hashes)
  if (sum(valid) < 2L) {
    return(out)
  }
  counts <- table(hashes[valid])
  shared <- names(counts[counts > 1L])
  out[valid] <- hashes[valid] %in% shared
  out
}

pv_compare_shared_external_target <- function(target_sources, target_hashes) {
  target_sources <- as.character(target_sources)
  target_hashes <- as.character(target_hashes)
  out <- rep(FALSE, length(target_hashes))
  valid <- !is.na(target_hashes) & nzchar(target_hashes) &
    target_sources == "external_brr_fay_rubin"
  if (sum(valid) < 2L) {
    return(out)
  }
  counts <- table(target_hashes[valid])
  shared <- names(counts[counts > 1L])
  out[valid] <- target_hashes[valid] %in% shared
  out
}

pv_compare_shares_reference_hash <- function(hashes, reference_label, labels) {
  hashes <- as.character(hashes)
  out <- rep(FALSE, length(hashes))
  idx <- match(reference_label, labels)
  if (is.na(idx) || idx < 1L || idx > length(hashes)) {
    return(out)
  }
  reference_hash <- hashes[[idx]]
  if (is.na(reference_hash) || !nzchar(reference_hash)) {
    return(out)
  }
  valid <- !is.na(hashes) & nzchar(hashes)
  out[valid] <- hashes[valid] == reference_hash
  out[[idx]] <- FALSE
  out
}

pv_compare_fit_center <- function(fit) {
  if (!identical(fit$method, "stack_direct")) {
    return(NA_character_)
  }
  center <- fit$control$center %||% fit$ccc$center
  if (!is.null(center)) {
    return(as.character(center))
  }
  NA_character_
}

pv_compare_shared_values <- function(values) {
  values <- as.character(values)
  valid <- !is.na(values) & nzchar(values)
  if (sum(valid) < 2L) {
    return(character())
  }
  counts <- table(values[valid])
  names(counts[counts > 1L])
}

pv_compare_target_overlap_diagnostics <- function(method_diagnostics) {
  shared_target_sources <- pv_compare_shared_values(method_diagnostics$target_source)
  list(
    shared_external_target = any(method_diagnostics$shared_external_target),
    shared_target_hash = any(method_diagnostics$shared_target_hash),
    shared_pooling_hash = any(method_diagnostics$shared_pooling_hash),
    shares_reference_target = any(method_diagnostics$shares_reference_target),
    shares_reference_pooling = any(method_diagnostics$shares_reference_pooling),
    shared_target_sources = shared_target_sources,
    independence_caveat_required = any(
      method_diagnostics$shared_external_target |
        method_diagnostics$shared_target_hash |
        method_diagnostics$shared_pooling_hash |
        method_diagnostics$shares_reference_target |
        method_diagnostics$shares_reference_pooling
    ) || length(shared_target_sources) > 0L,
    independence_caveat = paste(
      "Agreement bands are descriptive; shared target, pooling, or source",
      "metadata should not be read as independent corroboration."
    )
  )
}

pv_compare_fit_psis_status <- function(fit) {
  if (!is.null(fit$diagnostics$psis$status)) {
    return(fit$diagnostics$psis$status)
  }
  if (is.data.frame(fit$estimates) && "psis_status" %in% names(fit$estimates) && nrow(fit$estimates) > 0L) {
    statuses <- unique(as.character(fit$estimates$psis_status))
    statuses <- statuses[nzchar(statuses) & !is.na(statuses)]
    if (length(statuses) > 0L) {
      return(statuses[[1L]])
    }
  }
  NA_character_
}

pv_compare_fit_pareto_k_max <- function(fit) {
  if (!is.null(fit$diagnostics$psis$pareto_k_max)) {
    return(as.numeric(fit$diagnostics$psis$pareto_k_max))
  }
  if (is.data.frame(fit$estimates) && "pareto_k_max" %in% names(fit$estimates) && nrow(fit$estimates) > 0L) {
    values <- unique(as.numeric(fit$estimates$pareto_k_max))
    values <- values[is.finite(values)]
    if (length(values) > 0L) {
      return(values[[1L]])
    }
  }
  NA_real_
}

pv_compare_fit_field <- function(fit, field, type = c("character", "logical")) {
  type <- match.arg(type)
  if (!is.data.frame(fit$estimates) || !field %in% names(fit$estimates) || nrow(fit$estimates) == 0L) {
    return(if (identical(type, "logical")) NA else NA_character_)
  }
  values <- fit$estimates[[field]]
  values <- values[!is.na(values)]
  if (length(values) == 0L) {
    return(if (identical(type, "logical")) NA else NA_character_)
  }
  if (identical(type, "logical")) {
    return(all(as.logical(values)))
  }
  paste(unique(as.character(values)), collapse = ",")
}

pv_compare_row_field <- function(row, field, type = c("character", "logical")) {
  type <- match.arg(type)
  if (!is.data.frame(row) || nrow(row) != 1L || !field %in% names(row) || is.na(row[[field]][[1L]])) {
    return(if (identical(type, "logical")) NA else NA_character_)
  }
  if (identical(type, "logical")) {
    return(as.logical(row[[field]][[1L]]))
  }
  as.character(row[[field]][[1L]])
}

pv_compare_row_numeric <- function(row, field) {
  if (!is.data.frame(row) || nrow(row) != 1L ||
      !field %in% names(row) || is.na(row[[field]][[1L]])) {
    return(NA_real_)
  }
  as.numeric(row[[field]][[1L]])
}

pv_compare_one_row <- function(fit, label, term, reference, reference_label) {
  estimates <- fit$estimates
  row <- if (is.data.frame(estimates) && "term" %in% names(estimates)) {
    estimates[match(term, estimates$term), , drop = FALSE]
  } else {
    data.frame()
  }
  has_row <- nrow(row) == 1L && !is.na(row$term)
  ref_row <- reference[match(term, reference$term), , drop = FALSE]
  ref_has_row <- nrow(ref_row) == 1L && !is.na(ref_row$term)

  estimate <- if (has_row) row$estimate else NA_real_
  se <- if (has_row && "se" %in% names(row)) row$se else NA_real_
  df <- if (has_row && "df" %in% names(row)) row$df else NA_real_
  df_method <- pv_compare_row_field(row, "df_method", "character")
  df_complete <- pv_compare_row_numeric(row, "df_complete")
  conf_level <- pv_compare_row_numeric(row, "conf_level")
  interval_role <- pv_compare_row_field(row, "interval_role", "character")
  coverage_claim_allowed <- pv_compare_row_field(row, "coverage_claim_allowed", "logical")
  conf_low <- if (has_row && "conf_low" %in% names(row)) row$conf_low else NA_real_
  conf_high <- if (has_row && "conf_high" %in% names(row)) row$conf_high else NA_real_
  ref_estimate <- if (ref_has_row) ref_row$estimate else NA_real_
  ref_se <- if (ref_has_row && "se" %in% names(ref_row)) ref_row$se else NA_real_

  estimate_diff <- estimate - ref_estimate
  se_ratio <- se / ref_se
  abs_z_diff <- abs(estimate_diff) / sqrt(se^2 + ref_se^2)
  agreement_band <- if (identical(fit$status, "blocked") || !is.finite(abs_z_diff)) {
    "not_available"
  } else if (abs_z_diff < 0.1) {
    "close"
  } else if (abs_z_diff < 0.5) {
    "moderate"
  } else {
    "different"
  }

  data.frame(
    method = fit$method,
    method_label = label,
    term = term,
    status = fit$status,
    estimate = estimate,
    se = se,
    std.error = se,
    df = df,
    df_method = df_method,
    df_complete = df_complete,
    conf_level = conf_level,
    interval_role = interval_role,
    coverage_claim_allowed = coverage_claim_allowed,
    conf_low = conf_low,
    conf_high = conf_high,
    reference_method = reference_label,
    reference_estimate = ref_estimate,
    reference_se = ref_se,
    estimate_diff = estimate_diff,
    se_ratio = se_ratio,
    abs_z_diff = abs_z_diff,
    agreement_band = agreement_band,
    reason_codes = paste(fit$reason_codes, collapse = ","),
    stringsAsFactors = FALSE
  )
}

pv_compare_table <- function(fits, terms, reference_label) {
  reference <- fits[[reference_label]]$estimates
  rows <- list()
  idx <- 1L
  for (label in names(fits)) {
    fit <- fits[[label]]
    for (term in terms) {
      rows[[idx]] <- pv_compare_one_row(fit, label, term, reference, reference_label)
      idx <- idx + 1L
    }
  }
  do.call(rbind, rows)
}

pv_compare_diagnostics <- function(table, fits, reference_label, timings) {
  labels <- names(fits)
  status <- vapply(fits, `[[`, character(1), "status")
  methods <- vapply(fits, `[[`, character(1), "method")
  timing <- data.frame(
    method = methods,
    method_label = labels,
    status = status,
    n_fits = vapply(fits, pv_compare_fit_n, integer(1)),
    elapsed_seconds = vapply(seq_along(fits), function(i) {
      pv_compare_elapsed(labels[[i]], methods[[i]], timings)
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
  by_method <- lapply(labels, function(label) {
    x <- table[table$method_label == label, , drop = FALSE]
    finite_diff <- x$abs_z_diff[is.finite(x$abs_z_diff)]
    finite_ratio <- abs(log(x$se_ratio[is.finite(x$se_ratio) & x$se_ratio > 0]))
    data.frame(
      method = x$method[[1L]],
      method_label = x$method_label[[1L]],
      status = x$status[[1L]],
      max_abs_z_diff = if (length(finite_diff) == 0L) NA_real_ else max(finite_diff),
      max_abs_log_se_ratio = if (length(finite_ratio) == 0L) NA_real_ else max(finite_ratio),
      n_terms = length(unique(x$term)),
      n_available = sum(is.finite(x$estimate)),
      stringsAsFactors = FALSE
    )
  })
  agreement <- do.call(rbind, by_method)
  target_source <- vapply(fits, pv_compare_fit_target_source, character(1))
  target_hash <- vapply(fits, pv_compare_fit_target_hash, character(1))
  pooling_source <- vapply(fits, pv_compare_fit_pooling_source, character(1))
  pooling_hash <- vapply(fits, pv_compare_fit_pooling_hash, character(1))
  shared_target_hash <- pv_compare_shared_hash(target_hash)
  shared_pooling_hash <- pv_compare_shared_hash(pooling_hash)
  shared_external_target <- pv_compare_shared_external_target(target_source, target_hash)
  shares_reference_target <- pv_compare_shares_reference_hash(target_hash, reference_label, labels)
  shares_reference_pooling <- pv_compare_shares_reference_hash(pooling_hash, reference_label, labels)
  center <- vapply(fits, pv_compare_fit_center, character(1))
  method_diagnostics <- data.frame(
    method = methods,
    method_label = labels,
    status = status,
    reason_codes = vapply(fits, function(fit) paste(fit$reason_codes, collapse = ","), character(1)),
    warning_count = vapply(fits, function(fit) length(fit$warnings), integer(1)),
    n_terms = agreement$n_terms,
    n_available = agreement$n_available,
    df_method = vapply(fits, pv_compare_fit_field, character(1), field = "df_method", type = "character"),
    interval_role = vapply(fits, pv_compare_fit_field, character(1), field = "interval_role", type = "character"),
    coverage_claim_allowed = vapply(fits, pv_compare_fit_field, logical(1), field = "coverage_claim_allowed", type = "logical"),
    n_descriptive_intervals = vapply(fits, function(fit) {
      if (!is.data.frame(fit$estimates) || !"coverage_claim_allowed" %in% names(fit$estimates) || nrow(fit$estimates) == 0L) {
        return(0L)
      }
      sum(!is.na(fit$estimates$coverage_claim_allowed) & !as.logical(fit$estimates$coverage_claim_allowed))
    }, integer(1)),
    target_source = target_source,
    target_hash = target_hash,
    pooling_source = pooling_source,
    pooling_hash = pooling_hash,
    center = center,
    shared_target_hash = shared_target_hash,
    shared_pooling_hash = shared_pooling_hash,
    shared_external_target = shared_external_target,
    shares_reference_target = shares_reference_target,
    shares_reference_pooling = shares_reference_pooling,
    psis_status = vapply(fits, pv_compare_fit_psis_status, character(1)),
    pareto_k_max = vapply(fits, pv_compare_fit_pareto_k_max, numeric(1)),
    n_fits = timing$n_fits,
    elapsed_seconds = timing$elapsed_seconds,
    stringsAsFactors = FALSE
  )
  list(
    reference_method = reference_label,
    methods = methods,
    statuses = status,
    blocked_methods = labels[status == "blocked"],
    warning_methods = labels[status == "warning"],
    agreement = agreement,
    method_diagnostics = method_diagnostics,
    timing = timing,
    target_overlap = pv_compare_target_overlap_diagnostics(method_diagnostics)
  )
}

new_pvstackr_method_comparison <- function(table, diagnostics, fits, reference_method, include_fits = FALSE) {
  comparison <- list(
    table = table,
    estimate_table = table,
    diagnostics = diagnostics,
    diagnostic_table = diagnostics$method_diagnostics,
    agreement = diagnostics$agreement,
    timing = diagnostics$timing,
    fits = if (isTRUE(include_fits)) fits else NULL,
    reference_method = reference_method,
    methods = vapply(fits, `[[`, character(1), "method"),
    method_labels = names(fits),
    created_at = as.character(Sys.time()),
    schema_version = pv_schema_version(),
    provenance = pv_provenance("pv_compare_methods"),
    warnings = character()
  )
  class(comparison) <- c("pvstackr_method_comparison", "list")
  validate_pvstackr_method_comparison(comparison)
  comparison
}

pv_numeric_equal_na <- function(x, y, tol = 0) {
  if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y)) {
    return(FALSE)
  }
  both_na <- is.na(x) & is.na(y)
  both_value <- !is.na(x) & !is.na(y) & abs(x - y) <= tol
  all(both_na | both_value)
}

pv_validate_comparison_labels <- function(comparison) {
  labels <- comparison$method_labels
  methods <- comparison$methods
  if (!is.character(labels) || length(labels) < 2L ||
      any(is.na(labels)) || any(!nzchar(labels)) || anyDuplicated(labels)) {
    pv_abort("Method comparison `method_labels` must contain at least two unique non-empty labels.")
  }
  if (!is.character(methods) || length(methods) != length(labels) ||
      any(!methods %in% pv_allowed_methods())) {
    pv_abort("Method comparison `methods` must use implemented method IDs.")
  }
  if (!identical(names(methods), labels)) {
    pv_abort("Method comparison `methods` names must align with method labels.")
  }
  labels
}

pv_validate_comparison_table <- function(comparison, labels) {
  if (!is.data.frame(comparison$table) || nrow(comparison$table) == 0L) {
    pv_abort("Method comparison `table` must be a non-empty data frame.")
  }
  if (!identical(comparison$table, comparison$estimate_table)) {
    pv_abort("Method comparison `estimate_table` must match `table`.")
  }
  table_required <- c(
    "method", "method_label", "term", "status", "estimate", "se",
    "std.error", "df", "df_method", "df_complete", "conf_level",
    "interval_role", "coverage_claim_allowed", "conf_low", "conf_high",
    "reference_method", "reference_estimate", "reference_se",
    "estimate_diff", "se_ratio", "abs_z_diff", "agreement_band",
    "reason_codes"
  )
  missing_table <- setdiff(table_required, names(comparison$table))
  if (length(missing_table) > 0L) {
    pv_abort(sprintf("Method comparison table is missing required column(s): %s.", paste(missing_table, collapse = ", ")))
  }
  table <- comparison$table
  if (any(!as.character(table$method) %in% pv_allowed_methods())) {
    pv_abort("Method comparison table must use implemented method IDs.")
  }
  if (any(!as.character(table$status) %in% c("ok", "warning", "blocked"))) {
    pv_abort("Method comparison table status must be one of: ok, warning, blocked.")
  }
  if (!is.logical(table$coverage_claim_allowed)) {
    pv_abort("Method comparison table `coverage_claim_allowed` must be logical.")
  }
  if (!is.numeric(table$df_complete) || !is.numeric(table$conf_level)) {
    pv_abort("Method comparison table interval df and confidence-level metadata must be numeric.")
  }
  if (any(!is.na(table$conf_level) & (!is.finite(table$conf_level) |
      table$conf_level <= 0 | table$conf_level >= 1))) {
    pv_abort("Method comparison table `conf_level` must be in (0, 1) or NA.")
  }
  if (any(!as.character(table$method_label) %in% labels)) {
    pv_abort("Method comparison table contains unknown method labels.")
  }
  if (!all(as.character(table$reference_method) == comparison$reference_method)) {
    pv_abort("Method comparison table reference metadata must match `reference_method`.")
  }
  if (!pv_numeric_equal_na(table$std.error, table$se)) {
    pv_abort("Method comparison table `std.error` must match `se`.")
  }
  if (any(!as.character(table$agreement_band) %in% c("not_available", "close", "moderate", "different"))) {
    pv_abort("Method comparison table has unknown agreement bands.")
  }
  for (label in labels) {
    method_rows <- table[table$method_label == label, , drop = FALSE]
    if (nrow(method_rows) == 0L) {
      pv_abort("Method comparison table must contain rows for every method label.")
    }
    expected_method <- unname(comparison$methods[[label]])
    if (!all(as.character(method_rows$method) == expected_method)) {
      pv_abort("Method comparison table method IDs must align with method labels.")
    }
  }
  terms <- unique(as.character(table$term))
  if (length(terms) == 0L || any(is.na(terms)) || any(!nzchar(terms))) {
    pv_abort("Method comparison table must contain non-empty terms.")
  }
  grid_counts <- table(as.character(table$method_label), as.character(table$term))
  if (!setequal(rownames(grid_counts), labels) ||
      !setequal(colnames(grid_counts), terms) ||
      any(grid_counts != 1L)) {
    pv_abort("Method comparison table must contain exactly one row for each method label and term.")
  }
  invisible(table)
}

pv_validate_comparison_diagnostics <- function(comparison, labels) {
  pv_assert_named_list(comparison$diagnostics, "diagnostics")
  required <- c(
    "reference_method", "methods", "statuses", "blocked_methods",
    "warning_methods", "agreement", "method_diagnostics", "timing",
    "target_overlap"
  )
  missing <- setdiff(required, names(comparison$diagnostics))
  if (length(missing) > 0L) {
    pv_abort(sprintf("Method comparison diagnostics are missing required field(s): %s.", paste(missing, collapse = ", ")))
  }
  agreement <- comparison$diagnostics$agreement
  method_diagnostics <- comparison$diagnostics$method_diagnostics
  timing <- comparison$diagnostics$timing
  if (!is.data.frame(agreement) || nrow(agreement) == 0L ||
      !is.data.frame(method_diagnostics) || nrow(method_diagnostics) == 0L ||
      !is.data.frame(timing) || nrow(timing) == 0L) {
    pv_abort("Method comparison diagnostics must include non-empty agreement, method, and timing data frames.")
  }
  if (!identical(comparison$agreement, agreement) ||
      !identical(comparison$diagnostic_table, method_diagnostics) ||
      !identical(comparison$timing, timing)) {
    pv_abort("Method comparison top-level diagnostic fields must match nested diagnostics.")
  }
  if (!identical(as.character(timing$method_label), labels) ||
      !identical(as.character(agreement$method_label), labels) ||
      !identical(as.character(method_diagnostics$method_label), labels)) {
    pv_abort("Method comparison diagnostic tables must align with method labels.")
  }
  if (!all(as.character(timing$method) %in% pv_allowed_methods()) ||
      !all(as.character(agreement$method) %in% pv_allowed_methods()) ||
      !all(as.character(method_diagnostics$method) %in% pv_allowed_methods())) {
    pv_abort("Method comparison diagnostics must use implemented method IDs.")
  }
  method_required <- c(
    "df_method", "interval_role", "coverage_claim_allowed",
    "n_descriptive_intervals", "target_source", "target_hash",
    "pooling_source", "pooling_hash", "center", "shared_target_hash",
    "shared_pooling_hash", "shared_external_target",
    "shares_reference_target", "shares_reference_pooling"
  )
  missing_method <- setdiff(method_required, names(method_diagnostics))
  if (length(missing_method) > 0L) {
    pv_abort(sprintf("Method comparison diagnostics are missing interval metadata column(s): %s.", paste(missing_method, collapse = ", ")))
  }
  if (!is.logical(method_diagnostics$coverage_claim_allowed)) {
    pv_abort("Method comparison method diagnostics `coverage_claim_allowed` must be logical.")
  }
  shared_fields <- c(
    "shared_target_hash", "shared_pooling_hash", "shared_external_target",
    "shares_reference_target", "shares_reference_pooling"
  )
  for (field in shared_fields) {
    if (!is.logical(method_diagnostics[[field]])) {
      pv_abort(sprintf("Method comparison method diagnostics `%s` must be logical.", field))
    }
  }
  if (!is.integer(method_diagnostics$n_descriptive_intervals) && !is.numeric(method_diagnostics$n_descriptive_intervals)) {
    pv_abort("Method comparison method diagnostics `n_descriptive_intervals` must be numeric.")
  }
  if (any(!is.finite(method_diagnostics$n_descriptive_intervals) |
      method_diagnostics$n_descriptive_intervals < 0 |
      method_diagnostics$n_descriptive_intervals != floor(method_diagnostics$n_descriptive_intervals))) {
    pv_abort("Method comparison method diagnostics `n_descriptive_intervals` must contain non-negative integer counts.")
  }
  elapsed <- timing$elapsed_seconds
  if (!is.numeric(elapsed) || any(!is.na(elapsed) & (!is.finite(elapsed) | elapsed < 0))) {
    pv_abort("Method comparison timing diagnostics must use finite non-negative elapsed seconds or NA.")
  }
  n_fits <- timing$n_fits
  if (!is.integer(n_fits) && !is.numeric(n_fits)) {
    pv_abort("Method comparison timing diagnostics `n_fits` must be numeric.")
  }
  if (any(!is.na(n_fits) & (!is.finite(n_fits) | n_fits < 1 | n_fits != floor(n_fits)))) {
    pv_abort("Method comparison timing diagnostics `n_fits` must contain positive integer counts or NA.")
  }
  statuses <- unname(as.character(comparison$diagnostics$statuses))
  if (!identical(statuses, unname(as.character(comparison$diagnostic_table$status)))) {
    pv_abort("Method comparison diagnostic statuses must align.")
  }
  blocked <- labels[statuses == "blocked"]
  warning <- labels[statuses == "warning"]
  if (!identical(as.character(comparison$diagnostics$blocked_methods), blocked) ||
      !identical(as.character(comparison$diagnostics$warning_methods), warning)) {
    pv_abort("Method comparison blocked and warning method labels must align with statuses.")
  }
  expected_shared_target <- pv_compare_shared_hash(method_diagnostics$target_hash)
  expected_shared_pooling <- pv_compare_shared_hash(method_diagnostics$pooling_hash)
  expected_shared_external <- pv_compare_shared_external_target(
    method_diagnostics$target_source,
    method_diagnostics$target_hash
  )
  expected_reference_target <- pv_compare_shares_reference_hash(
    method_diagnostics$target_hash,
    comparison$reference_method,
    labels
  )
  expected_reference_pooling <- pv_compare_shares_reference_hash(
    method_diagnostics$pooling_hash,
    comparison$reference_method,
    labels
  )
  if (!identical(as.logical(method_diagnostics$shared_target_hash), expected_shared_target) ||
      !identical(as.logical(method_diagnostics$shared_pooling_hash), expected_shared_pooling) ||
      !identical(as.logical(method_diagnostics$shared_external_target), expected_shared_external) ||
      !identical(as.logical(method_diagnostics$shares_reference_target), expected_reference_target) ||
      !identical(as.logical(method_diagnostics$shares_reference_pooling), expected_reference_pooling)) {
    pv_abort("Method comparison shared-provenance flags must match target and pooling hashes.")
  }
  target_overlap <- comparison$diagnostics$target_overlap
  pv_assert_named_list(target_overlap, "diagnostics$target_overlap")
  target_overlap_required <- c(
    "shared_external_target", "shared_target_hash", "shared_pooling_hash",
    "shares_reference_target", "shares_reference_pooling",
    "shared_target_sources", "independence_caveat_required",
    "independence_caveat"
  )
  missing_target_overlap <- setdiff(target_overlap_required, names(target_overlap))
  if (length(missing_target_overlap) > 0L) {
    pv_abort(sprintf("Method comparison target-overlap diagnostics are missing required field(s): %s.", paste(missing_target_overlap, collapse = ", ")))
  }
  expected_target_overlap <- pv_compare_target_overlap_diagnostics(method_diagnostics)
  if (!identical(target_overlap$shared_external_target, expected_target_overlap$shared_external_target) ||
      !identical(target_overlap$shared_target_hash, expected_target_overlap$shared_target_hash) ||
      !identical(target_overlap$shared_pooling_hash, expected_target_overlap$shared_pooling_hash) ||
      !identical(target_overlap$shares_reference_target, expected_target_overlap$shares_reference_target) ||
      !identical(target_overlap$shares_reference_pooling, expected_target_overlap$shares_reference_pooling) ||
      !identical(as.character(target_overlap$shared_target_sources), as.character(expected_target_overlap$shared_target_sources)) ||
      !identical(target_overlap$independence_caveat_required, expected_target_overlap$independence_caveat_required)) {
    pv_abort("Method comparison target-overlap diagnostics must summarize method diagnostics.")
  }
  pv_assert_scalar_string(target_overlap$independence_caveat, "diagnostics$target_overlap$independence_caveat")
  for (label in labels) {
    rows <- comparison$table[comparison$table$method_label == label, , drop = FALSE]
    row <- agreement[agreement$method_label == label, , drop = FALSE]
    method_row <- method_diagnostics[method_diagnostics$method_label == label, , drop = FALSE]
    if (nrow(row) != 1L) {
      pv_abort("Method comparison agreement diagnostics must contain one row per method label.")
    }
    if (nrow(method_row) != 1L) {
      pv_abort("Method comparison method diagnostics must contain one row per method label.")
    }
    finite_diff <- rows$abs_z_diff[is.finite(rows$abs_z_diff)]
    finite_ratio <- abs(log(rows$se_ratio[is.finite(rows$se_ratio) & rows$se_ratio > 0]))
    expected_z <- if (length(finite_diff) == 0L) NA_real_ else max(finite_diff)
    expected_ratio <- if (length(finite_ratio) == 0L) NA_real_ else max(finite_ratio)
    available <- rows[is.finite(rows$estimate), , drop = FALSE]
    expected_df_method <- if (nrow(available) == 0L) NA_character_ else paste(unique(as.character(available$df_method[!is.na(available$df_method)])), collapse = ",")
    expected_interval_role <- if (nrow(available) == 0L) NA_character_ else paste(unique(as.character(available$interval_role[!is.na(available$interval_role)])), collapse = ",")
    if (!nzchar(expected_df_method)) expected_df_method <- NA_character_
    if (!nzchar(expected_interval_role)) expected_interval_role <- NA_character_
    coverage_values <- available$coverage_claim_allowed[!is.na(available$coverage_claim_allowed)]
    expected_coverage <- if (length(coverage_values) == 0L) NA else all(coverage_values)
    expected_descriptive <- if (nrow(available) == 0L) 0L else sum(!is.na(available$coverage_claim_allowed) & !available$coverage_claim_allowed)
    if (!pv_numeric_equal_na(row$max_abs_z_diff, expected_z, tol = 1e-12) ||
        !pv_numeric_equal_na(row$max_abs_log_se_ratio, expected_ratio, tol = 1e-12) ||
        !identical(as.integer(row$n_terms), length(unique(rows$term))) ||
        !identical(as.integer(row$n_available), sum(is.finite(rows$estimate)))) {
      pv_abort("Method comparison agreement diagnostics must match the comparison table.")
    }
    if (!identical(as.character(method_row$df_method), expected_df_method) ||
        !identical(as.character(method_row$interval_role), expected_interval_role) ||
        !identical(as.logical(method_row$coverage_claim_allowed), as.logical(expected_coverage)) ||
        !identical(as.integer(method_row$n_descriptive_intervals), as.integer(expected_descriptive))) {
      pv_abort("Method comparison method diagnostics must summarize interval metadata from the comparison table.")
    }
  }
  invisible(comparison)
}

validate_pvstackr_method_comparison <- function(comparison) {
  pv_assert_named_list(comparison, "comparison")
  required <- c(
    "table", "estimate_table", "diagnostics", "diagnostic_table",
    "agreement", "timing", "fits", "reference_method", "methods",
    "method_labels", "created_at", "schema_version", "provenance", "warnings"
  )
  missing <- setdiff(required, names(comparison))
  if (length(missing) > 0L) {
    pv_abort(sprintf("Method comparison is missing required field(s): %s.", paste(missing, collapse = ", ")))
  }
  labels <- pv_validate_comparison_labels(comparison)
  if (!comparison$reference_method %in% comparison$method_labels) {
    pv_abort("Method comparison reference method must be one of the method labels.")
  }
  pv_validate_comparison_table(comparison, labels)
  pv_validate_comparison_diagnostics(comparison, labels)
  if (!is.null(comparison$fits)) {
    pv_assert_named_list(comparison$fits, "fits")
    if (!identical(names(comparison$fits), labels)) {
      pv_abort("Retained comparison fits must align with method labels.")
    }
    for (label in labels) {
      fit <- comparison$fits[[label]]
      if (!inherits(fit, "pvstackr_fit")) {
        pv_abort("Retained comparison fits must be `pvstackr_fit` objects.")
      }
      validate_pvstackr_fit(fit)
      if (!identical(fit$method, unname(comparison$methods[[label]]))) {
        pv_abort("Retained comparison fit methods must align with method labels.")
      }
    }
  }
  pv_assert_scalar_string(comparison$created_at, "created_at")
  pv_validate_schema_version(comparison$schema_version)
  pv_validate_named_list_field(comparison$provenance, "provenance")
  pv_validate_character_field(comparison$warnings, "warnings")
  invisible(comparison)
}

#' Compare pvstackr Method Fits
#'
#' `pv_compare_methods()` compares already-created `pvstackr_fit` objects by
#' aligning their fixed-effect estimate tables, computing differences against a
#' reference method, and recording agreement and timing diagnostics.
#'
#' The returned object keeps blocked methods visible. Blocked fits receive rows
#' for every aligned fixed-effect term, but their reportable numeric comparison
#' fields are `NA` and their reason codes remain available.
#'
#' The comparison estimate table preserves interval metadata from each fit:
#' `df_method`, `df_complete`, `conf_level`, `interval_role`, and
#' `coverage_claim_allowed`. Default print methods keep their tables compact and
#' emit a one-line note when one or more compared intervals are descriptive
#' rather than coverage-claimable. In the method-level diagnostic table,
#' `coverage_claim_allowed` is `TRUE` only when all available reportable rows
#' for that method are coverage-claimable, `FALSE` when any available row is
#' descriptive, and `NA` when the method has no available estimate rows.
#'
#' Method diagnostics also carry target and pooling provenance where available:
#' `target_source`, `target_hash`, `pooling_source`, `pooling_hash`, and shared
#' hash flags. The nested `diagnostics$target_overlap` summary records whether
#' any compared methods share external targets, target hashes, pooling hashes,
#' reference-target hashes, or target-source families. Agreement diagnostics are
#' descriptive; close agreement is not automatic independent corroboration when
#' compared methods share a target hash, pooling hash, target source, or estimand
#' construction. `target_source` is provenance vocabulary, not necessarily a
#' formal target object: `stack_psis` rows use `stack_psis_rubin_pooling` even
#' though [get_target()] returns `NULL` for `stack_psis` fits.
#'
#' @param ... Two or more `pvstackr_fit` objects, or a single list of fits.
#' @param fits Optional explicit list of `pvstackr_fit` objects.
#' @param reference_method Optional method label or method id used as the
#'   comparison reference. Defaults to the first non-blocked `per_pv` fit when
#'   available, otherwise the first non-blocked fit.
#' @param timings Optional named numeric vector of elapsed seconds, named by
#'   method label or method id.
#' @param include_fits Whether to retain the original fit objects in the
#'   comparison object.
#'
#' @returns A `pvstackr_method_comparison` object with top-level
#'   `estimate_table`, `diagnostic_table`, `agreement`, and `timing` fields.
#'   `estimate_table` includes method-level interval metadata alongside the
#'   aligned fixed-effect estimates. `diagnostic_table` includes method-level
#'   target/pooling provenance and shared-provenance flags.
#'
#' @section Interpreting agreement:
#' Agreement diagnostics are **descriptive**, not independent corroboration.
#' Close agreement does **not** confirm a result when the compared methods share
#' a target hash, pooling hash, target source, or estimand construction: the
#' methods are then reading the same information, so concordance is expected and
#' uninformative about correctness. The nested `diagnostics$target_overlap`
#' summary flags this, and the print methods emit a one-line provenance note when
#' overlap is detected. Only `stack_direct` rows backed by the external
#' Rubin/BRR-Fay target are coverage-claimable; `per_pv` and `stack_psis` rows
#' are descriptive/reference even with Barnard-Rubin degrees of freedom.
#'
#' @examples
#' path <- system.file(
#'   "extdata", "examples", "pisa_tiny_stack_direct.rds", package = "pvstackr"
#' )
#' if (nzchar(path)) {
#'   fit_direct <- readRDS(path)$fit          # cached stack_direct pvstackr_fit
#'
#'   # Build a per_pv reference fit sharing the same fixed-effect names
#'   # (b_Intercept, b_x, b_female) from small injected posterior draw clouds.
#'   set.seed(1)
#'   fe_names <- c("b_Intercept", "b_x", "b_female")
#'   draw_block <- function() {
#'     matrix(rnorm(200 * 3), ncol = 3, dimnames = list(NULL, fe_names))
#'   }
#'   fit_per_pv <- pv_fit_reference(
#'     per_pv_draws = list(PV1READ = draw_block(), PV2READ = draw_block()),
#'     control      = pv_control(method = "per_pv")
#'   )
#'
#'   cmp <- pv_compare_methods(stack_direct = fit_direct, per_pv = fit_per_pv)
#'   cmp                          # descriptive agreement + provenance notes
#'   head(get_estimates(cmp))     # aligned per-method estimate rows
#' }
#' @family pvstackr-comparison
#' @seealso [pv_fit()], [pv_fit_direct()], [pv_fit_reference()],
#'   [pv_fit_stack_psis()], [get_estimates()], [get_diagnostics()]
#' @export
pv_compare_methods <- function(
  ...,
  fits = NULL,
  reference_method = NULL,
  timings = NULL,
  include_fits = FALSE
) {
  fits <- pv_compare_collect_fits(..., fits = fits)
  timings <- pv_compare_normalize_timings(timings)
  reference_label <- pv_compare_reference_label(fits, reference_method = reference_method)
  terms <- pv_compare_terms(fits)
  table <- pv_compare_table(fits, terms, reference_label)
  diagnostics <- pv_compare_diagnostics(table, fits, reference_label, timings)
  new_pvstackr_method_comparison(
    table = table,
    diagnostics = diagnostics,
    fits = fits,
    reference_method = reference_label,
    include_fits = include_fits
  )
}
