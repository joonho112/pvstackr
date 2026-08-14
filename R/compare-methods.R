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
    validate_pvstackr_fit(dots[[i]], tier = "deep")
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
  if (inherits(fit, "pvstackr_comparison_source_projection")) {
    return(as.integer(fit$n_fits))
  }
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
  if (inherits(fit, "pvstackr_comparison_source_projection")) {
    return(fit$target_hash)
  }
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
  if (inherits(fit, "pvstackr_comparison_source_projection")) {
    return(fit$target_source)
  }
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
  if (inherits(fit, "pvstackr_comparison_source_projection")) {
    return(fit$pooling_hash)
  }
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
  if (inherits(fit, "pvstackr_comparison_source_projection")) {
    return(fit$pooling_source)
  }
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
  if (inherits(fit, "pvstackr_comparison_source_projection")) {
    return(fit$center)
  }
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
  if (inherits(fit, "pvstackr_comparison_source_projection")) {
    return(fit$psis_status)
  }
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
  if (inherits(fit, "pvstackr_comparison_source_projection")) {
    return(fit$pareto_k_max)
  }
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

pv_compare_fit_weight_method <- function(fit) {
  if (inherits(fit, "pvstackr_comparison_source_projection")) {
    return(fit$weight_method)
  }
  if (identical(fit$method, "stack_psis") &&
      is.list(fit$diagnostics$psis) &&
      is.character(fit$diagnostics$psis$weight_method) &&
      length(fit$diagnostics$psis$weight_method) == 1L) {
    return(fit$diagnostics$psis$weight_method)
  }
  NA_character_
}

pv_compare_fit_psis_provenance <- function(fit, field) {
  if (inherits(fit, "pvstackr_comparison_source_projection")) {
    return(fit[[field]])
  }
  if (!identical(fit$method, "stack_psis") ||
      !is.list(fit$diagnostics$psis)) {
    return(NA_character_)
  }
  diagnostic_field <- switch(
    field,
    psis_source = "source",
    psis_producer = "producer",
    psis_producer_version = "producer_version",
    field
  )
  value <- fit$diagnostics$psis[[diagnostic_field]]
  if (!is.character(value) || length(value) != 1L) NA_character_ else value
}

pv_compare_fit_weight_stat <- function(
  fit,
  field,
  aggregate = c("min", "max")
) {
  aggregate <- match.arg(aggregate)
  if (inherits(fit, "pvstackr_comparison_source_projection")) {
    projection_field <- switch(
      field,
      weight_ess_iid = "weight_ess_iid_min",
      weight_ess_fraction = "weight_ess_fraction_min",
      max_normalized_weight = "max_normalized_weight_max",
      field
    )
    return(fit[[projection_field]])
  }
  if (!identical(fit$method, "stack_psis") ||
      !is.list(fit$diagnostics$psis)) {
    return(NA_real_)
  }
  value <- fit$diagnostics$psis[[field]]
  if (!is.numeric(value) || length(value) == 0L || any(!is.finite(value))) {
    return(NA_real_)
  }
  if (identical(aggregate, "min")) min(value) else max(value)
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
  psis_source <- pv_compare_fit_psis_provenance(fit, "psis_source")
  pareto_k_source <- pv_compare_fit_psis_provenance(fit, "pareto_k_source")
  weight_method <- pv_compare_fit_weight_method(fit)
  psis_producer <- pv_compare_fit_psis_provenance(fit, "psis_producer")
  psis_producer_version <- pv_compare_fit_psis_provenance(
    fit,
    "psis_producer_version"
  )
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
    psis_source = psis_source,
    pareto_k_source = pareto_k_source,
    weight_method = weight_method,
    psis_producer = psis_producer,
    psis_producer_version = psis_producer_version,
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
    psis_source = vapply(
      fits,
      pv_compare_fit_psis_provenance,
      character(1),
      field = "psis_source"
    ),
    pareto_k_source = vapply(
      fits,
      pv_compare_fit_psis_provenance,
      character(1),
      field = "pareto_k_source"
    ),
    weight_method = vapply(fits, pv_compare_fit_weight_method, character(1)),
    psis_producer = vapply(
      fits,
      pv_compare_fit_psis_provenance,
      character(1),
      field = "psis_producer"
    ),
    psis_producer_version = vapply(
      fits,
      pv_compare_fit_psis_provenance,
      character(1),
      field = "psis_producer_version"
    ),
    weight_ess_iid_min = vapply(
      fits,
      pv_compare_fit_weight_stat,
      numeric(1),
      field = "weight_ess_iid",
      aggregate = "min"
    ),
    weight_ess_fraction_min = vapply(
      fits,
      pv_compare_fit_weight_stat,
      numeric(1),
      field = "weight_ess_fraction",
      aggregate = "min"
    ),
    max_normalized_weight_max = vapply(
      fits,
      pv_compare_fit_weight_stat,
      numeric(1),
      field = "max_normalized_weight",
      aggregate = "max"
    ),
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

pv_comparison_schema_version <- function() {
  "0.2.0"
}

pv_comparison_validation_schema <- function() {
  "pvstackr_method_comparison_validation_v1"
}

pv_comparison_validation_sentinel <- function() {
  paste0("sha256:", strrep("0", 64L))
}

pv_comparison_source_fit_validation <- function(fits) {
  out <- lapply(fits, function(fit) {
    validate_pvstackr_fit(fit, tier = "deep")
    list(
      schema_version = fit$validation$schema_version,
      stamp = fit$validation$stamp
    )
  })
  names(out) <- names(fits)
  out
}

pv_comparison_source_fit_reportability <- function(fits) {
  out <- lapply(fits, function(fit) {
    if (identical(fit$method, "stack_psis")) {
      pv_fit_summary_reportability_fit(fit)
    } else {
      NULL
    }
  })
  names(out) <- names(fits)
  out
}

pv_comparison_authority_fits <- function(fits, reportability) {
  out <- fits
  for (label in names(out)) {
    if (!is.null(reportability[[label]])) {
      out[[label]] <- reportability[[label]]
    }
  }
  out
}

pv_comparison_validation_record <- function(
  stamp = pv_comparison_validation_sentinel()
) {
  list(
    schema_version = pv_comparison_validation_schema(),
    policy_id = "source_fit_stamps_plus_owned_payload_sha256_v1",
    canonicalizer_id = "r_xdr_v2_comparison_owned_payload_v1",
    stamp = stamp
  )
}

pv_comparison_validation_projection <- function(comparison) {
  projected <- comparison
  projected$validation$stamp <- pv_comparison_validation_sentinel()
  if (!is.null(projected$fits)) {
    projected$fits <- projected$source_fit_validation
  }
  projected
}

pv_comparison_validation_digest <- function(comparison) {
  payload <- serialize(
    pv_comparison_validation_projection(comparison),
    NULL,
    ascii = FALSE,
    xdr = TRUE,
    version = 2L
  )
  bytes <- c(
    charToRaw("pvstackr-method-comparison-validation-v1"),
    as.raw(0L),
    payload
  )
  paste0("sha256:", digest::digest(bytes, algo = "sha256", serialize = FALSE))
}

pv_comparison_issue_validation_stamp <- function(comparison) {
  comparison$validation <- pv_comparison_validation_record()
  comparison$validation$stamp <- pv_comparison_validation_digest(comparison)
  comparison
}

pv_comparison_contains_stack_psis <- function(comparison) {
  values <- character()
  if (is.character(comparison$methods)) {
    values <- c(values, unname(comparison$methods))
  }
  if (is.data.frame(comparison$table) && "method" %in% names(comparison$table)) {
    values <- c(values, as.character(comparison$table$method))
  }
  if (is.data.frame(comparison$diagnostic_table) &&
      "method" %in% names(comparison$diagnostic_table)) {
    values <- c(values, as.character(comparison$diagnostic_table$method))
  }
  any(values == "stack_psis", na.rm = TRUE)
}

pv_validate_comparison_source_fit_validation <- function(comparison, labels) {
  source <- comparison$source_fit_validation
  if (!is.list(source) || !identical(names(source), labels) ||
      !identical(attributes(source), list(names = labels))) {
    pv_abort("Current method comparison source-fit validation records must align with method labels.")
  }
  for (label in labels) {
    record <- source[[label]]
    fields <- c("schema_version", "stamp")
    if (!is.list(record) || !identical(names(record), fields) ||
        !identical(attributes(record), list(names = fields)) ||
        !identical(record$schema_version, pv_fit_validation_schema()) ||
        !is.character(record$stamp) || length(record$stamp) != 1L ||
        is.na(record$stamp) ||
        !grepl("^sha256:[0-9a-f]{64}$", record$stamp)) {
      pv_abort("Current method comparison source-fit validation records are noncanonical.")
    }
  }
  invisible(source)
}

pv_validate_comparison_validation <- function(comparison, labels) {
  validation <- comparison$validation
  fields <- c("schema_version", "policy_id", "canonicalizer_id", "stamp")
  if (!is.list(validation) || !identical(names(validation), fields) ||
      !identical(attributes(validation), list(names = fields)) ||
      !identical(
        validation,
        pv_comparison_validation_record(validation$stamp %||% "")
      ) || !grepl("^sha256:[0-9a-f]{64}$", validation$stamp)) {
    pv_abort("Current method comparison validation record is noncanonical.")
  }
  pv_validate_comparison_source_fit_validation(comparison, labels)
  if (!identical(comparison$schema_version, pv_comparison_schema_version())) {
    pv_abort("Current method comparison schema version is unsupported.")
  }
  if (!identical(
        comparison$provenance,
        pv_provenance(
          "pv_compare_methods",
          schema_version = pv_comparison_schema_version()
        )
      )) {
    pv_abort("Current method comparison provenance must match its schema authority.")
  }
  if (!identical(validation$stamp, pv_comparison_validation_digest(comparison))) {
    pv_abort("Method comparison validation stamp does not match the current owned payload.")
  }
  invisible(validation)
}

pv_validate_stack_psis_derived_tables <- function(
  estimate_table,
  diagnostic_table,
  context = "method comparison"
) {
  if (!is.data.frame(diagnostic_table) ||
      !all(c(
        "method", "method_label", "status", "reason_codes", "n_available",
        "target_source", "target_hash", "pooling_source", "pooling_hash",
        "df_method", "interval_role", "coverage_claim_allowed",
        "n_descriptive_intervals", "psis_status", "pareto_k_max",
        "psis_source", "pareto_k_source", "weight_method", "psis_producer",
        "psis_producer_version", "weight_ess_iid_min",
        "weight_ess_fraction_min", "max_normalized_weight_max"
      ) %in% names(diagnostic_table))) {
    pv_abort(sprintf("%s lacks the fields required by the PSIS derived-output firewall.", context))
  }
  psis_rows <- diagnostic_table[diagnostic_table$method == "stack_psis", , drop = FALSE]
  if (nrow(psis_rows) == 0L) {
    return(invisible(TRUE))
  }
  if (any(psis_rows$status == "warning")) {
    pv_abort(sprintf("Legacy warning-status stack_psis %s output is inspection-only.", context))
  }
  numeric_result_fields <- intersect(
    c(
      "estimate", "se", "std.error", "df", "df_complete", "conf_level",
      "conf_low", "conf_high", "estimate_diff", "se_ratio", "abs_z_diff"
    ),
    names(estimate_table)
  )
  for (i in seq_len(nrow(psis_rows))) {
    row <- psis_rows[i, , drop = FALSE]
    label <- as.character(row$method_label[[1L]])
    status <- as.character(row$status[[1L]])
    table_rows <- estimate_table[
      as.character(estimate_table$method_label) == label,
      ,
      drop = FALSE
    ]
    weight_method <- as.character(row$weight_method[[1L]])
    psis_source <- as.character(row$psis_source[[1L]])
    pareto_k_source <- as.character(row$pareto_k_source[[1L]])
    producer <- as.character(row$psis_producer[[1L]])
    producer_version <- as.character(row$psis_producer_version[[1L]])
    weight_values <- c(
      row$weight_ess_iid_min[[1L]],
      row$weight_ess_fraction_min[[1L]],
      row$max_normalized_weight_max[[1L]]
    )
    if (!psis_source %in% c(
      "supplied_psis_weights", "injected_psis_function",
      "self_normalized_log_ratios"
    ) || !pareto_k_source %in% c("supplied", "injected_function_output") ||
        !weight_method %in% c(
      "caller_declared_external_psis", "unspecified_external",
      "self_normalized_raw_importance"
    ) || any(!is.finite(weight_values)) || weight_values[[1L]] < 1 ||
        weight_values[[2L]] <= 0 || weight_values[[2L]] > 1 ||
        weight_values[[3L]] <= 0 || weight_values[[3L]] > 1) {
      pv_abort(sprintf("stack_psis %s weight provenance and concentration diagnostics are incoherent.", context))
    }
    if (identical(status, "ok")) {
      if (!identical(as.character(row$psis_status[[1L]]), "ok") ||
          !identical(weight_method, "caller_declared_external_psis") ||
          is.na(producer) || !nzchar(producer) ||
          is.na(producer_version) || !nzchar(producer_version) ||
          !is.finite(row$pareto_k_max[[1L]]) ||
          row$pareto_k_max[[1L]] >= pv_legacy_psis_hard_threshold() ||
          as.integer(row$n_available[[1L]]) < 1L) {
        pv_abort(sprintf("Reportable stack_psis %s output fails the immutable Pareto-k firewall.", context))
      }
      next
    }
    if (!identical(status, "blocked") ||
        !as.character(row$psis_status[[1L]]) %in% c(
          "failed", "not_evaluated", "unsmoothed", "provenance_incomplete"
        ) ||
        as.integer(row$n_available[[1L]]) != 0L ||
        !nzchar(as.character(row$reason_codes[[1L]]))) {
      pv_abort(sprintf("Non-reportable stack_psis %s status metadata are incoherent.", context))
    }
    absent_metadata <- c(
      row$target_source[[1L]], row$target_hash[[1L]],
      row$pooling_source[[1L]], row$pooling_hash[[1L]],
      row$df_method[[1L]], row$interval_role[[1L]]
    )
    if (any(!is.na(absent_metadata) & nzchar(as.character(absent_metadata))) ||
        !is.na(row$coverage_claim_allowed[[1L]]) ||
        as.integer(row$n_descriptive_intervals[[1L]]) != 0L ||
        any(vapply(table_rows[numeric_result_fields], function(value) {
          any(!is.na(value))
        }, logical(1)))) {
      pv_abort(sprintf("Blocked stack_psis %s output retains reportable or pooling metadata.", context))
    }
  }
  invisible(TRUE)
}

pv_comparison_source_projection_validation_schema <- function() {
  "pvstackr_comparison_source_projection_validation_v1"
}

pv_comparison_source_projection_validation_record <- function(
  stamp = pv_comparison_validation_sentinel()
) {
  list(
    schema_version = pv_comparison_source_projection_validation_schema(),
    policy_id = "validated_fit_reportability_projection_sha256_v1",
    stamp = stamp
  )
}

pv_comparison_source_projection_digest <- function(projection) {
  projected <- projection
  projected$validation$stamp <- pv_comparison_validation_sentinel()
  payload <- serialize(
    projected,
    NULL,
    ascii = FALSE,
    xdr = TRUE,
    version = 2L
  )
  bytes <- c(
    charToRaw("pvstackr-comparison-source-projection-v1"),
    as.raw(0L),
    payload
  )
  paste0("sha256:", digest::digest(bytes, algo = "sha256", serialize = FALSE))
}

pv_comparison_source_projection <- function(fit, source_validation = NULL) {
  validate_pvstackr_fit(fit, tier = "deep")
  source_validation <- source_validation %||% list(
    schema_version = fit$validation$schema_version,
    stamp = fit$validation$stamp
  )
  out <- list(
    method = fit$method,
    status = fit$status,
    reason_codes = fit$reason_codes,
    warnings = fit$warnings,
    estimates = fit$estimates,
    n_fits = pv_compare_fit_n(fit),
    target_source = pv_compare_fit_target_source(fit),
    target_hash = pv_compare_fit_target_hash(fit),
    pooling_source = pv_compare_fit_pooling_source(fit),
    pooling_hash = pv_compare_fit_pooling_hash(fit),
    center = pv_compare_fit_center(fit),
    psis_status = pv_compare_fit_psis_status(fit),
    pareto_k_max = pv_compare_fit_pareto_k_max(fit),
    psis_source = pv_compare_fit_psis_provenance(fit, "psis_source"),
    pareto_k_source = pv_compare_fit_psis_provenance(fit, "pareto_k_source"),
    weight_method = pv_compare_fit_weight_method(fit),
    psis_producer = pv_compare_fit_psis_provenance(fit, "psis_producer"),
    psis_producer_version = pv_compare_fit_psis_provenance(
      fit,
      "psis_producer_version"
    ),
    weight_ess_iid_min = pv_compare_fit_weight_stat(
      fit,
      "weight_ess_iid",
      "min"
    ),
    weight_ess_fraction_min = pv_compare_fit_weight_stat(
      fit,
      "weight_ess_fraction",
      "min"
    ),
    max_normalized_weight_max = pv_compare_fit_weight_stat(
      fit,
      "max_normalized_weight",
      "max"
    ),
    source_validation = source_validation,
    validation = NULL
  )
  class(out) <- c("pvstackr_comparison_source_projection", "list")
  out$validation <- pv_comparison_source_projection_validation_record()
  out$validation$stamp <- pv_comparison_source_projection_digest(out)
  pv_validate_comparison_source_projection(out)
  out
}

pv_comparison_source_fit_projection <- function(fits, source_validation) {
  out <- lapply(names(fits), function(label) {
    pv_comparison_source_projection(
      fits[[label]],
      source_validation[[label]]
    )
  })
  names(out) <- names(fits)
  out
}

pv_validate_comparison_source_projection <- function(projection) {
  fields <- c(
    "method", "status", "reason_codes", "warnings", "estimates", "n_fits",
    "target_source", "target_hash", "pooling_source", "pooling_hash",
    "center", "psis_status", "pareto_k_max", "psis_source",
    "pareto_k_source", "weight_method", "psis_producer",
    "psis_producer_version",
    "weight_ess_iid_min", "weight_ess_fraction_min",
    "max_normalized_weight_max", "source_validation",
    "validation"
  )
  attrs <- attributes(projection)
  if (!is.list(projection) || !identical(names(projection), fields) ||
      !identical(names(attrs), c("names", "class")) ||
      !identical(attrs$names, fields) ||
      !identical(
        attrs$class,
        c("pvstackr_comparison_source_projection", "list")
      )) {
    pv_abort("Comparison source projection fields, order, class, and attributes must be exact.")
  }
  pv_validate_method(projection$method)
  if (!is.character(projection$status) || length(projection$status) != 1L ||
      is.na(projection$status) ||
      !projection$status %in% c("ok", "warning", "blocked") ||
      !is.data.frame(projection$estimates) ||
      !is.numeric(projection$n_fits) || length(projection$n_fits) != 1L ||
      (!is.na(projection$n_fits) && projection$n_fits < 1L) ||
      !is.character(projection$reason_codes) ||
      !is.character(projection$warnings) ||
      !is.numeric(projection$pareto_k_max) ||
      length(projection$pareto_k_max) != 1L ||
      !is.character(projection$psis_source) ||
      length(projection$psis_source) != 1L ||
      !is.character(projection$pareto_k_source) ||
      length(projection$pareto_k_source) != 1L ||
      !is.character(projection$weight_method) ||
      length(projection$weight_method) != 1L ||
      !is.character(projection$psis_producer) ||
      length(projection$psis_producer) != 1L ||
      !is.character(projection$psis_producer_version) ||
      length(projection$psis_producer_version) != 1L ||
      !is.numeric(projection$weight_ess_iid_min) ||
      length(projection$weight_ess_iid_min) != 1L ||
      !is.numeric(projection$weight_ess_fraction_min) ||
      length(projection$weight_ess_fraction_min) != 1L ||
      !is.numeric(projection$max_normalized_weight_max) ||
      length(projection$max_normalized_weight_max) != 1L ||
      any(!vapply(
        projection[c(
          "method", "status", "reason_codes", "warnings", "n_fits",
          "target_source", "target_hash", "pooling_source", "pooling_hash",
          "center", "psis_status", "pareto_k_max", "psis_source",
          "pareto_k_source", "weight_method", "psis_producer",
          "psis_producer_version",
          "weight_ess_iid_min", "weight_ess_fraction_min",
          "max_normalized_weight_max"
        )],
        function(value) is.null(attributes(value)),
        logical(1)
      ))) {
    pv_abort("Comparison source projection scalar and estimate fields are noncanonical.")
  }
  pv_validate_summary_source(
    projection$source_validation,
    pv_fit_validation_schema()
  )
  validation <- projection$validation
  validation_fields <- c("schema_version", "policy_id", "stamp")
  if (!is.list(validation) ||
      !identical(names(validation), validation_fields) ||
      !identical(attributes(validation), list(names = validation_fields)) ||
      !identical(
        validation,
        pv_comparison_source_projection_validation_record(
          validation$stamp %||% ""
        )
      ) || !grepl("^sha256:[0-9a-f]{64}$", validation$stamp) ||
      !identical(
        validation$stamp,
        pv_comparison_source_projection_digest(projection)
      )) {
    pv_abort("Comparison source projection validation stamp is invalid.")
  }
  if (identical(projection$method, "stack_psis")) {
    if (identical(projection$status, "warning")) {
      pv_abort("Warning-status PSIS source projections are inspection-only.")
    } else if (identical(projection$status, "ok")) {
      if (!identical(projection$psis_status, "ok") ||
          !identical(
            projection$weight_method,
            "caller_declared_external_psis"
          ) ||
          is.na(projection$psis_producer) ||
          !nzchar(projection$psis_producer) ||
          is.na(projection$psis_producer_version) ||
          !nzchar(projection$psis_producer_version) ||
          !is.finite(projection$pareto_k_max) ||
          projection$pareto_k_max >= pv_legacy_psis_hard_threshold() ||
          nrow(projection$estimates) < 1L) {
        pv_abort("Reportable PSIS source projection fails the immutable gate.")
      }
    } else if (!projection$psis_status %in% c(
      "failed", "not_evaluated", "unsmoothed", "provenance_incomplete"
    ) ||
        nrow(projection$estimates) != 0L ||
        any(!is.na(c(
          projection$target_source, projection$target_hash,
          projection$pooling_source, projection$pooling_hash
        ))) || length(projection$reason_codes) == 0L) {
      pv_abort("Blocked PSIS source projection retains reportable metadata.")
    }
    weight_values <- c(
      projection$weight_ess_iid_min,
      projection$weight_ess_fraction_min,
      projection$max_normalized_weight_max
    )
    if (!projection$weight_method %in% c(
      "caller_declared_external_psis", "unspecified_external",
      "self_normalized_raw_importance"
    ) || any(!is.finite(weight_values)) || weight_values[[1L]] < 1 ||
        weight_values[[2L]] <= 0 || weight_values[[2L]] > 1 ||
        weight_values[[3L]] <= 0 || weight_values[[3L]] > 1) {
      pv_abort("PSIS source projection weight diagnostics are incoherent.")
    }
  }
  invisible(projection)
}

pv_validate_comparison_derivation <- function(comparison, labels) {
  projections <- comparison$source_fit_projection
  if (!is.list(projections) || !identical(names(projections), labels) ||
      !identical(attributes(projections), list(names = labels))) {
    pv_abort("Comparison source projections must align with method labels.")
  }
  for (label in labels) {
    projection <- projections[[label]]
    pv_validate_comparison_source_projection(projection)
    if (!identical(
          projection$source_validation,
          comparison$source_fit_validation[[label]]
        ) || !identical(
          projection$method,
          unname(comparison$methods[[label]])
        )) {
      pv_abort("Comparison source projection does not match its source-fit authority.")
    }
  }
  expected_terms <- pv_compare_terms(projections)
  expected_table <- pv_compare_table(
    projections,
    expected_terms,
    comparison$reference_method
  )
  timing_values <- stats::setNames(
    as.numeric(comparison$timing$elapsed_seconds),
    labels
  )
  expected_diagnostics <- pv_compare_diagnostics(
    expected_table,
    projections,
    comparison$reference_method,
    timing_values
  )
  if (!identical(comparison$table, expected_table) ||
      !identical(comparison$estimate_table, expected_table) ||
      !identical(comparison$diagnostics, expected_diagnostics) ||
      !identical(
        comparison$diagnostic_table,
        expected_diagnostics$method_diagnostics
      ) || !identical(comparison$agreement, expected_diagnostics$agreement) ||
      !identical(comparison$timing, expected_diagnostics$timing)) {
    pv_abort("Method comparison must exactly reproduce its validated source projections.")
  }
  invisible(projections)
}

pv_validate_comparison_source_reportability <- function(comparison, labels) {
  reportability <- comparison$source_fit_reportability
  if (!is.list(reportability) || !identical(names(reportability), labels) ||
      !identical(attributes(reportability), list(names = labels))) {
    pv_abort("Comparison source reportability fits must align with method labels.")
  }
  for (label in labels) {
    method <- unname(comparison$methods[[label]])
    source <- reportability[[label]]
    if (identical(method, "stack_psis")) {
      if (!inherits(source, "pvstackr_fit")) {
        pv_abort("Every comparison PSIS source requires a validated reportability fit.")
      }
      validate_pvstackr_fit(source, tier = "deep")
      if (!identical(source$method, "stack_psis") ||
          !identical(
            comparison$source_fit_validation[[label]],
            list(
              schema_version = source$validation$schema_version,
              stamp = source$validation$stamp
            )
          ) || !identical(
            comparison$source_fit_projection[[label]],
            pv_comparison_source_projection(
              source,
              comparison$source_fit_validation[[label]]
            )
          )) {
        pv_abort("Comparison PSIS source authority and reportability projection disagree.")
      }
    } else if (!is.null(source)) {
      pv_abort("Only stack_psis comparison sources may retain a reportability fit.")
    }
  }
  invisible(reportability)
}

new_pvstackr_method_comparison <- function(table, diagnostics, fits, reference_method, include_fits = FALSE) {
  source_fit_reportability <- pv_comparison_source_fit_reportability(fits)
  authority_fits <- pv_comparison_authority_fits(
    fits,
    source_fit_reportability
  )
  source_fit_validation <- pv_comparison_source_fit_validation(authority_fits)
  source_fit_projection <- pv_comparison_source_fit_projection(
    authority_fits,
    source_fit_validation
  )
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
    source_fit_validation = source_fit_validation,
    source_fit_reportability = source_fit_reportability,
    source_fit_projection = source_fit_projection,
    created_at = as.character(Sys.time()),
    schema_version = pv_comparison_schema_version(),
    provenance = pv_provenance(
      "pv_compare_methods",
      schema_version = pv_comparison_schema_version()
    ),
    validation = NULL,
    warnings = character()
  )
  class(comparison) <- c("pvstackr_method_comparison", "list")
  comparison <- pv_comparison_issue_validation_stamp(comparison)
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
    "interval_role", "coverage_claim_allowed", "psis_source",
    "pareto_k_source", "weight_method", "psis_producer",
    "psis_producer_version", "conf_low", "conf_high",
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
    "shares_reference_target", "shares_reference_pooling", "psis_status",
    "pareto_k_max", "psis_source", "pareto_k_source", "weight_method",
    "psis_producer", "psis_producer_version", "weight_ess_iid_min",
    "weight_ess_fraction_min", "max_normalized_weight_max"
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
  legacy_required <- c(
    "table", "estimate_table", "diagnostics", "diagnostic_table",
    "agreement", "timing", "fits", "reference_method", "methods",
    "method_labels", "created_at", "schema_version", "provenance", "warnings"
  )
  current_required <- c(
    "table", "estimate_table", "diagnostics", "diagnostic_table",
    "agreement", "timing", "fits", "reference_method", "methods",
    "method_labels", "source_fit_validation", "source_fit_reportability",
    "source_fit_projection", "created_at", "schema_version", "provenance",
    "validation", "warnings"
  )
  is_current <- all(c("source_fit_validation", "validation") %in% names(comparison))
  required <- if (is_current) current_required else legacy_required
  missing <- setdiff(required, names(comparison))
  if (length(missing) > 0L) {
    pv_abort(sprintf("Method comparison is missing required field(s): %s.", paste(missing, collapse = ", ")))
  }
  if (!is_current && pv_comparison_contains_stack_psis(comparison)) {
    pv_abort(
      paste(
        "Legacy method comparisons containing stack_psis are inspection-only",
        "and must be rebuilt from current validated fits."
      )
    )
  }
  if (is_current) {
    root_attrs <- attributes(comparison)
    if (!identical(names(comparison), current_required) ||
        !identical(names(root_attrs), c("names", "class")) ||
        !identical(root_attrs$names, current_required) ||
        !identical(
          root_attrs$class,
          c("pvstackr_method_comparison", "list")
        )) {
      pv_abort("Current method comparison fields, order, class, and root attributes must be exact.")
    }
  }
  labels <- pv_validate_comparison_labels(comparison)
  if (!comparison$reference_method %in% comparison$method_labels) {
    pv_abort("Method comparison reference method must be one of the method labels.")
  }
  pv_validate_comparison_table(comparison, labels)
  pv_validate_comparison_diagnostics(comparison, labels)
  pv_validate_stack_psis_derived_tables(
    comparison$estimate_table,
    comparison$diagnostic_table
  )
  if (is_current) {
    pv_validate_comparison_source_reportability(comparison, labels)
    pv_validate_comparison_derivation(comparison, labels)
    pv_validate_comparison_validation(comparison, labels)
  }
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
      validate_pvstackr_fit(fit, tier = "cheap")
      if (is_current) {
        authority_fit <- if (identical(fit$method, "stack_psis")) {
          pv_fit_summary_reportability_fit(fit)
        } else {
          fit
        }
        expected_reportability <- if (identical(fit$method, "stack_psis")) {
          authority_fit
        } else {
          NULL
        }
        expected_validation <- list(
          schema_version = authority_fit$validation$schema_version,
          stamp = authority_fit$validation$stamp
        )
        if (!identical(
              comparison$source_fit_reportability[[label]],
              expected_reportability
            ) || !identical(
              comparison$source_fit_validation[[label]],
              expected_validation
            ) || !identical(
              comparison$source_fit_projection[[label]],
              pv_comparison_source_projection(
                authority_fit,
                expected_validation
              )
            )) {
          pv_abort("Retained comparison fits must reproduce their source reportability authorities.")
        }
      }
      if (!identical(fit$method, unname(comparison$methods[[label]]))) {
        pv_abort("Retained comparison fit methods must align with method labels.")
      }
    }
  }
  pv_assert_scalar_string(comparison$created_at, "created_at")
  if (!is_current) {
    pv_validate_schema_version(comparison$schema_version)
  }
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
#' For `stack_psis`, the aligned estimate rows additionally retain
#' `psis_source`, `pareto_k_source`, `weight_method`, `psis_producer`, and
#' `psis_producer_version`. The method diagnostic table adds the minimum
#' per-PV Kish-style iid weight ESS and ESS fraction, plus the maximum
#' normalized weight. These fields preserve the source fit's bounded
#' provenance and concentration diagnostics; they are not MCMC ESS or a
#' substitute for Pareto-k gating.
#'
#' Current comparison objects record deep-valid compact PSIS source fits,
#' canonical per-source reportability projections, each source authority's
#' validation stamp, and an owned-payload SHA-256 stamp. Derived tables and
#' diagnostics are recomputed from those projections during validation. A
#' pre-marker serialized comparison containing
#' `stack_psis` is inspection-only and must be rebuilt from current validated
#' fits; its saved numeric table is not grandfathered. Independently of the
#' stamp, warning-status PSIS rows are forbidden and blocked PSIS rows must
#' carry no reportable numeric or pooling metadata.
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
#'   `estimate_table` includes method-level interval and PSIS provenance
#'   metadata alongside the aligned fixed-effect estimates. `diagnostic_table`
#'   includes method-level target/pooling provenance, PSIS source and bounded
#'   weight-concentration summaries, and shared-provenance flags.
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
